#!/usr/bin/env bash

#source activate DRAM

# Check command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta_file> <output_directory>"
    exit 1
fi

# Preserve pipeline-level OUTPUT_DIR if it exists (used later for reporting)
PIPELINE_OUTPUT_DIR="${OUTPUT_DIR:-}"

# Accept command-line arguments
INPUT_FASTA=$1
OUTPUT_DIR=$2
if [ -z "$PIPELINE_OUTPUT_DIR" ]; then
    PIPELINE_OUTPUT_DIR="$OUTPUT_DIR"
fi
mkdir -p "${PIPELINE_OUTPUT_DIR}/Log"

# Check if output file already exists
if [ -f "$OUTPUT_DIR/DRAM_annotations.tsv" ]; then
    echo "Output file '$OUTPUT_DIR/DRAM_annotations.tsv' already exists. Skipping analysis."
    exit 0
fi

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/split_files"

echo -e "\n\n\n# Performing DRAM analysis!!!\n\n\n"
pwd

# Get the base filename of the input file (without path)
BASE_INPUT_FASTA=$(basename "$INPUT_FASTA")

# Split the input fasta file into smaller files, each containing 1000 sequences
awk -v output_dir="$OUTPUT_DIR/split_files" -v base_name="$BASE_INPUT_FASTA" 'BEGIN {n_seq=0;} 
     /^>/ {
        if (n_seq % 1000 == 0) {
            file = sprintf("%s/%s_%d.fna", output_dir, base_name, n_seq);
        } 
        print >> file; 
        n_seq++; 
        next;
     } 
     { print >> file; }' "$INPUT_FASTA"

# Change to the split files directory
cd "$OUTPUT_DIR/split_files" || exit

# List all split fna files
ls *.fna > DRAM

# Run the Python script for DRAM annotation
python "${ScriptDir}/run_DRAM.py"

all_tasks_completed=false

# Monitor task completion
while [ "$all_tasks_completed" == "false" ]; do
    sleep 30
    all_tasks_completed=true

    # Iterate over all directories ending with _DRAMAnnot
    for dir in *_DRAMAnnot; do
        if [ ! -f "$dir/annotations.tsv" ]; then
            echo "DRAM annotation still in progress in $dir."
            all_tasks_completed=false
            break
        fi
    done

    # If not completed, wait another 30 seconds
    if [ "$all_tasks_completed" == "false" ]; then
        sleep 30
    fi
done

echo "All DRAM annotations completed."

# Merge all annotation results into the output directory
awk 'FNR==1 && NR!=1{next;} {print}' "$OUTPUT_DIR"/split_files/*_DRAMAnnot/annotations.tsv > "$OUTPUT_DIR/DRAM_annotations.tsv"

echo "Annotation complete. Results combined and saved to $OUTPUT_DIR/DRAM_annotations.tsv"

# Aggregate gene FASTA sequences for downstream abundance calculations
GENE_FNA="$OUTPUT_DIR/DRAM_genes.fna"
if [ -f "$GENE_FNA" ]; then
    echo "Aggregated gene FASTA already exists at $GENE_FNA. Skipping concatenation."
else
    shopt -s nullglob
    gene_sources=("$OUTPUT_DIR"/split_files/*_DRAMAnnot/genes.fna)
    shopt -u nullglob
    if [ "${#gene_sources[@]}" -eq 0 ]; then
        echo "Warning: No genes.fna files were produced by DRAM; gene abundance calculation will be skipped."
        GENE_FNA=""
    else
        cat "${gene_sources[@]}" > "$GENE_FNA"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to aggregate DRAM gene FASTA files."
            exit 1
        fi
        echo "Aggregated $(printf '%s\n' "${gene_sources[@]}" | wc -l | tr -d ' ') gene FASTA files into $GENE_FNA"
    fi
fi

# Gene abundance calculation (mirrors vOTU TPM workflow)
if [ -n "$GENE_FNA" ] && [ -s "$GENE_FNA" ]; then
    if [ -z "${RAW_SEQ_DIR:-}" ] || [ -z "${FILES:-}" ]; then
        echo "Warning: RAW_SEQ_DIR or FILES environment variables are missing; skipping gene abundance calculation."
    else
        GENE_ABUNDANCE_CSV="$OUTPUT_DIR/DRAM_Gene.Abundance.csv"
        if [ -f "$GENE_ABUNDANCE_CSV" ]; then
            echo "Gene abundance table already exists at $GENE_ABUNDANCE_CSV. Skipping recalculation."
        else
            echo "Starting gene abundance calculation using DRAM gene predictions..."
            GENE_TEMP_DIR="$OUTPUT_DIR/GeneTPMTemp"
            GENE_COVERAGE_DIR="$OUTPUT_DIR/GeneCoverage"
            mkdir -p "$GENE_TEMP_DIR" "$GENE_COVERAGE_DIR"
            if [ $? -ne 0 ]; then
                echo "Error: Failed to create temporary directories for gene abundance calculation."
                exit 1
            fi

            GENE_INDEX_PREFIX="$GENE_TEMP_DIR/GeneIndex"
            if [ ! -f "${GENE_INDEX_PREFIX}.sa" ]; then
                echo "[🔄] Building BWA index for DRAM genes..."
                bwa index -b "100000000" -p "$GENE_INDEX_PREFIX" "$GENE_FNA"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to build BWA index for DRAM genes."
                    exit 1
                fi
            else
                echo "[✅] BWA index for DRAM genes already exists. Skipping index build."
            fi

            for FILE in $FILES; do
                BASENAME=$(basename "$FILE" .fa)
                BASENAME=${BASENAME%.fasta}
                if [ -z "$BASENAME" ]; then
                    echo "Warning: Unable to derive sample name from $FILE. Skipping."
                    continue
                fi

                COVERAGE_FILE="$GENE_COVERAGE_DIR/${BASENAME}_gene_coverage.tsv"
                if [ -f "$COVERAGE_FILE" ]; then
                    echo "[✅] Gene coverage for ${BASENAME} already exists. Skipping."
                    continue
                fi

                Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
                Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
                if [ -z "$Read1" ] || [ -z "$Read2" ]; then
                    echo "Warning: Missing paired reads for ${BASENAME}. Skipping gene abundance for this sample."
                    continue
                fi

                SORTED_BAM="$GENE_TEMP_DIR/${BASENAME}_sorted_gene.bam"
                if [ ! -f "$SORTED_BAM" ]; then
                    echo "[🔄] Aligning reads for ${BASENAME} against DRAM genes..."
                    bwa mem -t "${THREADS}" "$GENE_INDEX_PREFIX" "${Read1}" "${Read2}" 2>> "${PIPELINE_OUTPUT_DIR}/Log/DRAM_gene_abundance.log" | \
                        sambamba view -S -f bam -t "${THREADS}" /dev/stdin | \
                        sambamba sort -t "${THREADS}" -o "$SORTED_BAM" /dev/stdin
                    if [ $? -ne 0 ]; then
                        echo "Error: Alignment pipeline failed for ${BASENAME} during gene abundance calculation."
                        exit 1
                    fi
                    sambamba index -t "${THREADS}" "$SORTED_BAM"
                    if [ $? -ne 0 ]; then
                        echo "Error: Failed to index BAM for ${BASENAME} during gene abundance calculation."
                        exit 1
                    fi
                else
                    echo "[✅] Alignment already completed for ${BASENAME}. Skipping alignment."
                fi

                echo "[🔄] Calculating gene TPM for ${BASENAME}..."
                coverm contig --methods tpm --bam-files "$SORTED_BAM" -t "${THREADS}" -o "$COVERAGE_FILE"
                if [ $? -ne 0 ]; then
                    echo "Error: coverm contig failed for ${BASENAME} during gene abundance calculation."
                    exit 1
                fi
            done

            if [ -d "$GENE_COVERAGE_DIR" ] && [ "$(ls -1 "$GENE_COVERAGE_DIR"/*.tsv 2>/dev/null | wc -l)" -gt 0 ]; then
                python ${ScriptDir}/TPM_caculate.py "$GENE_COVERAGE_DIR" "$GENE_ABUNDANCE_CSV" Gene
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to merge gene TPM tables."
                    exit 1
                fi
                echo "[✅] Gene abundance table written to $GENE_ABUNDANCE_CSV"
            else
                echo "Warning: No gene coverage files were generated; skipping merge step."
            fi
        fi
    fi
fi

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/split_files"
rm -rf "$OUTPUT_DIR/DRAM_results"/*_DRAMAnnot

echo "Cleanup complete."
#conda deactivate
