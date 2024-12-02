# ViOTUcluster Setup and Usage Guide（Still working...）

This guide will help you set up and use the ViOTUcluster environment, as well as the vRhyme environment if needed. Follow the steps below to configure and run your analysis.

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)
- [mamba](https://github.com/mamba-org/mamba) (for faster package management)

## Setup ViOTUcluster Environment

1. **Create and activate the ViOTUcluster environment**

    ```bash
    mamba create -n ViOTUcluster -c conda-forge -c bioconda "python=3.8" 
    mamba activate ViOTUcluster
    ```

2. **Download and setup ViOTUcluster**

    ```bash
    git clone https://github.com/liusihang/ViOTUcluster.git
    cd ViOTUcluster-master
    bash setup.sh "/path/to/db" "num"
    ```
    Replace `/path/to/db` with the desired database directory and `num` with the number of threads you wish to use.

   
## Setup vRhyme Environment

If you already have vRhyme environment, you can skip this section.

1. **Create and activate the vRhyme environment**

    ```bash
    mamba env create -f vRhyme.yaml
    mamba activate vRhyme
    ```

2. **Clone the vRhyme repository and install**

    ```bash
    git clone https://github.com/AnantharamanLab/vRhyme
    cd vRhyme
    gunzip vRhyme/models/vRhyme_machine_model_ET.sav.gz
    pip install .
    ```

## Troubleshooting
If you encounter any issues during the setup or usage of these environments, please refer to the following steps:

Ensure that you have the latest version of mamba and conda.
Check that all dependencies are correctly installed.
Verify that the paths and environment variables are set correctly.
For additional help, feel free to open an issue on the respective GitHub repositories.


## How to Use

To run the pipeline, use the following command structure:

1. **Create and activate the vRhyme environment**

    ```bash
    ViOTUcluster -i <input_path_to_contigs> -r <input_path_raw_seqs> -o <output_path> -d <database_path> -n <threads> --non-con/--con [--reassemble]
    ```

2. **Start with raw fastq files**
    ```bash
    ViOTUcluster_AllinOne -r <input_path_raw_seqs> -o <output_path> -d <database_path> -t <sample_type> -a <assembly_software> --non-con/--con [--reassemble]"
    ```
## Parameters

- **`-i <input_path_to_contigs>`**: Specifies the directory containing the assembled contig files in FASTA format (e.g., `example1.fasta`). Each contig file should have corresponding raw sequencing FASTQ files in the raw sequence directory, sharing the same prefix.

- **`-r <input_path_raw_seqs>`**: Specifies the directory with raw sequencing data in FASTQ format. The FASTQ files must have the same prefix as the corresponding contigs file. For example, if the contigs file is `example1.fasta`, the FASTQ files should be named `example1_R1.fq` and `example1_R2.fq`.

- **`-o <output_path>`**: Defines the output directory for storing the processed results. This will include filtered sequences, prediction outcomes, binning results, and the final dereplicated viral contigs.

- **`-d <database_path>`**: Points to the required database for performing viral prediction, binning, and dereplication steps.

- **`--non-con/--con`**: Specifies the viral prediction criteria based on the sample preparation method. Use `--non-con` for samples that were not enriched using viral-particle concentration methods, typically containing a low viral proportion. Use `--con` for samples subjected to concentration methods, which are expected to have a medium to high viral proportion.

- **`--reassemble`**: (Optional) Enables reassembly of bins after the initial binning process to enhance the accuracy and quality of the final contigs. This feature is still in beta and can significantly increase runtime.

- **`-a <assembly_software>`**: (For `ViOTUcluster_AllinOne` only) Specifies the assembly software used during the raw sequence processing. Accepted values are `-a megahit` or `-a metaspades`.
### File Structure Example

Below is a tree list of how the file structure should be organized, assuming the prefix for the example files is `example1`:

```plaintext
<project_directory>/
│
├── input_contigs/
│   ├── example1.fasta
│   └── example2.fasta
│
├── input_fastq/
│   ├── example1_R1.fq
│   ├── example1_R2.fq
│   ├── example2_R1.fq
│   └── example2_R2.fq
│
├── output_path/
│   ├── Summary/
│   │   ├── SeperateRes
│   │   │   ├── example1
│   │   │   │    ├── CheckVRes
│   │   │   │    └── example1_ViralList.fasta
│   │   │   └── example2
│   │   │        ├── CheckVRes
│   │   │        └── example2_ViralList.fasta
│   │   └── vOTU
│   │        ├── vOTU.fasta
│   │        ├── vOTU.Abundance.csv
│   │        ├── vOTU.Taxonomy.csv
│   │        ├── CheckVRes
│   │        ├── DRAMRes(Optional)
│   │        └── iPhopRes(Optional)
│   └── (IntermediateFile....)
│
└── databases/
    ├── db/                # VirSorter2 database
    ├── ViralVerify/       # ViralVerify database
    ├── checkv-db-v1.5/    # CheckV database (version 1.5)
    └── genomad_db/        # Genomad database
```

In this structure:
- `input_contigs/` contains the assembled contigs (e.g., `example1.fasta`).
- `input_fastq/` contains the corresponding FASTQ files (e.g., `example1_R1.fq` and `example1_R2.fq`).
- `output_results/` is the directory where all output files will be stored, organized into various subdirectories such as `FilteredSeqs`, `SeprateFile`, `Summary`, and `CheckRes`.
- `databases/` contains the required databases for the analysis, including:
  - `db/`: The VirSorter2 database.
  - `ViralVerify/`: The ViralVerify database.
  - `checkv-db-v1.5/`: The CheckV database (version 1.5).
  - `genomad_db/`: The Genomad database.

### Final Output

The processed data is organized under the specified `output_path/`, with the following structure:

- **`output_path/Summary`**: Contains the final results and summaries for all processed samples, organized into the following subdirectories:
  - **`SeprateFile`**: Holds individual directories for each sample (e.g., `example1`, `example2`), with subdirectories for specific result types, including:
    - **`CheckVRes`**: Stores quality assessment results from CheckV for each sample.
    - **`<sample>_ViralList.fasta`**: The list of predicted viral contigs for the respective sample.
  - **`vOTU/`**: Contains the final processed viral OTU (vOTU) results across all samples:
    - **`vOTU.fasta`**: The final dereplicated viral contigs after clustering from all samples.
    - **`vOTU.Abundance.csv`**: Abundance data of the vOTUs across samples.
    - **`vOTU.Taxonomy.csv`**: Taxonomic assignments for the vOTUs, if available.
    - **`CheckVRes`**: Summarized CheckV quality assessments for all vOTUs.
    - **`DRAMRes (Optional)`**: Optional results from DRAM annotation if included in the workflow.
    - **`iPhopRes (Optional)`**: Optional results from iPhop annotation if included in the workflow.

- **`output_path/IntermediateFile`**: This directory holds intermediate files generated during the processing pipeline, such as filtered sequences and any temporary data.

- **`databases/`**: Stores the necessary databases used for various stages of the analysis:
  - **`db/`**: The VirSorter2 database.
  - **`ViralVerify/`**: The ViralVerify database, used for viral prediction.
  - **`checkv-db-v1.5/`**: The CheckV database (version 1.5) for quality control of viral sequences.
  - **`genomad_db/`**: The Genomad database for viral identification and dereplication.


## License
This project is licensed under the MIT License - see the LICENSE file for details.
