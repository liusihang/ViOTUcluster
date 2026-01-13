#!/usr/bin/env python3
"""
Gene catalog workflow:
  - Prodigal gene prediction per sample
  - Filter genes (exclude partial=11)
  - Salmon quantification per sample
  - Merge TPM / counts
  - mmseqs2 clustering on protein sequences
"""

from __future__ import annotations

import argparse
import csv
import glob
import logging
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .config import FASTA_EXTENSIONS, FASTQ_EXTENSIONS


LOGGER = logging.getLogger(__name__)


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _check_tools(required: Sequence[str]) -> None:
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(f"Missing required tools: {', '.join(missing)}")


def _find_fasta_files(directory: str) -> List[str]:
    files: List[str] = []
    for ext in FASTA_EXTENSIONS:
        files.extend(glob.glob(os.path.join(directory, f"*{ext}")))
    return sorted(set(files))


def _discover_contigs(output_dir: str, contigs_dir: Optional[str]) -> List[str]:
    if contigs_dir:
        return _find_fasta_files(contigs_dir)

    filtered_dir = os.path.join(output_dir, "FilteredSeqs")
    if os.path.isdir(filtered_dir):
        filtered = _find_fasta_files(filtered_dir)
        if filtered:
            return filtered

    return []


def _sample_name_from_path(path: str) -> str:
    base = os.path.basename(path)
    for ext in FASTA_EXTENSIONS:
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
            break
    for suffix in ("_viralseqs", "_unbined", "_bins"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return base or os.path.basename(path)


def _find_read_pair(raw_seq_dir: str, sample_name: str) -> Optional[Tuple[str, str]]:
    for ext in FASTQ_EXTENSIONS:
        r1 = os.path.join(raw_seq_dir, f"{sample_name}_R1{ext}")
        r2 = os.path.join(raw_seq_dir, f"{sample_name}_R2{ext}")
        if os.path.exists(r1) and os.path.exists(r2):
            return r1, r2
    for ext in FASTQ_EXTENSIONS:
        r1 = os.path.join(raw_seq_dir, f"{sample_name}_1{ext}")
        r2 = os.path.join(raw_seq_dir, f"{sample_name}_2{ext}")
        if os.path.exists(r1) and os.path.exists(r2):
            return r1, r2
    return None


def _parse_gene_ids_from_faa(faa_path: str) -> List[str]:
    """
    Collect gene IDs from the protein FASTA headers, excluding partial=11.

    Prodigal headers typically look like:
      >contig_1_1 # 132 # 293 # -1 # ID=1_1;partial=01;...
    """
    gene_ids: List[str] = []
    with open(faa_path, "r", encoding="utf-8") as handle:
        for header, _seq in SimpleFastaParser(handle):
            header = header.strip()
            # Attribute block is after the last '#'
            attr_block = ""
            if "#" in header:
                attr_block = header.split("#")[-1].strip()
            attr_dict = {}
            if attr_block:
                for item in attr_block.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attr_dict[key.strip()] = value.strip()
            if attr_dict.get("partial") == "11":
                continue
            gene_id = header.split()[0]
            if gene_id:
                gene_ids.append(gene_id)
    return gene_ids


def _filter_fasta_by_ids(input_fasta: str, output_fasta: str, target_ids: Iterable[str]) -> int:
    target_set = set(target_ids)
    count = 0
    with open(output_fasta, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in target_set:
                SeqIO.write(record, out_handle, "fasta")
                count += 1
    return count


def _run_prodigal(contig_fasta: str, gene_pred_dir: str) -> Tuple[str, str, str]:
    _ensure_dir(gene_pred_dir)
    gff_path = os.path.join(gene_pred_dir, "genes.gff")
    fna_path = os.path.join(gene_pred_dir, "genes.fna")
    faa_path = os.path.join(gene_pred_dir, "genes.faa")

    if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in (gff_path, fna_path, faa_path)):
        return gff_path, fna_path, faa_path

    cmd = [
        "prodigal",
        "-p",
        "meta",
        "-i",
        contig_fasta,
        "-o",
        gff_path,
        "-a",
        faa_path,
        "-d",
        fna_path,
        "-q",
    ]
    subprocess.run(cmd, check=True)
    return gff_path, fna_path, faa_path


def _run_salmon(
    complete_fna: str,
    read1: str,
    read2: str,
    threads: int,
    abundance_dir: str,
) -> str:
    _ensure_dir(abundance_dir)
    index_dir = os.path.join(abundance_dir, "salmon_index")
    quant_dir = os.path.join(abundance_dir, "quant_result")
    quant_file = os.path.join(quant_dir, "quant.sf")

    if not os.path.exists(complete_fna) or os.path.getsize(complete_fna) == 0:
        raise RuntimeError(f"No complete gene sequences found at {complete_fna}")

    if not os.path.isdir(index_dir):
        cmd_index = [
            "salmon",
            "index",
            "-t",
            complete_fna,
            "-i",
            index_dir,
            "-p",
            str(threads),
        ]
        subprocess.run(cmd_index, check=True)

    if not os.path.isfile(quant_file):
        cmd_quant = [
            "salmon",
            "quant",
            "-i",
            index_dir,
            "-l",
            "A",
            "-1",
            read1,
            "-2",
            read2,
            "-p",
            str(threads),
            "--meta",
            "--validateMappings",
            "-o",
            quant_dir,
        ]
        subprocess.run(cmd_quant, check=True)

    return quant_file


def _merge_quant_files(
    quant_files: Sequence[str],
    output_tpm: str,
    output_counts: str,
    add_sample_prefix: bool,
) -> None:
    all_tpm = []
    all_counts = []

    for qfile in quant_files:
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(qfile)))
        df = pd.read_csv(qfile, sep="\t")
        df = df[["Name", "TPM", "NumReads"]]

        if add_sample_prefix:
            df["Name"] = f"{sample_name}|" + df["Name"].astype(str)

        series_tpm = df.set_index("Name")["TPM"]
        series_tpm.name = sample_name
        all_tpm.append(series_tpm)

        series_counts = df.set_index("Name")["NumReads"]
        series_counts.name = sample_name
        all_counts.append(series_counts)

    merged_tpm = pd.concat(all_tpm, axis=1, sort=False).fillna(0)
    merged_tpm.index.name = "GeneID"
    merged_counts = pd.concat(all_counts, axis=1, sort=False).fillna(0)
    merged_counts.index.name = "GeneID"

    merged_tpm.to_csv(output_tpm)
    merged_counts.to_csv(output_counts)


def _merge_proteins_with_prefix(
    faa_files: Sequence[Tuple[str, str]],
    merged_faa: str,
) -> None:
    with open(merged_faa, "w", encoding="utf-8") as out_handle:
        for sample_name, faa_path in faa_files:
            with open(faa_path, "r", encoding="utf-8") as in_handle:
                for header, seq in SimpleFastaParser(in_handle):
                    original_id = header.split()[0]
                    out_handle.write(f">{sample_name}_{original_id}\n{seq}\n")


def _write_mmseqs_csvs(
    rep_fasta: str,
    cluster_tsv: str,
    out_seq_csv: str,
    out_cluster_csv: str,
) -> None:
    if os.path.exists(rep_fasta):
        with open(out_seq_csv, "w", newline="", encoding="utf-8") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["prot_id", "seq"])
            with open(rep_fasta, "r", encoding="utf-8") as f_in:
                for header, seq in SimpleFastaParser(f_in):
                    clean_id = header.split()[0]
                    writer.writerow([clean_id, seq])

    if os.path.exists(cluster_tsv):
        with open(out_cluster_csv, "w", newline="", encoding="utf-8") as f_out:
            writer = csv.writer(f_out)
            writer.writerow(["representative_id", "member_id"])
            with open(cluster_tsv, "r", encoding="utf-8") as f_in:
                tsv_reader = csv.reader(f_in, delimiter="\t")
                for row in tsv_reader:
                    if row:
                        writer.writerow(row)


def _run_mmseqs_cluster(
    faa_files: Sequence[Tuple[str, str]],
    mmseqs_dir: str,
    min_id: float,
    cov: float,
) -> None:
    _ensure_dir(mmseqs_dir)
    merged_faa = os.path.join(mmseqs_dir, "all_samples_merged.faa")
    mmseqs_prefix = os.path.join(mmseqs_dir, "mmseqs_result")
    tmp_dir = os.path.join(mmseqs_dir, "tmp")
    rep_fasta = f"{mmseqs_prefix}_rep_seq.fasta"
    cluster_tsv = f"{mmseqs_prefix}_cluster.tsv"
    out_seq_csv = os.path.join(mmseqs_dir, "deduplicated_seqs.csv")
    out_cluster_csv = os.path.join(mmseqs_dir, "cluster_info.csv")

    if not (os.path.exists(rep_fasta) and os.path.exists(cluster_tsv)):
        _merge_proteins_with_prefix(faa_files, merged_faa)
        cmd = [
            "mmseqs",
            "easy-cluster",
            merged_faa,
            mmseqs_prefix,
            tmp_dir,
            "--min-seq-id",
            str(min_id),
            "-c",
            str(cov),
            "--cov-mode",
            "1",
        ]
        subprocess.run(cmd, check=True)

    _write_mmseqs_csvs(rep_fasta, cluster_tsv, out_seq_csv, out_cluster_csv)


def _process_sample(
    sample_name: str,
    contig_fasta: str,
    raw_seq_dir: str,
    output_root: str,
    salmon_threads: int,
) -> str:
    sample_dir = os.path.join(output_root, "Samples", sample_name)
    gene_pred_dir = os.path.join(sample_dir, "03_GenePred")
    abundance_dir = os.path.join(sample_dir, "04_Abundance")
    _ensure_dir(gene_pred_dir)
    _ensure_dir(abundance_dir)

    read_pair = _find_read_pair(raw_seq_dir, sample_name)
    if not read_pair:
        return f"{sample_name}: skipped (paired reads not found)"

    gff_path, fna_path, faa_path = _run_prodigal(contig_fasta, gene_pred_dir)
    gene_ids = _parse_gene_ids_from_faa(faa_path)
    if not gene_ids:
        return f"{sample_name}: skipped (no genes after filtering partial=11)"

    complete_fna = os.path.join(abundance_dir, "complete_genes.fna")
    complete_faa = os.path.join(abundance_dir, "complete_genes.faa")

    if not os.path.exists(complete_fna):
        _filter_fasta_by_ids(fna_path, complete_fna, gene_ids)
    if not os.path.exists(complete_faa):
        _filter_fasta_by_ids(faa_path, complete_faa, gene_ids)

    _run_salmon(complete_fna, read_pair[0], read_pair[1], salmon_threads, abundance_dir)
    return f"{sample_name}: done"


def run_gene_catalog(
    output_dir: str,
    raw_seq_dir: str,
    contigs_dir: Optional[str] = None,
    max_parallel: int = 1,
    salmon_threads: int = 1,
    mmseqs_min_id: float = 0.9,
    mmseqs_cov: float = 0.8,
    add_sample_prefix: bool = True,
) -> bool:
    _check_tools(["prodigal", "salmon", "mmseqs"])

    contig_files = _discover_contigs(output_dir, contigs_dir)
    if not contig_files:
        LOGGER.error("No contigs found for gene catalog.")
        return False

    gene_root = os.path.join(output_dir, "Summary", "GeneCatalog")
    _ensure_dir(gene_root)
    _ensure_dir(os.path.join(gene_root, "Samples"))

    if max_parallel < 1:
        max_parallel = 1
    if salmon_threads < 1:
        salmon_threads = 1

    LOGGER.info("Gene catalog inputs: %d contig files", len(contig_files))
    LOGGER.info("Gene catalog output: %s", gene_root)

    errors = []
    if max_parallel == 1:
        for contig_fasta in contig_files:
            sample_name = _sample_name_from_path(contig_fasta)
            try:
                msg = _process_sample(sample_name, contig_fasta, raw_seq_dir, gene_root, salmon_threads)
                LOGGER.info(msg)
            except Exception as exc:
                errors.append(f"{sample_name}: {exc}")
    else:
        with ProcessPoolExecutor(max_workers=max_parallel) as executor:
            future_map = {}
            for contig_fasta in contig_files:
                sample_name = _sample_name_from_path(contig_fasta)
                future = executor.submit(
                    _process_sample,
                    sample_name,
                    contig_fasta,
                    raw_seq_dir,
                    gene_root,
                    salmon_threads,
                )
                future_map[future] = sample_name
            for future in as_completed(future_map):
                sample_name = future_map[future]
                try:
                    msg = future.result()
                    LOGGER.info(msg)
                except Exception as exc:
                    errors.append(f"{sample_name}: {exc}")

    if errors:
        for err in errors:
            LOGGER.error("Gene catalog error: %s", err)
        return False

    quant_files = glob.glob(
        os.path.join(gene_root, "Samples", "*", "04_Abundance", "quant_result", "quant.sf")
    )
    if quant_files:
        output_tpm = os.path.join(gene_root, "Gene.Abundance.TPM.csv")
        output_counts = os.path.join(gene_root, "Gene.Abundance.Counts.csv")
        _merge_quant_files(quant_files, output_tpm, output_counts, add_sample_prefix)
        LOGGER.info("Merged abundance tables written to %s", gene_root)
    else:
        LOGGER.warning("No salmon quant files found; skipping merge step.")

    faa_files = []
    for faa_path in glob.glob(
        os.path.join(gene_root, "Samples", "*", "04_Abundance", "complete_genes.faa")
    ):
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(faa_path)))
        faa_files.append((sample_name, faa_path))

    if faa_files:
        mmseqs_dir = os.path.join(gene_root, "mmseqs")
        _run_mmseqs_cluster(faa_files, mmseqs_dir, mmseqs_min_id, mmseqs_cov)
        LOGGER.info("mmseqs2 outputs written to %s", mmseqs_dir)
    else:
        LOGGER.warning("No protein FASTA found for mmseqs2 clustering.")

    return True


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="ViOTUcluster gene catalog workflow.")
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="ViOTUcluster output directory (used for Summary/GeneCatalog).",
    )
    parser.add_argument(
        "-r",
        "--raw-seqs",
        required=True,
        help="Directory with raw paired-end reads.",
    )
    parser.add_argument(
        "--contigs-dir",
        default=None,
        help="Optional contigs directory (default: output/Summary/SeperateRes).",
    )
    parser.add_argument(
        "--max-parallel",
        type=int,
        default=1,
        help="Max parallel samples for Prodigal/Salmon (default: 1).",
    )
    parser.add_argument(
        "--salmon-threads",
        type=int,
        default=1,
        help="Threads per Salmon task (default: 1).",
    )
    parser.add_argument(
        "--mmseqs-min-id",
        type=float,
        default=0.9,
        help="mmseqs2 min sequence identity (default: 0.9).",
    )
    parser.add_argument(
        "--mmseqs-cov",
        type=float,
        default=0.8,
        help="mmseqs2 coverage threshold (default: 0.8).",
    )
    parser.add_argument(
        "--no-sample-prefix",
        action="store_true",
        help="Do not prefix gene IDs with sample name when merging abundance.",
    )

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    ok = run_gene_catalog(
        output_dir=os.path.abspath(args.output_dir),
        raw_seq_dir=os.path.abspath(args.raw_seqs),
        contigs_dir=os.path.abspath(args.contigs_dir) if args.contigs_dir else None,
        max_parallel=args.max_parallel,
        salmon_threads=args.salmon_threads,
        mmseqs_min_id=args.mmseqs_min_id,
        mmseqs_cov=args.mmseqs_cov,
        add_sample_prefix=not args.no_sample_prefix,
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
