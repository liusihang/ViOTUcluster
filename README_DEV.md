# ViOTUcluster Developer Quickread

Fast orientation for maintainers.

## What it is
- Python + Bash hybrid pipeline for viral OTU processing from metagenomic reads.
- Entrypoints: `ViOTUcluster` (contigs + raw reads), `ViOTUcluster_AllinOne` (raw reads → assembly → full pipeline), optional gene catalog sidecar.

## Primary flows
1) **Filtering** (`filter_contigs.py`): min-length gate, outputs `FilteredSeqs/`.
2) **Viral prediction** (`viral_prediction_module.sh` → `viralprediction.py`): viralverify / virsorter2 / genomad in parallel.
3) **Cross validation + CheckV** (`cross_validation_module.sh`): `CrossValid.py`, `FilterRawResSeqs.py`, `check_removal.py`.
4) **Binning** (`binning_merge_module.sh`): vRhyme on per-sample contigs, optional reassembly, outputs bins + unbinned FASTA.
5) **dRep** (`drep_module.sh`): deduplicate across samples.
6) **Summary** (`summary_module.sh`): CheckV on merged vOTUs, coverm TPM via BWA+sambamba, taxonomy via genomad.
7) **Advanced (optional)**: DRAM + iPhop.

## Gene catalog (optional or standalone)
- Module: `ViOTUcluster/gene_catalog.py`.
- CLI flags: `--gene-catalog` (after summary) or `--gene-catalog-only` (skip main flow).
- Steps: Prodigal → filter complete genes → Salmon quant → merge TPM/counts → mmseqs2 clustering.
- Outputs: `Summary/GeneCatalog/` with TPM/Counts matrices and mmseqs CSVs.

## Key files
- Python: `ViOTUcluster/pipeline.py` (orchestration), `viotucluster_cli.py`, `viotucluster_allinone_cli.py`, `config.py`.
- Bash modules: `Modules/*.sh` (prediction, binning, summary, etc.), wrappers `Modules/ViOTUcluster`, `Modules/ViOTUcluster_AllinOne`.
- Setup: `setup.py` exposes console scripts, ships modules.

## Running (minimal)
- Contigs + reads: `ViOTUcluster -i <contigs_dir> -r <reads_dir> -o <out> -d <db> --con|--non-con [tuning flags] [--gene-catalog]`.
- All-in-one: `ViOTUcluster_AllinOne -r <raw_reads> -o <out> -d <db> -a megahit|metaspades --con|--non-con [tuning flags]`.
- Gene catalog only: `ViOTUcluster -r <reads_dir> -o <out> --gene-catalog-only [--gene-contigs-dir <dir>]`.

## Concurrency knobs
- Prediction: `--max-prediction-tasks` (VirSorter2/ViralVerify/genomad).
- TPM/BAM: `--tpm-tasks`.
- Assembly (all-in-one): `--assemble-jobs`.
- Threads per heavy task: `-n/--threads`.

## Outputs (core)
- `Summary/SeperateRes/`: per-sample viral FASTA, bins/unbinned.
- `Summary/vOTU/`: merged vOTU FASTA, abundance, taxonomy, CheckV.
- `Summary/GeneCatalog/`: gene abundance matrices + mmseqs clustering (if enabled).

## Databases expected
- Under `DATABASE`: `db` (VirSorter2), `ViralVerify/`, `checkv-db-v1.5/`, `genomad_db/`, plus `DRAM/`, `Aug_2023_pub_rw/` for advanced steps.

## Notes
- Bash wrappers set envs for modules; Python pipeline uses `.Modules` lookup via `ScriptDir`.
- Keep ASCII; avoid altering user-created files in `tempscript/`.
