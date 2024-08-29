# vOTUcluster Setup and Usage Guide

This guide will help you set up and use the vOTUcluster environment, as well as the vRhyme environment if needed. Follow the steps below to configure and run your analysis.

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)
- [mamba](https://github.com/mamba-org/mamba) (for faster package management)

## Setup vOTUcluster Environment

1. **Create and activate the vOTUcluster environment**

    ```bash
    mamba create -n vOTUcluster -c conda-forge -c bioconda "python=3.8" 
    mamba activate vOTUcluster
    ```

2. **Download and setup vOTUcluster**

    ```bash
    git clone https://github.com/liusihang/vOTUcluster.git
    cd vOTUcluster-master
    bash setup.sh "/path/to/db" "num"
    ```
    Replace `/path/to/db` with the desired database directory and `num` with the number of threads you wish to use.

3. **Download viralverify database**

   Download viralverify database from https://figshare.com/s/f897d463b31a35ad7bf0.

   Unzip it to `/path/to/db/viralverify`

   The finial folder for viralverify should looks like `/path/to/db/viralverify/nbc_hmms.hmm`
   
## Setup vRhyme Environment

If you already have vRhyme, you can skip this section.

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

```bash
mamba activate vOTUcluster
```


## Parameters
- **`-d <database_path>`**: Path to the directory where the VirSorter database will be set up.

- **`-j <Number of Threads>`**: Number of threads to use for setting up the database and running analyses.
## Troubleshooting
If you encounter any issues during the setup or usage of these environments, please refer to the following steps:

Ensure that you have the latest version of mamba and conda.
Check that all dependencies are correctly installed.
Verify that the paths and environment variables are set correctly.
For additional help, feel free to open an issue on the respective GitHub repositories.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements
Special thanks to the developers of the following tools and libraries used in this project:

mamba
conda-forge
bioconda


## How to Use

To run the pipeline, use the following command structure:

```markdown
vOTUcluster -i <input_path_to_contigs> -r <input_path_fastq> -o <output_path> -d <database_path> -t <sample_type> [--reassemble]
```

## Parameters

- **`-i <input_path_to_contigs>`**: Specifies the input directory containing the assembled contigs files. These files should be in FASTA format (e.g., `example1.fasta`). Each contig file should have corresponding FASTQ files with the same prefix in the raw sequence directory.

- **`-r <input_path_fastq>`**: Specifies the input directory containing the raw sequencing files (FASTQ format). The FASTQ files should have the same prefix as the corresponding contigs file. For instance, if the contigs file is named `example1.fasta`, the FASTQ files should be named `example1_R1.fq` and `example1_R2.fq`.

- **`-o <output_path>`**: Specifies the output directory where the processed results will be stored. This directory will include filtered sequences, prediction results, binning outcomes, and final dereplicated viral contigs.

- **`-d <database_path>`**: Specifies the path to the database required for various analysis steps, such as viral prediction, binning, and dereplication.

- **`-t <sample_type>`**: Specifies the type of sample being analyzed. Possible values are `DNA`, `RNA`, or `Mix`, which will determine the viral groups used during the prediction step.

- **`--reassemble`**: (Optional) If included, this option enables the reassembly of bins after the initial binning process. This step can improve the accuracy and quality of the final contigs. Notably, this function is still in beta stage. Enable `--reassemble` will significantly increase running time.

### File Structure Example

Below is a tree list of how the file structure should be organized, assuming the prefix for the example files is `example1`:

```plaintext
<project_directory>/
│
├── input_contigs/
│   └── example1.fasta
│
├── input_fastq/
│   ├── example1_R1.fq
│   └── example1_R2.fq
│
├── output_results/
│   ├── FilteredSeqs/
│   ├── SeprateFile/
│   ├── Summary/
│   │   ├── unbined/
│   │   ├── bins/
│   │   └── Viralcontigs/
│   └── CheckRes/
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

The processed data is organized in the following directory structure under the specified `OUTPUT_DIR`:

- **`OUTPUT_DIR/FilteredSeqs`**: Contains the filtered sequences that are longer than 2000bp.
- **`OUTPUT_DIR/SeprateFile`**: Holds individual directories for each sample, including:
  - **`RoughViralPrediction`**: Contains the results from ViralVerify, VirSorter2, and Genomad predictions.
  - **`Binning`**: Holds the results from vRhyme binning, including reassembly (if enabled).
  - **`Summary`**: Contains the summarized results from binning, reassembly, and cross-validation.
  
- **`OUTPUT_DIR/Summary`**: Contains the final summary of all processed sequences, organized as follows:
  - **`bins`**: Contains dereplicated bins from all samples.
  - **`Viralcontigs`**: Final viral contigs (`vOTU.fasta`) combining both bins and unbinned sequences after clustering.
  - **`CheckRes`**: Quality summary results from CheckV, assessing the quality of the final viral contigs.

