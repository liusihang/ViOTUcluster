# vOTUcluster Setup and Usage Guide

This guide will help you set up and use the vOTUcluster environment, as well as the vRhyme environment if needed. Follow the steps below to configure and run your analysis.

## Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)
- [mamba](https://github.com/mamba-org/mamba) (for faster package management)

## Setup vOTUcluster Environment

1. **Create and activate the vOTUcluster environment**

    ```bash
    mamba create -n vOTUcluster -c conda-forge -c bioconda "python=3.8" scikit-learn=0.22.1 imbalanced-learn pandas seaborn pyhmmer==0.10.14 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click "conda-package-handling<=1.9"
    mamba activate vOTUcluster
    ```

2. **Download and setup vOTUcluster**

    ```bash
    gzip -d vOTUcluster-master.zip
    cd vOTUcluster-master
    bash setup.sh "/path/to/db" "num"
    ```
    Replace `/path/to/db` with the desired database directory and `num` with the number of threads you wish to use.
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
Database Directory (-d): Path to the directory where the VirSorter database will be set up.
Number of Threads (-j): Number of threads to use for setting up the database and running analyses.
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
