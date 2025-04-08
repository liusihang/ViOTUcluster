#!/usr/bin/env bash

# Exit on any error
set -e

# Function to check if running in a terminal
is_tty() {
    [ -t 0 ] || [ -t 1 ] || [ -t 2 ]
}

# Redirect echo to stderr if no terminal, otherwise stdout
echo_msg() {
    if is_tty; then
        echo "$@"  # Output to stdout in terminal
    else
        echo "$@" >&2  # Output to stderr in non-TTY
    fi
}

# Get the base installation path of Conda
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo_msg "Error: Conda not found. Please install Conda first."; exit 1; }

# Check if CONDA_BASE is set
if [ -z "$CONDA_BASE" ]; then
    echo_msg "Error: Could not determine Conda base path."
    exit 1
fi

# Source Conda initialization if available
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    . "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo_msg "Error: Conda initialization script not found at $CONDA_BASE/etc/profile.d/conda.sh"
    exit 1
fi

# Create ViOTUcluster environment directory
echo_msg "Creating ViOTUcluster environment directory..."
mkdir -p "$CONDA_BASE/envs/ViOTUcluster" || { echo_msg "Error: Failed to create directory."; exit 1; }

# Download necessary packages
echo_msg "Downloading ViOTUcluster, vRhyme, and DRAM packages..."
wget -q https://zenodo.org/records/15108141/files/ViOTUcluster.tar.gz -O "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" || { echo_msg "Error: Failed to download ViOTUcluster.tar.gz"; exit 1; }
#wget -q https://zenodo.org/records/15108141/files/vRhyme.tar.gz -O "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz" || { echo_msg "Error: Failed to download vRhyme.tar.gz"; exit 1; }
wget -q https://zenodo.org/records/15108141/files/DRAM.tar.gz -O "$CONDA_BASE/envs/ViOTUcluster/DRAM.tar.gz" || { echo_msg "Error: Failed to download DRAM.tar.gz"; exit 1; }

# Extract and unpack ViOTUcluster environment
echo_msg "Extracting and unpacking ViOTUcluster environment..."
tar -xzf "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster" || { echo_msg "Error: Failed to extract ViOTUcluster.tar.gz"; exit 1; }
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || { echo_msg "Warning: conda activate failed, but proceeding with unpack."; }
conda-unpack 2>/dev/null || { echo_msg "Error: Failed to unpack ViOTUcluster environment"; exit 1; }

# Configure SSL certificate verification
echo_msg "Configuring SSL certificate verification..."
conda config --env --set ssl_verify "$CONDA_BASE/envs/ViOTUcluster/ssl/cacert.pem" 2>/dev/null || echo_msg "Warning: Failed to set ssl_verify."
conda env config vars set SSL_CERT_FILE=$(python -c "import certifi; print(certifi.where())" 2>/dev/null) 2>/dev/null || echo_msg "Warning: Failed to set SSL_CERT_FILE."

# Create and unpack vRhyme environment
echo_msg "Creating and unpacking vRhyme environment..."
mkdir -p "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme" || { echo_msg "Error: Failed to create vRhyme directory."; exit 1; }
tar -xzf "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme" || { echo_msg "Error: Failed to extract vRhyme.tar.gz"; exit 1; }
conda activate "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme" 2>/dev/null || { echo_msg "Warning: vRhyme conda activate failed, but proceeding."; }
conda-unpack 2>/dev/null || { echo_msg "Error: Failed to unpack vRhyme environment"; exit 1; }

# Create DRAM environment
echo_msg "Creating and unpacking DRAM environment..."
mkdir -p "$CONDA_BASE/envs/ViOTUcluster/envs/DRAM" || { echo_msg "Error: Failed to create DRAM directory."; exit 1; }

wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml || { echo_msg "Error: Failed to download environment.yaml"; exit 1; }
mamba env create -f environment.yaml -p "$CONDA_BASE/envs/ViOTUcluster/envs/DRAM" || { echo_msg "Error: Failed to create DRAM environment"; exit 1; }
conda run -p "$CONDA_BASE/envs/ViOTUcluster/envs/DRAM" bash -c "echo 'DRAM environment activated'" || { echo_msg "Warning: DRAM conda activate failed, but proceeding."; }

echo_msg "DRAM environment setup completed."

# Create iPhop environment
echo_msg "Creating iPhop environment..."

# Try to activate ViOTUcluster environment
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || {
    echo_msg "Warning: Failed to activate ViOTUcluster, proceeding."
}
# Determine whether to use mamba or conda
if command -v mamba >/dev/null 2>&1; then
    mamba_cmd="mamba"
else
    echo_msg "Warning: Mamba not found, falling back to conda"
    mamba_cmd="conda"
fi
# Create iPhop environment
"$mamba_cmd" create -y -c conda-forge -p "$CONDA_BASE/envs/ViOTUcluster/envs/iPhop" python=3.8 mamba || {
    echo_msg "Error: Failed to create iPhop environment"
    exit 1
}
# Activate iPhop environment
conda activate "$CONDA_BASE/envs/ViOTUcluster/envs/iPhop" 2>/dev/null || {
    echo_msg "Warning: iPhop conda activate failed, proceeding."
}
# Install iPhop
"$mamba_cmd" install -y -c conda-forge -c bioconda iphop || {
    echo_msg "Error: Failed to install iphop"
    exit 1
}

# Final activation of ViOTUcluster
echo_msg "Finalizing setup..."
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || { echo_msg "Warning: Final conda activate failed, environment should still be usable."; }

echo_msg "[✅] ViOTUcluster Setup complete."
echo_msg "Current version: 0.5.2"

# Clean up downloaded files
rm -f "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" 2>/dev/null
rm -f "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz" 2>/dev/null
rm -f "$CONDA_BASE/envs/ViOTUcluster/DRAM.tar.gz" 2>/dev/null
