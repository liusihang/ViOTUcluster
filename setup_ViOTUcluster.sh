#!/usr/bin/env bash

# Exit on any error
set -e

# --- Configuration for Download URLs ---
DEFAULT_VIOTUCLUSTER_URL="https://zenodo.org/records/15108141/files/ViOTUcluster.tar.gz"
DEFAULT_VRHYME_URL="https://zenodo.org/records/15108141/files/vRhyme.tar.gz"

CHINA_VIOTUCLUSTER_URL="https://china.scidb.cn/download?fileId=b7fafcfff327e317438e35c4e05a0050&traceId=d641a200-44c5-4f88-beab-e3c92322a486"
CHINA_VRHYME_URL="https://china.scidb.cn/download?fileId=97d5f6f255685b037c1b196799d4351d&traceId=d641a200-44c5-4f88-beab-e3c92322a486"

VIOTUCLUSTER_URL="$DEFAULT_VIOTUCLUSTER_URL"
VRHYME_URL="$DEFAULT_VRHYME_URL"
# --- End Configuration ---

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

# --- Option Parsing for China Mirror ---
if [[ "$1" == "china" || "$1" == "--china" ]]; then
    echo_msg "Using China download mirror (china.scidb.cn)."
    VIOTUCLUSTER_URL="$CHINA_VIOTUCLUSTER_URL"
    VRHYME_URL="$CHINA_VRHYME_URL"
    # Shift argument if you plan to parse more options later
    # shift
elif [[ -n "$1" && "$1" != "china" && "$1" != "--china" ]]; then
    echo_msg "Usage: $0 [china|--china]"
    echo_msg "  Pass 'china' or '--china' to use download mirrors in China."
    echo_msg "  Ignoring unrecognized argument: $1"
    echo_msg "Proceeding with default download URLs from Zenodo."
elif [[ -z "$1" ]]; then
    echo_msg "Proceeding with default download URLs from Zenodo."
    echo_msg "Tip: You can use '$0 china' or '$0 --china' to use download mirrors in China."
fi
# --- End Option Parsing ---


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
echo_msg "Downloading ViOTUcluster from: $VIOTUCLUSTER_URL"
wget -q "$VIOTUCLUSTER_URL" -O "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" || { echo_msg "Error: Failed to download ViOTUcluster.tar.gz"; exit 1; }

echo_msg "Downloading vRhyme from: $VRHYME_URL"
wget -q "$VRHYME_URL" -O "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz" || { echo_msg "Error: Failed to download vRhyme.tar.gz"; exit 1; }
#wget -q https://zenodo.org/records/15108141/files/DRAM.tar.gz -O "$CONDA_BASE/envs/ViOTUcluster/DRAM.tar.gz" || { echo_msg "Error: Failed to download DRAM.tar.gz"; exit 1; }

# Extract and unpack ViOTUcluster environment
echo_msg "Extracting and unpacking ViOTUcluster environment..."
tar -xzf "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster" || { echo_msg "Error: Failed to extract ViOTUcluster.tar.gz"; exit 1; }
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || { echo_msg "Warning: conda activate failed for ViOTUcluster, but proceeding with unpack."; }
conda-unpack 2>/dev/null || { echo_msg "Error: Failed to unpack ViOTUcluster environment"; exit 1; }

# Configure SSL certificate verification
echo_msg "Configuring SSL certificate verification..."
conda config --env --set ssl_verify "$CONDA_BASE/envs/ViOTUcluster/ssl/cacert.pem" 2>/dev/null || echo_msg "Warning: Failed to set ssl_verify."
# Attempt to get certifi path, suppress error if python or certifi not found yet in base
CERTIFI_PATH=$(python -c "import certifi; print(certifi.where())" 2>/dev/null)
if [ -n "$CERTIFI_PATH" ]; then
    conda env config vars set SSL_CERT_FILE="$CERTIFI_PATH" 2>/dev/null || echo_msg "Warning: Failed to set SSL_CERT_FILE."
else
    echo_msg "Warning: certifi not found or Python error, skipping SSL_CERT_FILE setting for now."
fi


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
#conda run -p "$CONDA_BASE/envs/ViOTUcluster/envs/DRAM" bash -c "echo 'DRAM environment activated'" || { echo_msg "Warning: DRAM conda activate failed, but proceeding."; }
echo_msg "Downloading and replacing DRAM-setup.py script..."
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/scripts/DRAM-setup.py -O "$CONDA_BASE/envs/ViOTUcluster/envs/DRAM/DRAM-setup.py" || { echo_msg "Error: Failed to download DRAM-setup.py"; exit 1; }
echo_msg "DRAM environment setup completed."


# Create iPhop environment
echo_msg "Creating iPhop environment..."

# Try to activate ViOTUcluster environment first
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || {
    echo_msg "Warning: Failed to activate ViOTUcluster before iPhop setup, proceeding."
}
# Determine whether to use mamba or conda
if command -v mamba >/dev/null 2>&1; then
    mamba_cmd_iphop="mamba"
else
    echo_msg "Warning: Mamba not found, falling back to conda for iPhop environment creation."
    mamba_cmd_iphop="conda"
fi
# Create iPhop environment
echo_msg "Creating iPhop environment using $mamba_cmd_iphop..."
"$mamba_cmd_iphop" create -y -c conda-forge -p "$CONDA_BASE/envs/ViOTUcluster/envs/iPhop" python=3.8 mamba || { # Added mamba here as it's often useful in sub-envs too
    echo_msg "Error: Failed to create iPhop environment"
    exit 1
}
# Activate iPhop environment
conda activate "$CONDA_BASE/envs/ViOTUcluster/envs/iPhop" 2>/dev/null || {
    echo_msg "Warning: iPhop conda activate failed, proceeding."
}
# Install iPhop
echo_msg "Installing iPhop using $mamba_cmd_iphop..."
# Ensure the iPhop env is active for the install command or use -p
conda run -p "$CONDA_BASE/envs/ViOTUcluster/envs/iPhop" "$mamba_cmd_iphop" install -y -c conda-forge -c bioconda iphop || {
    echo_msg "Error: Failed to install iphop"
    exit 1
}

# Final activation of ViOTUcluster
echo_msg "Finalizing setup..."
conda activate "$CONDA_BASE/envs/ViOTUcluster" 2>/dev/null || { echo_msg "Warning: Final conda activate failed, environment should still be usable."; }

echo_msg "[✅] ViOTUcluster Setup complete."
echo_msg "Current version: 0.5.2.1 (Script version, not necessarily package version)"

# Clean up downloaded files
echo_msg "Cleaning up downloaded archive files..."
rm -f "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz" 2>/dev/null
rm -f "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz" 2>/dev/null
# rm -f "$CONDA_BASE/envs/ViOTUcluster/DRAM.tar.gz" 2>/dev/null # This was commented out in original

echo_msg "Setup finished."
echo_msg "To activate the main environment, run: conda activate $CONDA_BASE/envs/ViOTUcluster"
echo_msg "Sub-environments like vRhyme, DRAM, iPhop are nested and typically managed by ViOTUcluster scripts."
