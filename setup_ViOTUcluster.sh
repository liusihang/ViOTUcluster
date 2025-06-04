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

# --- Script Variables ---
INSTALL_PREFIX="" # Custom installation prefix
USE_CHINA_MIRROR=false
VIOTUCLUSTER_ENV_NAME="ViOTUcluster" # Name of the main conda environment to be created

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

# --- Option Parsing ---
show_help() {
    echo_msg "Usage: $0 [options]"
    echo_msg "Options:"
    echo_msg "  --china                        Use download mirrors in China (china.scidb.cn)."
    echo_msg "  -p, --prefix, --install-path <path>  Specify a custom base directory for installation."
    echo_msg "                                 ViOTUcluster will be installed into <path>/$VIOTUCLUSTER_ENV_NAME."
    echo_msg "  -h, --help                     Show this help message."
}

ORIGINAL_ARGS=("$@")
UNRECOGNIZED_ARGS=()

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --china)
        USE_CHINA_MIRROR=true
        shift
        ;;
        -p|--prefix|--install-path)
        if [[ -z "$2" || "$2" == -* ]]; then
            echo_msg "Error: Argument for $1 is missing or looks like another option."
            show_help
            exit 1
        fi
        INSTALL_PREFIX="$2"
        shift
        shift
        ;;
        -h|--help)
        show_help
        exit 0
        ;;
        *)
        UNRECOGNIZED_ARGS+=("$1")
        shift
        ;;
    esac
done

if [ ${#UNRECOGNIZED_ARGS[@]} -gt 0 ]; then
    echo_msg "Ignoring unrecognized arguments: ${UNRECOGNIZED_ARGS[*]}"
fi

if [ "$USE_CHINA_MIRROR" = true ]; then
    echo_msg "Using China download mirror (china.scidb.cn)."
    VIOTUCLUSTER_URL="$CHINA_VIOTUCLUSTER_URL"
    VRHYME_URL="$CHINA_VRHYME_URL"
else
    echo_msg "Proceeding with default download URLs from Zenodo."
    echo_msg "Tip: You can use '$0 --china' to use download mirrors in China."
fi
# --- End Option Parsing ---


CONDA_BASE_DETECTED=$(conda info --base 2>/dev/null) || { echo_msg "Error: Conda not found. Please install Conda first."; exit 1; }

if [ -z "$CONDA_BASE_DETECTED" ]; then
    echo_msg "Error: Could not determine Conda executable base path."
    exit 1
fi

if [ -f "$CONDA_BASE_DETECTED/etc/profile.d/conda.sh" ]; then
    . "$CONDA_BASE_DETECTED/etc/profile.d/conda.sh"
else
    echo_msg "Error: Conda initialization script not found at $CONDA_BASE_DETECTED/etc/profile.d/conda.sh"
    exit 1
fi

if [ -n "$INSTALL_PREFIX" ]; then
    mkdir -p "$INSTALL_PREFIX" || { echo_msg "Error: Failed to create custom install prefix directory '$INSTALL_PREFIX'."; exit 1; }
    INSTALL_PREFIX_ABS=$(cd "$INSTALL_PREFIX"; pwd)
    VIOTUCLUSTER_ENV_ROOT="$INSTALL_PREFIX_ABS/$VIOTUCLUSTER_ENV_NAME"
    echo_msg "Custom installation target: $VIOTUCLUSTER_ENV_ROOT"
else
    VIOTUCLUSTER_ENV_ROOT="$CONDA_BASE_DETECTED/envs/$VIOTUCLUSTER_ENV_NAME"
    echo_msg "Default Conda installation target: $VIOTUCLUSTER_ENV_ROOT"
fi

# Marker files for idempotency
# Using a dot-directory inside the env root for markers to keep it tidy
MARKER_DIR="$VIOTUCLUSTER_ENV_ROOT/.install_markers"
mkdir -p "$MARKER_DIR"

VIOTUCLUSTER_DOWNLOADED_MARKER="$MARKER_DIR/viotucluster_downloaded"
VIOTUCLUSTER_UNPACKED_MARKER="$MARKER_DIR/viotucluster_unpacked"
SSL_CONFIGURED_MARKER="$MARKER_DIR/ssl_configured"
VRHYME_DOWNLOADED_MARKER="$MARKER_DIR/vrhyme_downloaded"
VRHYME_UNPACKED_MARKER="$MARKER_DIR/vrhyme_unpacked" # Path specific to vRhyme env for its marker
DRAM_ENV_YAML_DOWNLOADED_MARKER="$MARKER_DIR/dram_env_yaml_downloaded"
DRAM_ENV_CREATED_MARKER="$MARKER_DIR/dram_env_created"
DRAM_SETUP_PY_DOWNLOADED_MARKER="$MARKER_DIR/dram_setup_py_downloaded"
IPHOP_ENV_SHELL_CREATED_MARKER="$MARKER_DIR/iphop_env_shell_created"
IPHOP_PACKAGE_INSTALLED_MARKER="$MARKER_DIR/iphop_package_installed"


echo_msg "Creating ViOTUcluster environment directory at $VIOTUCLUSTER_ENV_ROOT..."
mkdir -p "$VIOTUCLUSTER_ENV_ROOT" || { echo_msg "Error: Failed to create directory '$VIOTUCLUSTER_ENV_ROOT'."; exit 1; }

echo_msg "Downloading ViOTUcluster, vRhyme, and DRAM packages..." # Original message placement

# --- ViOTUcluster Download ---
VIOTUCLUSTER_ARCHIVE_PATH="$VIOTUCLUSTER_ENV_ROOT/ViOTUcluster.tar.gz"
if [ -f "$VIOTUCLUSTER_DOWNLOADED_MARKER" ] && [ -f "$VIOTUCLUSTER_ARCHIVE_PATH" ]; then
    echo_msg "Skipping ViOTUcluster download: $VIOTUCLUSTER_ARCHIVE_PATH already exists and marked as downloaded."
else
    echo_msg "Downloading ViOTUcluster from: $VIOTUCLUSTER_URL"
    wget -q "$VIOTUCLUSTER_URL" -O "$VIOTUCLUSTER_ARCHIVE_PATH" || { echo_msg "Error: Failed to download ViOTUcluster.tar.gz"; exit 1; }
    touch "$VIOTUCLUSTER_DOWNLOADED_MARKER"
fi

# --- vRhyme Download ---
VRHYME_ARCHIVE_PATH="$VIOTUCLUSTER_ENV_ROOT/vRhyme.tar.gz"
if [ -f "$VRHYME_DOWNLOADED_MARKER" ] && [ -f "$VRHYME_ARCHIVE_PATH" ]; then
    echo_msg "Skipping vRhyme download: $VRHYME_ARCHIVE_PATH already exists and marked as downloaded."
else
    echo_msg "Downloading vRhyme from: $VRHYME_URL"
    wget -q "$VRHYME_URL" -O "$VRHYME_ARCHIVE_PATH" || { echo_msg "Error: Failed to download vRhyme.tar.gz"; exit 1; }
    touch "$VRHYME_DOWNLOADED_MARKER"
fi

# --- ViOTUcluster Extract and Unpack ---
# Check for a key file/dir that indicates successful unpack, like conda-meta/history or bin directory
VIOTUCLUSTER_CONDA_META_HISTORY="$VIOTUCLUSTER_ENV_ROOT/conda-meta/history"
if [ -f "$VIOTUCLUSTER_UNPACKED_MARKER" ] && [ -f "$VIOTUCLUSTER_CONDA_META_HISTORY" ]; then
    echo_msg "Skipping ViOTUcluster extract and unpack: Already marked as unpacked and history file exists."
else
    if [ ! -d "$VIOTUCLUSTER_ENV_ROOT/bin" ] || [ ! -f "$VIOTUCLUSTER_CONDA_META_HISTORY" ]; then # If not already extracted (basic check)
        echo_msg "Extracting ViOTUcluster environment..."
        tar -xzf "$VIOTUCLUSTER_ARCHIVE_PATH" -C "$VIOTUCLUSTER_ENV_ROOT" || { echo_msg "Error: Failed to extract ViOTUcluster.tar.gz"; exit 1; }
    fi
    # Activate the ViOTUcluster environment BEFORE SSL configuration and other operations that depend on it being active
    conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || { echo_msg "Warning: conda activate failed for ViOTUcluster, but proceeding with unpack."; }
    conda-unpack 2>/dev/null || { echo_msg "Error: Failed to unpack ViOTUcluster environment"; exit 1; }
    touch "$VIOTUCLUSTER_UNPACKED_MARKER"
fi


# --- SSL Configuration ---
# Check if SSL was already configured (e.g., based on a successful run of these commands)
# A simple marker is sufficient here as the commands themselves are idempotent if certs are already set.
if [ -f "$SSL_CONFIGURED_MARKER" ]; then
    echo_msg "Skipping SSL certificate verification configuration: Already marked as configured."
else
    echo_msg "Configuring SSL certificate verification..."
    conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || { echo_msg "Warning: conda activate failed for ViOTUcluster, but proceeding with SSL config."; } # Original warning text
    conda config --env --set ssl_verify "$VIOTUCLUSTER_ENV_ROOT/ssl/cacert.pem" 2>/dev/null || echo_msg "Warning: Failed to set ssl_verify."
    conda env config vars set SSL_CERT_FILE=$(python -c "import certifi; print(certifi.where())" 2>/dev/null) 2>/dev/null || echo_msg "Warning: Failed to set SSL_CERT_FILE."
    touch "$SSL_CONFIGURED_MARKER"
fi
# --- End SSL Configuration ---


# --- vRhyme Environment ---
VRHYME_ENV_PATH="$VIOTUCLUSTER_ENV_ROOT/envs/vRhyme"
VRHYME_CONDA_META_HISTORY="$VRHYME_ENV_PATH/conda-meta/history"
if [ -f "$VRHYME_UNPACKED_MARKER" ] && [ -f "$VRHYME_CONDA_META_HISTORY" ]; then
    echo_msg "Skipping vRhyme environment creation and unpack: Already marked and history file exists."
else
    echo_msg "Creating and unpacking vRhyme environment at $VRHYME_ENV_PATH..."
    mkdir -p "$VRHYME_ENV_PATH" || { echo_msg "Error: Failed to create vRhyme directory."; exit 1; }
    if [ ! -d "$VRHYME_ENV_PATH/bin" ] || [ ! -f "$VRHYME_CONDA_META_HISTORY" ]; then # If not already extracted
        tar -xzf "$VRHYME_ARCHIVE_PATH" -C "$VRHYME_ENV_PATH" || { echo_msg "Error: Failed to extract vRhyme.tar.gz"; exit 1; }
    fi
    # Temporarily activate vRhyme for unpack
    (
      conda activate "$VRHYME_ENV_PATH" 2>/dev/null || echo_msg "Warning: vRhyme conda activate failed for unpack, but proceeding."
      conda-unpack 2>/dev/null || { echo_msg "Error: Failed to unpack vRhyme environment"; exit 1; } # exit 1 will exit subshell only
    ) || { echo_msg "Error during vRhyme unpack subshell."; exit 1; } # Propagate error
    touch "$VRHYME_UNPACKED_MARKER"
fi
# Ensure main env is active again
conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || echo_msg "Warning: Failed to reactivate $VIOTUCLUSTER_ENV_NAME after vRhyme."


# --- DRAM Environment ---
DRAM_ENV_PATH="$VIOTUCLUSTER_ENV_ROOT/envs/DRAM"
DRAM_ENV_YAML_FILENAME="environment.yaml" # Original script downloads to current dir
DRAM_SETUP_PY_PATH="$DRAM_ENV_PATH/DRAM-setup.py"
DRAM_CONDA_META_DIR="$DRAM_ENV_PATH/conda-meta"

echo_msg "Creating DRAM environment at $DRAM_ENV_PATH..." # Original message placement
mkdir -p "$DRAM_ENV_PATH" || { echo_msg "Error: Failed to create DRAM directory."; exit 1; }

# Ensure ViOTUcluster is active (as in original script via deactivate then activate)
conda deactivate 2>/dev/null # Suppress errors if not active
conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || { echo_msg "Warning: conda activate failed for ViOTUcluster before DRAM setup."; }


if [ -f "$DRAM_ENV_CREATED_MARKER" ] && [ -d "$DRAM_CONDA_META_DIR" ]; then
    echo_msg "Skipping DRAM environment creation: Already marked as created and conda-meta exists."
else
    if [ -f "$DRAM_ENV_YAML_DOWNLOADED_MARKER" ] && [ -f "$DRAM_ENV_YAML_FILENAME" ]; then
        echo_msg "Skipping DRAM environment.yaml download: File exists and marked."
    else
        wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml -O "$DRAM_ENV_YAML_FILENAME" || { echo_msg "Error: Failed to download environment.yaml"; exit 1; }
        touch "$DRAM_ENV_YAML_DOWNLOADED_MARKER"
    fi

    # Determine mamba or conda for DRAM creation
    DRAM_CREATE_CMD="conda"
    if command -v mamba >/dev/null 2>&1; then
        DRAM_CREATE_CMD="mamba"
    fi
    echo_msg "Using $DRAM_CREATE_CMD to create DRAM environment..." # Added for clarity
    "$DRAM_CREATE_CMD" env create -f "$DRAM_ENV_YAML_FILENAME" -p "$DRAM_ENV_PATH" || { echo_msg "Error: Failed to create DRAM environment"; exit 1; }
    touch "$DRAM_ENV_CREATED_MARKER"
fi

if [ -f "$DRAM_SETUP_PY_DOWNLOADED_MARKER" ] && [ -f "$DRAM_SETUP_PY_PATH" ]; then
    echo_msg "Skipping DRAM-setup.py download: File exists and marked."
else
    echo_msg "Downloading and replacing DRAM-setup.py script..."
    wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/scripts/DRAM-setup.py -O "$DRAM_SETUP_PY_PATH" || { echo_msg "Error: Failed to download DRAM-setup.py"; exit 1; }
    touch "$DRAM_SETUP_PY_DOWNLOADED_MARKER"
fi
echo_msg "DRAM environment setup completed."


# --- iPhop Environment ---
IPHOP_ENV_PATH="$VIOTUCLUSTER_ENV_ROOT/envs/iPhop"
IPHOP_CONDA_META_DIR="$IPHOP_ENV_PATH/conda-meta"
# A way to check if iphop package is installed: check for a key binary or use conda list
IPHOP_INSTALLED_CHECK_FILE="$IPHOP_ENV_PATH/bin/iphop" # Example check

echo_msg "Creating iPhop environment at $IPHOP_ENV_PATH..." # Original message placement

# Ensure main ViOTUcluster env is active
conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || { echo_msg "Warning: conda activate failed for ViOTUcluster before iPhop setup."; }

# Determine mamba or conda for iPhop operations (as per original script)
mamba_cmd_iphop="conda"
if command -v mamba >/dev/null 2>&1; then
    mamba_cmd_iphop="mamba"
else
    echo_msg "Warning: Mamba not found, falling back to conda for iPhop environment creation." # Original message
fi


if [ -f "$IPHOP_ENV_SHELL_CREATED_MARKER" ] && [ -d "$IPHOP_CONDA_META_DIR" ]; then
    echo_msg "Skipping iPhop environment shell creation: Already marked and conda-meta exists."
else
    echo_msg "Creating iPhop environment shell using $mamba_cmd_iphop..."
    "$mamba_cmd_iphop" create -y -c conda-forge -p "$IPHOP_ENV_PATH" python=3.8 mamba || {
        echo_msg "Error: Failed to create iPhop environment shell"
        exit 1
    }
    touch "$IPHOP_ENV_SHELL_CREATED_MARKER"
fi


if [ -f "$IPHOP_PACKAGE_INSTALLED_MARKER" ] && [ -f "$IPHOP_INSTALLED_CHECK_FILE" ]; then
    echo_msg "Skipping iPhop package installation: Already marked and check file exists."
    # For a more robust check, one might run `conda run -p "$IPHOP_ENV_PATH" conda list iphop`
    # but for simplicity and adherence to "don't change other things", a file check is used.
else
    echo_msg "Installing iPhop into $IPHOP_ENV_PATH using $mamba_cmd_iphop..."
    conda activate "$IPHOP_ENV_PATH" 2>/dev/null || { echo_msg "Warning: conda activate failed for iPhop environment, but proceeding with iphop install."; } # Original had different warning here
    "$mamba_cmd_iphop" install -y -c conda-forge -c bioconda iphop || {
        echo_msg "Error: Failed to install iphop package"
        exit 1
    }
    touch "$IPHOP_PACKAGE_INSTALLED_MARKER"
fi


echo_msg "Finalizing setup..."
# Ensure ViOTUcluster is the active environment at the end
conda activate "$VIOTUCLUSTER_ENV_ROOT" 2>/dev/null || { echo_msg "Warning: Final conda activate failed, environment should still be usable."; }

echo_msg "[✅] ViOTUcluster Setup complete."
echo_msg "Current version: 0.5.2.2"

echo_msg "Cleaning up downloaded archive files..."
if [ -f "$VIOTUCLUSTER_UNPACKED_MARKER" ]; then
    rm -f "$VIOTUCLUSTER_ARCHIVE_PATH" 2>/dev/null
fi
if [ -f "$VRHYME_UNPACKED_MARKER" ]; then
    rm -f "$VRHYME_ARCHIVE_PATH" 2>/dev/null
fi
if [ -f "$DRAM_ENV_CREATED_MARKER" ] && [ -f "$DRAM_ENV_YAML_DOWNLOADED_MARKER" ] ; then # If DRAM env was created using the yaml
    rm -f "$DRAM_ENV_YAML_FILENAME" 2>/dev/null
fi


echo_msg "Setup finished."
echo_msg "To activate the main environment, run: conda activate $VIOTUCLUSTER_ENV_ROOT"
