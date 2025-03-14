#!/usr/bin/env bash

# 获取 Conda 基本安装路径
CONDA_BASE=$(conda info --base)

# 创建 ViOTUcluster 环境目录
echo "Creating ViOTUcluster environment directory..."
mkdir -p "$CONDA_BASE/envs/ViOTUcluster"

# 下载必要的包（示例：ViOTUcluster, vRhyme, iphop, DRAM）
# 在这里将下载链接替换为实际 URL
echo "Downloading ViOTUcluster, vRhyme, iphop, DRAM packages..."
#wget -q XXXX -O "$CONDA_BASE/envs/ViOTUcluster/ViOTUcluster.tar.gz"
#wget -q XXXX -O "$CONDA_BASE/envs/ViOTUcluster/vRhyme.tar.gz"

# 解压下载的文件到指定目录
echo "Extracting files..."
tar -xzf "ViOTUcluster.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster"
#tar -xzf "vRhyme.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster"

# 激活 ViOTUcluster 环境
echo "Activating ViOTUcluster environment..."
#source $(conda info --base)/etc/profile.d/conda.sh
source "$CONDA_BASE/envs/ViOTUcluster/bin/activate"

# 解压并准备 Conda 环境
echo "Unpacking Conda environment..."
conda unpack

# 设置 SSL 证书验证路径
echo "Configuring SSL certificate verification..."
conda config --env --set ssl_verify "$CONDA_PREFIX/ssl/cacert.pem"
conda env config vars set SSL_CERT_FILE=$(python -c "import certifi; print(certifi.where())")

# 如果需要重新配置证书，也可以选择下面的命令
# conda env config vars set ssl_verify "$CONDA_PREFIX/ssl/cacert.pem"

# 创建 vRhyme 环境目录
echo "Creating vRhyme environment directory..."
mkdir -p "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme"

# 解压 vRhyme 包到指定目录
echo "Extracting vRhyme package..."
tar -xzf "vRhyme.tar.gz" -C "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme"

# 激活 vRhyme 环境
echo "Activating vRhyme environment..."
source "$CONDA_BASE/envs/ViOTUcluster/envs/vRhyme/bin/activate"

# 解压并准备 vRhyme 环境
echo "Unpacking vRhyme Conda environment..."
conda unpack

# 最后重新激活 ViOTUcluster 环境
source $(conda info --base)/etc/profile.d/conda.sh

echo "[✅] ViOTUcluster Setup complete."
echo "Current version: 0.3.6"