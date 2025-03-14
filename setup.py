#!/usr/bin/env python3
import sys
import subprocess
import glob
from setuptools import setup, find_packages
from packaging.version import Version

# Check Python version, must be 3.8.x
if sys.version_info[:2] != (3, 8):
    sys.exit("This package requires Python 3.8.")

try:
    output = subprocess.check_output(["mamba", "--version"], universal_newlines=True)
    version_str = output.strip().split()[-1]
    required_mamba_version = Version("1.5.1")
    installed_mamba_version = Version(version_str)
    if installed_mamba_version < required_mamba_version:
        sys.exit("mamba version must be = 1.5.1")
except Exception as e:
    sys.exit("mamba not detected or version check failed: " + str(e))

setup(
    name="ViOTUcluster",
    version="0.4.4",
    packages=find_packages(),
    include_package_data=True,
    scripts=glob.glob("Modules/*"),
    license="GPL-2.0",
    license_files=["LICENSE"],
    python_requires=">=3.8, <3.9",
    author="Sihang Liu",
    description="ViOTUcluster: A high-speed, all-in-one solution that streamlines the entire virome analysis workflow",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown", 
    package_data={'': ['Modules/*']},
)
