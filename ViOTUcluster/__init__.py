"""
ViOTUcluster - A high-speed, all-in-one solution for viromic analysis.

This package provides tools for processing metagenomic data to identify
and cluster viral operational taxonomic units (vOTUs).
"""

from .config import VERSION

__version__ = VERSION
__author__ = "Sihang Liu"
__email__ = "liusihang@tongji.edu.cn"

__all__ = [
    'config',
    'validation',
    'pipeline',
    'filter_contigs',
    'viotucluster_cli',
]
