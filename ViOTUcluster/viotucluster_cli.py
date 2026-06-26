#!/usr/bin/env python3
"""
ViOTUcluster - Command Line Interface

A high-speed, all-in-one solution for viromic analysis from metagenomic data.

This is the main entry point for the ViOTUcluster pipeline.
"""

import argparse
import sys
import os
import multiprocessing

# Add parent directory to path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster.config import (
    VERSION,
    DEFAULT_MIN_LENGTH,
    DEFAULT_MAX_PREDICTION_TASKS,
    DEFAULT_TPM_TASKS,
    DEFAULT_ASSEMBLE_JOBS,
    DEFAULT_MODULE_TIMEOUT_SECONDS,
)


BANNER = r"""
####################################################################################################
██╗   ██╗██╗ ██████╗ ████████╗██╗   ██╗ ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ 
██║   ██║██║██╔═══██╗╚══██╔══╝██║   ██║██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗
██║   ██║██║██║   ██║   ██║   ██║   ██║██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝
╚██╗ ██╔╝██║██║   ██║   ██║   ██║   ██║██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗
 ╚████╔╝ ██║╚██████╔╝   ██║   ╚██████╔╝╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║
  ╚═══╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝  ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝
####################################################################################################
"""


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='ViOTUcluster',
        description='ViOTUcluster: A high-speed, all-in-one solution for viromic analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with assembled contigs
  ViOTUcluster -i contigs/ -r raw_reads/ -o output/ -d database/ --con

  # With custom thread count and minimum length
  ViOTUcluster -i contigs/ -r raw_reads/ -o output/ -d database/ -n 32 -m 3000 --non-con

  # Disable binning for problematic samples
  ViOTUcluster -i contigs/ -r raw_reads/ -o output/ -d database/ --con --disable-binning

For more information, visit: https://github.com/liusihang/ViOTUcluster
        """
    )
    
    # Required arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument(
        '-i', '--input',
        dest='input_dir',
        required=True,
        help='Input directory containing assembled contig files (FASTA format)'
    )
    required.add_argument(
        '-r', '--raw-seqs',
        dest='raw_seq_dir',
        required=True,
        help='Directory containing raw paired-end FASTQ files'
    )
    required.add_argument(
        '-o', '--output',
        dest='output_dir',
        required=True,
        help='Output directory for results'
    )
    required.add_argument(
        '-d', '--database',
        dest='database',
        required=True,
        help='Path to the database directory'
    )
    
    # Concentration type (mutually exclusive)
    conc_group = parser.add_mutually_exclusive_group(required=True)
    conc_group.add_argument(
        '--con',
        dest='concentration',
        action='store_true',
        help='Use for samples enriched via viral-particle concentration methods'
    )
    conc_group.add_argument(
        '--non-con',
        dest='non_concentration',
        action='store_true',
        help='Use for samples NOT enriched (low viral proportion)'
    )
    
    # Optional arguments
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument(
        '-n', '--threads',
        type=int,
        default=0,
        help=f'Number of threads to use (default: all available = {multiprocessing.cpu_count()})'
    )
    optional.add_argument(
        '-m', '--min-length',
        type=int,
        default=DEFAULT_MIN_LENGTH,
        help=f'Minimum sequence length in bp (default: {DEFAULT_MIN_LENGTH})'
    )
    optional.add_argument(
        '--reassemble',
        action='store_true',
        help='Enable reassembly of bins (beta feature, increases runtime)'
    )
    optional.add_argument(
        '--disable-binning',
        action='store_true',
        help='Skip the vRhyme binning stage'
    )
    optional.add_argument(
        '--save-sambamba-intermediate',
        action='store_true',
        help='Keep Sambamba intermediate BAM files (for debugging)'
    )
    
    # Concurrency controls
    perf = parser.add_argument_group('Performance Tuning')
    perf.add_argument(
        '-P', '--max-prediction-tasks',
        type=int,
        default=DEFAULT_MAX_PREDICTION_TASKS,
        help=f'Max concurrent prediction jobs (default: {DEFAULT_MAX_PREDICTION_TASKS})'
    )
    perf.add_argument(
        '-T', '--tpm-tasks',
        type=int,
        default=DEFAULT_TPM_TASKS,
        help=f'Max concurrent BAM/TPM processing samples (default: {DEFAULT_TPM_TASKS})'
    )
    perf.add_argument(
        '-A', '--assemble-jobs',
        type=int,
        default=DEFAULT_ASSEMBLE_JOBS,
        help=f'Max concurrent assembly samples (default: {DEFAULT_ASSEMBLE_JOBS})'
    )
    perf.add_argument(
        '--module-timeout-hours',
        type=int,
        default=DEFAULT_MODULE_TIMEOUT_SECONDS // 3600,
        help='Abort a pipeline stage if it runs longer than this many hours (0 disables timeout)'
    )
    
    # Version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'ViOTUcluster {VERSION}'
    )
    
    return parser


def validate_args(args: argparse.Namespace) -> bool:
    """
    Perform quick validation of parsed arguments.
    
    Returns:
        True if valid, False otherwise
    """
    errors = []
    
    # Check input directory
    if not os.path.isdir(args.input_dir):
        errors.append(f"Input directory not found: {args.input_dir}")
    
    # Check raw sequences directory
    if not os.path.isdir(args.raw_seq_dir):
        errors.append(f"Raw sequences directory not found: {args.raw_seq_dir}")
    
    # Check database directory
    if not os.path.isdir(args.database):
        errors.append(f"Database directory not found: {args.database}")
    
    # Validate numeric parameters
    if args.min_length < 0:
        errors.append(f"min-length must be non-negative, got: {args.min_length}")
    
    if args.threads < 0:
        errors.append(f"threads must be non-negative, got: {args.threads}")
    
    if args.max_prediction_tasks < 1:
        errors.append(f"max-prediction-tasks must be >= 1, got: {args.max_prediction_tasks}")
    
    if args.tpm_tasks < 1:
        errors.append(f"tpm-tasks must be >= 1, got: {args.tpm_tasks}")
    
    if args.assemble_jobs < 1:
        errors.append(f"assemble-jobs must be >= 1, got: {args.assemble_jobs}")

    if args.module_timeout_hours < 0:
        errors.append(
            f"module-timeout-hours must be non-negative, got: {args.module_timeout_hours}"
        )
    
    if errors:
        print("Error: Invalid arguments", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return False
    
    return True


def main(argv=None) -> int:
    """
    Main entry point for ViOTUcluster.
    
    Args:
        argv: Command line arguments (defaults to sys.argv)
        
    Returns:
        Exit code (0 for success)
    """
    # Parse arguments
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Show banner
    print(BANNER)
    print(f"ViOTUcluster v{VERSION}")
    print()
    
    # Validate arguments
    if not validate_args(args):
        return 1
    
    # Determine concentration type
    concentration_type = "concentration" if args.concentration else "non-concentration"
    
    # Import pipeline (delayed to speed up --help)
    from ViOTUcluster.pipeline import run_pipeline
    
    # Run the pipeline
    return run_pipeline(
        input_dir=args.input_dir,
        raw_seq_dir=args.raw_seq_dir,
        output_dir=args.output_dir,
        database=args.database,
        threads=args.threads,
        min_length=args.min_length,
        concentration_type=concentration_type,
        reassemble=args.reassemble,
        disable_binning=args.disable_binning,
        save_sambamba_intermediate=args.save_sambamba_intermediate,
        max_prediction_tasks=args.max_prediction_tasks,
        tpm_tasks=args.tpm_tasks,
        assemble_jobs=args.assemble_jobs,
        module_timeout_seconds=(
            None if args.module_timeout_hours == 0
            else args.module_timeout_hours * 3600
        ),
    )


if __name__ == '__main__':
    sys.exit(main())
