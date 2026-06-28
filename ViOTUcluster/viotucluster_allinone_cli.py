#!/usr/bin/env python3
"""
ViOTUcluster_AllinOne - All-in-One Command Line Interface

Starts from raw FASTQ files, handles preprocessing (fastp + assembly),
then runs the complete ViOTUcluster pipeline.
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


BANNER = """
####################################################################################################
 ██╗   ██╗██╗ ██████╗ ████████╗██╗   ██╗ ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ 
 ██║   ██║██║██╔═══██╗╚══██╔══╝██║   ██║██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗
 ██║   ██║██║██║   ██║   ██║   ██║   ██║██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝
 ╚██╗ ██╔╝██║██║   ██║   ██║   ██║   ██║██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗
  ╚████╔╝ ██║╚██████╔╝   ██║   ╚██████╔╝╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║
   ╚═══╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝  ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝
####################################################################################################
                                      ALL-IN-ONE MODE
"""


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='ViOTUcluster_AllinOne',
        description='ViOTUcluster All-in-One: From raw reads to vOTU analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with MEGAHIT assembler
  ViOTUcluster_AllinOne -r raw_reads/ -o output/ -d database/ -a megahit --con

  # Using metaSPAdes with custom threads
  ViOTUcluster_AllinOne -r raw_reads/ -o output/ -d database/ -a metaspades -n 32 --non-con

For more information, visit: https://github.com/liusihang/ViOTUcluster
        """
    )
    
    # Required arguments
    required = parser.add_argument_group('Required Arguments')
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
        help='Output directory for all results'
    )
    required.add_argument(
        '-d', '--database',
        dest='database',
        required=True,
        help='Path to the database directory'
    )
    required.add_argument(
        '-a', '--assembler',
        choices=['megahit', 'metaspades'],
        required=True,
        help='Assembly software to use'
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
        version=f'ViOTUcluster_AllinOne {VERSION}'
    )
    
    return parser


def validate_args(args: argparse.Namespace) -> bool:
    """
    Perform quick validation of parsed arguments.
    
    Returns:
        True if valid, False otherwise
    """
    errors = []
    
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


def preprocessing_outputs_ready(output_dir: str) -> bool:
    """Return True only when preprocessing produced at least one FASTA contig file."""
    contigs_dir = os.path.join(output_dir, "Contigs")
    if not os.path.isdir(contigs_dir):
        return False

    for entry in os.listdir(contigs_dir):
        if entry.lower().endswith((".fa", ".fasta", ".fna")):
            fasta_path = os.path.join(contigs_dir, entry)
            if os.path.isfile(fasta_path) and os.path.getsize(fasta_path) > 0:
                return True
    return False


def run_preprocessing(args) -> bool:
    """
    Run preprocessing: fastp cleaning and assembly.
    
    Returns:
        True if successful
    """
    import logging
    from ViOTUcluster.ContigsPreprocess import main as preprocess_main
    
    logger = logging.getLogger(__name__)
    logger.info("[🔄] Starting preprocessing (fastp + assembly)...")
    
    # Build arguments for ContigsPreprocess
    preprocess_args = [
        '-i', args.raw_seq_dir,
        '-o', args.output_dir,
        '-c', str(args.threads) if args.threads > 0 else str(multiprocessing.cpu_count()),
        '-a', args.assembler,
        '--asm_concurrency', str(args.assemble_jobs),
    ]
    
    try:
        exit_code = preprocess_main(preprocess_args)
        if exit_code != 0:
            logger.error(f"[❌] Preprocessing failed with code: {exit_code}")
            return False
        logger.info("[✅] Preprocessing completed")
        return True
    except SystemExit as e:
        if e.code == 0:
            return True
        logger.error(f"[❌] Preprocessing failed with code: {e.code}")
        return False
    except Exception as e:
        logger.error(f"[❌] Preprocessing failed: {e}")
        return False


def main(argv=None) -> int:
    """
    Main entry point for ViOTUcluster_AllinOne.
    
    Args:
        argv: Command line arguments (defaults to sys.argv)
        
    Returns:
        Exit code (0 for success)
    """
    import logging
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)
    
    # Parse arguments
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Show banner
    print(BANNER)
    print(f"ViOTUcluster_AllinOne v{VERSION}")
    print()
    
    # Validate arguments
    if not validate_args(args):
        return 1
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: Run preprocessing
    if not run_preprocessing(args):
        return 1
    
    # Determine paths after preprocessing
    contigs_dir = os.path.join(args.output_dir, "Contigs")
    clean_reads_dir = os.path.join(args.output_dir, "Cleanreads")
    
    # Check if preprocessing produced expected outputs
    if not preprocessing_outputs_ready(args.output_dir):
        logger.error(
            f"[❌] Preprocessing finished without usable contigs under: {contigs_dir}"
        )
        return 1
    
    # Determine concentration type
    concentration_type = "concentration" if args.concentration else "non-concentration"
    
    # Step 2: Run main pipeline
    from ViOTUcluster.pipeline import run_pipeline
    
    return run_pipeline(
        input_dir=contigs_dir,
        raw_seq_dir=clean_reads_dir,
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
