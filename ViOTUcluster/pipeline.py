#!/usr/bin/env python3
"""
ViOTUcluster Pipeline Orchestration Module.

This module handles the execution flow of the ViOTUcluster pipeline,
coordinating between Python processing steps and Shell module calls.
"""

import os
import sys
import subprocess
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any

from .config import (
    DEFAULT_MIN_LENGTH,
    DEFAULT_MAX_PREDICTION_TASKS,
    DEFAULT_TPM_TASKS,
    DEFAULT_ASSEMBLE_JOBS,
    DB_VIRSORTER,
    DB_VIRALVERIFY,
    DB_CHECKV,
    DB_GENOMAD,
    get_threads,
    get_virsorter_groups,
)


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class PipelineError(Exception):
    """Custom exception for pipeline errors."""
    pass


class ViOTUclusterPipeline:
    """
    Main pipeline orchestrator for ViOTUcluster.
    
    Handles argument validation, environment setup, and module execution.
    """
    
    def __init__(
        self,
        input_dir: str,
        raw_seq_dir: str,
        output_dir: str,
        database: str,
        threads: int = 0,
        min_length: int = DEFAULT_MIN_LENGTH,
        concentration_type: str = "non-concentration",
        reassemble: bool = False,
        disable_binning: bool = False,
        save_sambamba_intermediate: bool = False,
        max_prediction_tasks: int = DEFAULT_MAX_PREDICTION_TASKS,
        tpm_tasks: int = DEFAULT_TPM_TASKS,
        assemble_jobs: int = DEFAULT_ASSEMBLE_JOBS,
        sample_type: str = "Mix",
    ):
        """
        Initialize the pipeline with configuration.
        
        Args:
            input_dir: Directory containing assembled contig files
            raw_seq_dir: Directory with raw sequencing FASTQ files
            output_dir: Output directory for results
            database: Path to database directory
            threads: Number of threads (0 = auto-detect)
            min_length: Minimum sequence length for filtering
            concentration_type: 'concentration' or 'non-concentration'
            reassemble: Enable reassembly of bins
            disable_binning: Skip vRhyme binning stage
            save_sambamba_intermediate: Keep intermediate BAM files
            max_prediction_tasks: Max concurrent prediction jobs
            tpm_tasks: Max concurrent BAM/TPM processing samples
            assemble_jobs: Max concurrent assembly samples
            sample_type: 'DNA', 'RNA', or 'Mix'
        """
        self.input_dir = os.path.abspath(input_dir)
        self.raw_seq_dir = os.path.abspath(raw_seq_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.database = os.path.abspath(database)
        self.threads = get_threads(threads)
        self.min_length = min_length
        self.concentration_type = concentration_type
        self.reassemble = reassemble
        self.disable_binning = disable_binning
        self.save_sambamba_intermediate = save_sambamba_intermediate
        self.max_prediction_tasks = max_prediction_tasks
        self.tpm_tasks = tpm_tasks
        self.assemble_jobs = assemble_jobs
        self.sample_type = sample_type
        self.group = get_virsorter_groups(sample_type)
        
        # Derived paths
        self.log_dir = os.path.join(self.output_dir, "Log")
        self.log_file = os.path.join(self.output_dir, "pipeline.log")
        self.filtered_seqs_dir = os.path.join(self.output_dir, "FilteredSeqs")
        
        # Get script directory (where Shell modules are located)
        self.script_dir = self._find_script_dir()
        
        # Track start time
        self.start_time = None
        
    def _find_script_dir(self) -> str:
        """Find the directory containing Shell modules."""
        # Try relative to this file
        this_dir = Path(__file__).parent
        
        # Check if Modules directory exists at expected locations
        possible_paths = [
            this_dir.parent / "Modules",  # Development layout
            this_dir / "Modules",
            Path(sys.prefix) / "Modules",  # Installed layout
        ]
        
        for path in possible_paths:
            if path.exists() and (path / "viral_prediction_module.sh").exists():
                return str(path)
        
        # Fallback: use environment variable if set
        if "ScriptDir" in os.environ:
            return os.environ["ScriptDir"]
        
        raise PipelineError(
            "Cannot find Shell modules directory. "
            "Ensure ViOTUcluster is properly installed."
        )
    
    def _setup_environment(self) -> Dict[str, str]:
        """Setup environment variables for Shell modules."""
        env = os.environ.copy()
        
        env.update({
            "INPUT_DIR": self.input_dir,
            "RAW_SEQ_DIR": self.raw_seq_dir,
            "OUTPUT_DIR": self.output_dir,
            "DATABASE": self.database,
            "THREADS": str(self.threads),
            "MIN_LENGTH": str(self.min_length),
            "CONCENTRATION_TYPE": self.concentration_type,
            "Group": self.group,
            "ScriptDir": self.script_dir,
            "REASSEMBLE": str(self.reassemble).lower(),
            "DISABLE_BINNING": str(self.disable_binning).lower(),
            "SAMBAMBA_SAVE_INTERMEDIATE": str(self.save_sambamba_intermediate).lower(),
            "MAX_PredictionTASKS": str(self.max_prediction_tasks),
            "MAX_TASKS": str(self.max_prediction_tasks),
            "TPM_tasks": str(self.tpm_tasks),
            "Assemble_jobs": str(self.assemble_jobs),
        })
        
        # Find and set FILES variable (list of FASTA files)
        files = self._find_fasta_files(self.filtered_seqs_dir)
        env["FILES"] = " ".join(files)
        
        return env
    
    def _find_fasta_files(self, directory: str) -> list:
        """Find all FASTA files in a directory."""
        if not os.path.exists(directory):
            return []
        
        files = []
        for f in os.listdir(directory):
            if f.lower().endswith(('.fa', '.fasta', '.fna')):
                files.append(os.path.join(directory, f))
        return sorted(files)
    
    def _run_module(
        self, 
        name: str, 
        script_path: str, 
        env: Dict[str, str],
        extra_args: str = ""
    ) -> bool:
        """
        Run a Shell module with logging.
        
        Args:
            name: Human-readable module name
            script_path: Path to the Shell script
            env: Environment variables
            extra_args: Additional command line arguments
            
        Returns:
            True if successful, False otherwise
        """
        logger.info(f"[🔄] Starting module: {name}")
        
        module_log = os.path.join(self.log_dir, f"{name.replace(' ', '_')}.log")
        
        cmd = f"bash {script_path} {extra_args}"
        
        try:
            with open(module_log, 'w') as log_file:
                result = subprocess.run(
                    cmd,
                    shell=True,
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=self.output_dir,
                )
            
            if result.returncode != 0:
                logger.error(f"[❌] Module {name} failed. Check log: {module_log}")
                return False
            
            logger.info(f"[✅] Module {name} completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"[❌] Error running module {name}: {e}")
            return False
    
    def _run_python_step(
        self,
        name: str,
        func: callable,
        *args,
        **kwargs
    ) -> bool:
        """
        Run a Python processing step with error handling.
        
        Args:
            name: Step name for logging
            func: Function to call
            *args, **kwargs: Arguments for the function
            
        Returns:
            True if successful
        """
        logger.info(f"[🔄] Running: {name}")
        
        try:
            func(*args, **kwargs)
            logger.info(f"[✅] {name} completed")
            return True
        except Exception as e:
            logger.error(f"[❌] {name} failed: {e}")
            return False
    
    def validate_inputs(self) -> bool:
        """Validate all input paths and parameters."""
        from .validation import (
            validate_directory,
            validate_database_structure,
            validate_positive_integer,
        )
        
        errors = []
        
        # Validate directories
        try:
            validate_directory(self.input_dir, description="Input directory")
        except Exception as e:
            errors.append(str(e))
        
        try:
            validate_directory(self.raw_seq_dir, description="Raw sequences directory")
        except Exception as e:
            errors.append(str(e))
        
        try:
            validate_directory(self.database, description="Database directory")
        except Exception as e:
            errors.append(str(e))
        
        # Validate min_length
        if self.min_length < 0:
            errors.append(f"min_length must be non-negative, got: {self.min_length}")
        
        # Check concentration type
        if self.concentration_type not in ("concentration", "non-concentration"):
            errors.append(
                f"concentration_type must be 'concentration' or 'non-concentration', "
                f"got: {self.concentration_type}"
            )
        
        if errors:
            for err in errors:
                logger.error(f"[❌] {err}")
            return False
        
        logger.info("[✅] Input validation passed")
        return True
    
    def setup_directories(self) -> bool:
        """Create required output directories."""
        dirs_to_create = [
            self.output_dir,
            self.log_dir,
            self.filtered_seqs_dir,
            os.path.join(self.output_dir, "SeprateFile"),
            os.path.join(self.output_dir, "Summary"),
            os.path.join(self.output_dir, "Summary", "SeperateRes"),
            os.path.join(self.output_dir, "Summary", "SeperateRes", "bins"),
            os.path.join(self.output_dir, "Summary", "SeperateRes", "unbined"),
        ]
        
        try:
            for d in dirs_to_create:
                os.makedirs(d, exist_ok=True)
            logger.info("[✅] Output directories created")
            return True
        except Exception as e:
            logger.error(f"[❌] Failed to create directories: {e}")
            return False
    
    def run_filter_contigs(self) -> bool:
        """Run contig filtering step."""
        from .filter_contigs import filter_sequences_flexible
        
        logger.info(f"[🔄] Filtering sequences (min length: {self.min_length}bp)")
        
        try:
            filter_sequences_flexible(
                self.min_length,
                self.input_dir,
                self.filtered_seqs_dir
            )
            logger.info("[✅] Sequence filtering completed")
            return True
        except Exception as e:
            logger.error(f"[❌] Sequence filtering failed: {e}")
            return False
    
    def run(self) -> int:
        """
        Execute the complete pipeline.
        
        Returns:
            Exit code (0 for success, non-zero for failure)
        """
        self.start_time = time.time()
        
        logger.info("=" * 60)
        logger.info("ViOTUcluster Pipeline Starting")
        logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 60)
        
        # Step 1: Validate inputs
        if not self.validate_inputs():
            return 1
        
        # Step 2: Setup directories
        if not self.setup_directories():
            return 1
        
        # Step 3: Filter contigs
        if not self.run_filter_contigs():
            return 1
        
        # Setup environment for Shell modules
        env = self._setup_environment()
        
        # Step 4: Viral Prediction
        if not self._run_module(
            "Viral Prediction",
            os.path.join(self.script_dir, "viral_prediction_module.sh"),
            env
        ):
            return 1
        
        # Step 5: Cross Validation
        if not self._run_module(
            "Cross Validation",
            os.path.join(self.script_dir, "cross_validation_module.sh"),
            env,
            f"--{self.concentration_type}"
        ):
            return 1
        
        # Step 6: Binning (or skip)
        if self.disable_binning:
            logger.info("[ℹ️] Binning disabled, preparing unbinned contigs")
            # Copy filtered contigs to unbined directory
            self._prepare_unbinned_from_contigs()
        else:
            if not self._run_module(
                "Binning and Merge",
                os.path.join(self.script_dir, "binning_merge_module.sh"),
                env
            ):
                return 1
        
        # Step 7: dRep
        if not self._run_module(
            "dRep",
            os.path.join(self.script_dir, "drep_module.sh"),
            env
        ):
            return 1
        
        # Step 8: Summary
        if not self._run_module(
            "Summary",
            os.path.join(self.script_dir, "summary_module.sh"),
            env
        ):
            return 1
        
        # Calculate runtime
        total_runtime = int(time.time() - self.start_time)
        
        logger.info("=" * 60)
        logger.info("[✅][✅][✅] All basic analysis completed successfully!")
        logger.info(f"Results available in: {self.output_dir}/Summary")
        logger.info(f"Total runtime: {total_runtime} seconds")
        logger.info("=" * 60)
        
        return 0
    
    def _prepare_unbinned_from_contigs(self) -> bool:
        """Copy filtered contigs as unbinned when binning is disabled."""
        import shutil
        
        unbined_dir = os.path.join(
            self.output_dir, "Summary", "SeperateRes", "unbined"
        )
        os.makedirs(unbined_dir, exist_ok=True)
        
        files = self._find_fasta_files(self.filtered_seqs_dir)
        for f in files:
            basename = os.path.basename(f).replace('.fasta', '').replace('.fa', '')
            dest = os.path.join(unbined_dir, f"{basename}_unbined.fasta")
            shutil.copy(f, dest)
        
        logger.info(f"[✅] Copied {len(files)} files as unbinned contigs")
        return True


def run_pipeline(**kwargs) -> int:
    """
    Convenience function to create and run the pipeline.
    
    Args:
        **kwargs: Arguments passed to ViOTUclusterPipeline
        
    Returns:
        Exit code
    """
    try:
        pipeline = ViOTUclusterPipeline(**kwargs)
        return pipeline.run()
    except Exception as e:
        logger.error(f"[❌] Pipeline failed with error: {e}")
        return 1
