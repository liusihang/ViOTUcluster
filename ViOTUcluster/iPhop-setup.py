#!/usr/bin/env python
import os
import subprocess
import sys

def download_iphop_db(db_dir):
    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if not conda_prefix:
            raise EnvironmentError("Conda environment is not activated.")
        
        env_path = os.path.join(conda_prefix, "envs", "iPhop")
        
        # Define the command to download the iPhop database
        iphops_cmd = [
            'conda', 'run',
            '-p', env_path,
            'iphop', 'download',
            '--no_prompt',
            '--db_dir', db_dir
        ]
        # Run the command to download the database
        subprocess.run(iphops_cmd, check=True)
        print(f"[✅] iPhop database download completed. Output directory: {db_dir}")
    
    except subprocess.CalledProcessError as e:
        print(f"[❌] iPhop database download failed: {e}")
    except Exception as e:
        print(f"[❌] Error: {e}")

# Example usage
if len(sys.argv) != 2:
    print("Usage: iPhop_setup.py <db_directory>")
    sys.exit(1)

db_directory = sys.argv[1]  # The directory where the database will be downloaded

download_iphop_db(db_directory)
