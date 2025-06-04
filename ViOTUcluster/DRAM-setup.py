#!/usr/bin/env python3
import os
import subprocess
import sys
import urllib.request
import stat

def run_dram_setup(action, file_dir):
    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if not conda_prefix:
            raise EnvironmentError("Conda environment is not activated.")
        
        env_path = os.path.join(conda_prefix, "envs", "DRAM")
        bin_path = os.path.join(env_path, "bin") 
        target   = os.path.join(bin_path, "DRAM-setup.py")
        url = "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/scripts/DRAM-setup.py"
        with urllib.request.urlopen(url) as resp, open(target, "wb") as f:
            f.write(resp.read())

        current = os.stat(target).st_mode
        os.chmod(target, current | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        # Define the command based on the action type
        if action == "download":
            # Command for downloading the database
            dram_cmd = [
                'conda', 'run',
                '-p', env_path,
                'DRAM-setup.py', 'prepare_databases',
                '--output_dir', file_dir
            ]
        elif action == "import_config":
            # Command for importing configuration
            dram_cmd = [
                'conda', 'run',
                '-p', env_path,
                'DRAM-setup.py', 'import_config',
                '--config_loc', file_dir
            ]
        else:
            raise ValueError(f"Unknown action '{action}' provided. Use 'download' or 'import_config'.")

        # Run the command
        subprocess.run(dram_cmd, check=True)
        print(f"[✅] DRAM {action.replace('_', ' ')} completed. Output directory: {file_dir}")
    
    except subprocess.CalledProcessError as e:
        print(f"[❌] DRAM setup failed: {e}")
    except Exception as e:
        print(f"[❌] Error: {e}")

# Example usage
if len(sys.argv) != 3:
    print("Usage: DRAM_setup.py <download|import_config> <file_directory>")
    sys.exit(1)

action = sys.argv[1]  # The action ('download' or 'import_config')
file_directory = sys.argv[2]  # The directory (either for downloading or for config file)

run_dram_setup(action, file_directory)
