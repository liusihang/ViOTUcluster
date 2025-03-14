import os
import shutil
import stat

def replace_and_copy_scripts():
    # Get the current Conda environment directory
    conda_prefix = os.environ.get('CONDA_PREFIX')
    
    if not conda_prefix:
        print("Error: No Conda environment is currently active.")
        return

    # Define source and target directories
    source_dir = os.path.join(os.getcwd(), 'bin')
    target_dir = os.path.join(conda_prefix, 'bin')

    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Copy and replace viralverify file
    source_viralverify = os.path.join(source_dir, 'viralverify')
    if not os.path.exists(source_viralverify):
        print(f"Error: Source file '{source_viralverify}' does not exist.")
        return
    
    for root, dirs, files in os.walk(conda_prefix):
        for file in files:
            if file == 'viralverify':
                target_file = os.path.join(root, file)
                shutil.copy(source_viralverify, target_file)
                # Set executable permissions
                st = os.stat(target_file)
                os.chmod(target_file, st.st_mode | stat.S_IEXEC)
                print(f"Replaced and made executable {target_file} with {source_viralverify}")

    # Copy bin scripts to target directory and set executable permissions
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy(source_file, target_file)
            # Set executable permissions
            st = os.stat(target_file)
            os.chmod(target_file, st.st_mode | stat.S_IEXEC)
            print(f"Copied and made executable {source_file} to {target_file}")

    source_dir = os.path.join(os.getcwd(), 'Modules')
    # Copy module scripts to target directory and set executable permissions
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy(source_file, target_file)
            # Set executable permissions
            st = os.stat(target_file)
            os.chmod(target_file, st.st_mode | stat.S_IEXEC)
            print(f"Copied and made executable {source_file} to {target_file}")

    print("All scripts, including viralverify, have been replaced/copied and made executable in the Conda environment's bin directory.")

if __name__ == "__main__":
    replace_and_copy_scripts()