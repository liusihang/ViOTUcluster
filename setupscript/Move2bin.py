import os
import shutil
import stat

def copy_scripts_to_conda_bin():
    # 获取当前Conda环境的目录
    conda_prefix = os.environ.get('CONDA_PREFIX')
    
    if not conda_prefix:
        print("Error: No Conda environment is currently active.")
        return

    # 确定源和目标目录
    source_dir = os.path.join(os.getcwd(), 'bin')
    target_dir = os.path.join(conda_prefix, 'bin')

    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # 复制所有文件到目标目录，并设置可执行权限
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy(source_file, target_file)
            # 设置可执行权限
            st = os.stat(target_file)
            os.chmod(target_file, st.st_mode | stat.S_IEXEC)
            print(f"Copied and made executable {source_file} to {target_file}")

    print("All scripts have been copied to the Conda environment's bin directory and made executable.")

if __name__ == "__main__":
    copy_scripts_to_conda_bin()