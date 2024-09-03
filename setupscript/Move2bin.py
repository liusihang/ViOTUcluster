import os
import shutil
import stat

def replace_and_copy_scripts():
    # 获取当前Conda环境的目录
    conda_prefix = os.environ.get('CONDA_PREFIX')
    
    if not conda_prefix:
        print("Error: No Conda environment is currently active.")
        return

    # 确定源目录和目标目录
    source_dir = os.path.join(os.getcwd(), 'bin')
    target_dir = os.path.join(conda_prefix, 'bin')

    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # 复制和替换viralverify文件
    source_viralverify = os.path.join(source_dir, 'viralverify')
    if not os.path.exists(source_viralverify):
        print(f"Error: Source file '{source_viralverify}' does not exist.")
        return
    
    for root, dirs, files in os.walk(conda_prefix):
        for file in files:
            if file == 'viralverify':
                target_file = os.path.join(root, file)
                shutil.copy(source_viralverify, target_file)
                # 设置可执行权限
                st = os.stat(target_file)
                os.chmod(target_file, st.st_mode | stat.S_IEXEC)
                print(f"Replaced and made executable {target_file} with {source_viralverify}")

    # 复制Bin脚本文件到目标目录并设置可执行权限
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy(source_file, target_file)
            # 设置可执行权限
            st = os.stat(target_file)
            os.chmod(target_file, st.st_mode | stat.S_IEXEC)
            print(f"Copied and made executable {source_file} to {target_file}")

    source_dir = os.path.join(os.getcwd(), 'Modules')
    # 复制Source脚本文件到目标目录并设置可执行权限
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)

        if os.path.isfile(source_file):
            shutil.copy(source_file, target_file)
            # 设置可执行权限
            st = os.stat(target_file)
            os.chmod(target_file, st.st_mode | stat.S_IEXEC)
            print(f"Copied and made executable {source_file} to {target_file}")

    print("All scripts, including viralverify, have been replaced/copied and made executable in the Conda environment's bin directory.")

if __name__ == "__main__":
    replace_and_copy_scripts()