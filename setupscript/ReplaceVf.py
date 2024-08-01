import os
import shutil

def replace_viralverify_files():
    # 获取当前Conda环境的目录
    conda_prefix = os.environ.get('CONDA_PREFIX')
    
    if not conda_prefix:
        print("Error: No Conda environment is currently active.")
        return

    # 确定源文件路径
    source_file = os.path.join(os.getcwd(), 'bin', 'viralverify')

    if not os.path.exists(source_file):
        print(f"Error: Source file '{source_file}' does not exist.")
        return

    # 遍历Conda环境目录，查找并替换所有viralverify文件
    for root, dirs, files in os.walk(conda_prefix):
        for file in files:
            if file == 'viralverify':
                target_file = os.path.join(root, file)
                shutil.copy(source_file, target_file)
                print(f"Replaced {target_file} with {source_file}")

    print("All viralverify files have been replaced in the Conda environment.")

if __name__ == "__main__":
    replace_viralverify_files()