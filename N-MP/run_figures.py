"""
运行此脚本以生成第8章所有图表
Run this script to generate all figures for Chapter 8
"""
import os
import subprocess
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(script_dir, 'figures'), exist_ok=True)

scripts = [
    'fig8_1_radar.py',
    'fig8_2_biodistribution.py',
    'fig8_3_barriers.py',
    'fig8_4_pbpk.py',
    'fig8_5_timeline.py',
]

print("=== 第8章图表生成脚本 ===")
print(f"输出目录: {os.path.join(script_dir, 'figures')}\n")

for script in scripts:
    path = os.path.join(script_dir, 'figures', script)
    print(f"正在生成: {script} ...", end=' ')
    try:
        result = subprocess.run([sys.executable, path],
                                capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print("✓ 完成")
        else:
            print(f"✗ 错误:\n{result.stderr}")
    except subprocess.TimeoutExpired:
        print("✗ 超时")
    except FileNotFoundError:
        print(f"✗ 文件未找到: {path}")

print("\n=== 图表生成完毕 ===")
print("PNG和PDF格式文件已保存至 figures/ 目录")
