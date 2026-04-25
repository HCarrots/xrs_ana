#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
XRStools 石墨 X 射线拉曼散射数据处理脚本
功能：数据读取 → 能量校准 → 光谱合并 → 边缘提取 → 结果保存
"""

# ===================== 1. 基础库导入与环境配置 =====================
import sys
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt

# 屏蔽非关键警告（如数值溢出、除零等，不影响最终结果）
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# 设置 matplotlib 后端（非交互式，适合脚本运行；如需交互式可改为 'Qt5Agg'）
plt.switch_backend('agg')
# 设置绘图样式
plt.rcParams['figure.figsize'] = (10, 5)
plt.rcParams['font.size'] = 12
plt.rcParams['savefig.dpi'] = 150

# ===================== 2. 路径与参数配置（请根据实际情况修改） =====================
# 项目根目录（xrsana 文件夹的上一级）
PROJECT_ROOT = "/home/hushiqi/work/xrs_ana"
# 实验数据目录
DATA_PATH = "/home/hushiqi/work/forked_xrstools/XRStools/LERIX-Guide/graphite/"
# 结果保存目录（自动创建）
RESULT_PATH = os.path.join(DATA_PATH, "results")
os.makedirs(RESULT_PATH, exist_ok=True)

# 选定的分析器编号
ANALYZERS = [2, 5, 8, 9, 10]
# 用于平均的分析器
AVG_ANALYZERS = [2, 5, 8]

# ===================== 3. 导入 XRStools 处理模块 =====================
# 添加项目路径到 Python 搜索路径
sys.path.append(PROJECT_ROOT)
try:
    from xrsana import xrs_read, xrs_extraction
    print("✅ XRStools 模块导入成功")
except ImportError as e:
    print(f"❌ XRStools 模块导入失败: {e}")
    print("请检查 PROJECT_ROOT 路径是否正确")
    sys.exit(1)

# ===================== 4. 数据读取与校准 =====================
print("\n==================== 开始数据处理 ====================")

# 4.1 初始化 LERIX 数据读取器
print("📂 正在读取实验数据...")
cnt = xrs_read.read_lerix(
    exp_dir=DATA_PATH,
    elastic_name='elastic',
    nixs_name='NIXS',
    wide_name='wide'
)

# 4.2 加载弹性峰并校准
print("⚡ 正在加载弹性峰并校准能量...")
cnt.load_elastics(analyzers='all')
cnt.load_nixs()
cnt.load_wides()

# 4.3 合并 NIXS 谱和宽谱
print("🔗 正在合并 NIXS 谱与宽谱...")
cnt.join_nixs_wide(scaling='auto')

# ===================== 5. 数据可视化 =====================
print("📊 正在生成光谱图...")

# 绘制单个分析器光谱
plt.figure()
plt.plot(cnt.eloss, cnt.signals[:, 4], label=f'Analyzer {ANALYZERS[4]}')
plt.xlabel('Energy Loss (eV)')
plt.ylabel('Intensity (a.u.)')
plt.title('Single Analyzer Spectrum')
plt.legend()
plt.grid(alpha=0.3)
plt.savefig(os.path.join(RESULT_PATH, 'single_analyzer_spectrum.png'))
plt.close()

# 绘制多个分析器平均前光谱
plt.figure()
for ana in AVG_ANALYZERS:
    plt.plot(cnt.eloss, cnt.signals[:, ana], label=f'Analyzer {ana}', alpha=0.7)
plt.xlabel('Energy Loss (eV)')
plt.ylabel('Intensity (a.u.)')
plt.title('Multiple Analyzer Spectra')
plt.legend()
plt.grid(alpha=0.3)
plt.savefig(os.path.join(RESULT_PATH, 'multiple_analyzer_spectra.png'))
plt.close()

# ===================== 6. 保存中间数据 =====================
print("💾 正在保存中间 HDF5 数据...")
cnt.save_H5(H5name=os.path.join(RESULT_PATH, 'PS_CNT_369_TAKE2.h5'))

# 更新能量中心
print("🔄 正在更新能量中心...")
cnt.update_cenom(analyzers='all')

# ===================== 7. 边缘提取与物理分析 =====================
print("🔬 正在进行碳 K 边边缘提取...")

# 初始化边缘提取器
cnt_ex = xrs_extraction.edge_extraction(cnt, ['C'], [1.0], {'C': ['K']})

# 分析器信号平均
print("📈 正在平均分析器信号...")
cnt_ex.analyzerAverage(AVG_ANALYZERS, errorweighing=False)

# 扣除芯层电子贡献（Pearson 背景）
print("🧹 正在扣除芯层电子背景...")
cnt_ex.removeCorePearsonAv(
    'C', 'K',
    [150, 280.0],   # 预边拟合范围
    [330.0, 450.0],  #  post边拟合范围
    weights=[2, 1],
    HFcore_shift=5.6,
    scaling=5.9
)

# ===================== 8. 保存最终结果（已修复属性名） =====================
print("🎉 正在保存最终处理结果...")

# 保存动态结构因子 S(q,ω)（这是核心功能，必须保留）
output_file = os.path.join(RESULT_PATH, 'PSCNT_data[369]_extract_BETTERSUBTRACT.dat')
cnt_ex.save_average_Sqw(output_file, emin=275, emax=340)

# 绘制最终 S(q,ω) 光谱（已修复属性名，并增加容错处理）
print("📊 正在生成最终光谱图...")
try:
    # 尝试获取能量轴（优先尝试 elossav，若失败则使用 eloss）
    if hasattr(cnt_ex, 'elossav'):
        x_axis = cnt_ex.elossav
    elif hasattr(cnt_ex, 'eloss'):
        x_axis = cnt_ex.eloss
    else:
        # 如果都没有，使用 sqwav 的索引作为 x 轴
        x_axis = np.arange(len(cnt_ex.sqwav))
        print("⚠️  未找到明确的能量轴，使用数组索引代替")

    plt.figure()
    plt.plot(x_axis, cnt_ex.sqwav, 'k-', linewidth=1.5)
    plt.xlabel('Energy Loss (eV)')
    plt.ylabel('S(q,ω) (a.u.)')
    plt.title('Final XRS Spectrum - Carbon K Edge')
    
    # 仅在有明确能量轴时设置 x 范围
    if hasattr(cnt_ex, 'elossav') or hasattr(cnt_ex, 'eloss'):
        plt.xlim(275, 340)
        
    plt.grid(alpha=0.3)
    plt.savefig(os.path.join(RESULT_PATH, 'final_XRS_spectrum.png'))
    plt.close()
    print("✅ 最终光谱图已保存")
    
except Exception as e:
    print(f"⚠️  绘图过程中出现错误，已跳过: {e}")
    print("   但核心数据文件已成功保存！")
# ===================== 9. 处理完成总结 =====================
print("\n==================== 数据处理完成 ====================")
print(f"📁 所有结果已保存至: {RESULT_PATH}")
print("   - 中间数据: PS_CNT_369_TAKE2.h5")
print("   - 最终光谱: PSCNT_data[369]_extract_BETTERSUBTRACT.dat")
print("   - 可视化图: *.png")
print("\n✅ 石墨 XRS 数据处理流程结束！")