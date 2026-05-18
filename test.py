import os, sys
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-xrsana")
import matplotlib.pyplot as plt
import warnings
PROJECT_ROOT = "/home/hushiqi/work/xrs_ana"
SCAN_NAME = "Ho"
DATA_PATH = os.path.join(PROJECT_ROOT, "ex_space", "analysis", "data", SCAN_NAME)
RESULT_PATH = os.path.join(PROJECT_ROOT, "ex_space", "analysis", "result")
RESULT_DIR = os.path.join(RESULT_PATH, SCAN_NAME)
SQW_PATH = os.path.join(RESULT_DIR, "result.dat")

os.makedirs(RESULT_DIR, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-xrsana")
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, PROJECT_ROOT)

from xrsana import xrs_read
print("Loading HEPS ID33 reduced data :", DATA_PATH)
data = xrs_read.read_heps_id33(DATA_PATH, q_range=(9.2, 10))
print("E0 (keV)        :", data.E0)
print("tth (deg)       :", data.tth)
print("q values        :", data.q)
print("analyzer key    :", data.key)
print("selected ROIs   :", len(data.selected_rois))
from xrsana import xrs_process

processor = xrs_process.XRSProcess(data, ["Ho"], [1], {"Ho": ["N4"]})
processor.xrs_remove_elastic()
processor.xrs_remove_stray_background()

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    processor.eloss,
    processor.signals,
    color='#1f77b4',
    linewidth=1.2,
    label='Raw Signal (Mask Region)'
)

ax.set_title('Inelastic peak full spectrum', pad=15)
ax.set_xlabel('Energy Loss [eV]')
ax.set_ylabel('Intensity [arb. units]')

ax.grid(True, linestyle='--', alpha=0.5, color='lightgray', zorder=0)
ax.set_axisbelow(True)
#ax.legend(frameon=False)

fig.savefig('/home/hushiqi/report/2026-5-5/pic/pictemp/mask.pdf')
fig.savefig('/home/hushiqi/report/2026-5-5/pic/pictemp/mask.png')
plt.close(fig)


I0 = 1e13
sample_ro = 8.79
M = 164.93
sample_thickness = 0.1
pz = processor.xrs_energy_correction(0,8.79,sample_thickness,whichq=1)


fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    processor.eloss,
    processor.signals[:, 1],       
    color='#1f77b4',
    linewidth=1.2,
    label=f'q-channel 1 (tth={processor.tth[1]:.1f}°)'
)


ax.set_title('Energy Loss vs Intensity — q-channel 1', pad=15)
ax.set_xlabel('Energy Loss [eV]')
ax.set_ylabel('Intensity [arb. units]')


ax.legend(frameon=False) 
ax.grid(True, linestyle='--', alpha=0.5, color='lightgray', zorder=0)
ax.set_axisbelow(True)    

fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/q-channel-1.pdf')
fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/q-channel-1.png')
plt.close(fig)





plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

fig, ax = plt.subplots(figsize=(8, 5))


ax.plot(
    pz[:, 0],
    processor.signals[:, 1],
    color='#1f77b4',
    linewidth=1.2,
    label=f'q-channel 1 (tth={processor.tth[1]:.1f}°)'
)


ax.set_title('q-channel 1 in Momentum Space Signal', pad=15)
ax.set_xlabel('pz [a.u.]')
ax.set_ylabel('Intensity [arb. units]')


ax.legend(frameon=False)
ax.grid(True, linestyle='--', alpha=0.5, color='lightgray', zorder=0)
ax.set_axisbelow(True)


fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/Momentum.pdf')
fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/Momentum.png')

plt.close(fig)





processor.xrs_remove_poly_core(
        whichq = 1,
        polyregion = (200,600),
        coreregion = (150,170),
        weights=(1.0, 1.0),
        polyorder=2,
        scale=6,
        hfcoreshift=0.0,
        plot=True,
        ewindow=100.0,
        save_result=True,)


processor.xrs_remove_poly_core_2(
        whichq = 1,
        polyregion = [[100,150],[600,800]],
        coreregion = (150.2,170),
        weights=(1.0, 1.0),
        polyorder=2,
        scale=15,
        hfcoreshift=-3,
        plot=True,
        ewindow=100.0,
        save_result=True,)

processor.extractval(
    whichq=1,
    mirror=False,
    linrange1=(10, 20),
    linrange2=(80, 100),
    edge_pz=0.8,
    make_plots=False,
)



plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

fig, ax = plt.subplots(figsize=(8, 5))

# 绘制曲线（统一线条宽度，区分颜色提升可读性）
ax.plot(
    processor.eloss,
    processor.valence,
    color='#1f77b4',  # 标准Matplotlib蓝色
    linewidth=1.2,
    label="valence"
)
ax.plot(
    processor.eloss,
    processor.valasymmetry,
    color='#ff7f0e',  # 标准Matplotlib橙色
    linewidth=1.2,
    label="valence asymmetry"
)

# 标题与坐标轴（学术规范：标题增加间距避免拥挤）
ax.set_title('Extracted Valence & Asymmetry Profiles', pad=15)
ax.set_xlabel('Energy Loss [eV]')
ax.set_ylabel('Intensity [arb. units]')

# 图例与网格（网格置于底层，不遮挡数据）
#ax.legend(frameon=False)  # 无边框图例，期刊标准样式
ax.grid(True, linestyle='--', alpha=0.5, color='lightgray', zorder=0)
ax.set_axisbelow(True)

# 紧凑布局 + 显示 + 释放内存
plt.tight_layout()
# 按图名保存，避免覆盖
fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/valence_asymmetry_profiles.pdf')
fig.savefig('/home/hushiqi/work/xrs_ana/pictemp/valence_asymmetry_profiles.png')
plt.close(fig)




processor.get_all_valprof(
    1,
    0,
    make_plots= False,
    wait_for_input =False,
    interp_kind ="linear",
    q_epsilon = 1e-12,
    return_components = False)

fit = processor.remv_alence_prof(
    whichq=[0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],     
    eoffset=0.0,
    fit_shift=True,
    initial_scale=1.0,
    initial_shift=0.0,
    fit_start=None,
    make_plots=False,
    wait_for_input=True,
    update_valence=True,
    interp_kind="linear",
    return_fit=True,
)


result = processor.averageqs(
    whichq=[0, 2, 3],      
    error_weighting=True,
    min_error=1.0,
    legacy_unweighted_sum=False,
    return_result=True,
)

processor.plotsqwav(
    emin=140,
    emax=200,
    show_error=True,
    savepath="/home/hushiqi/work/xrs_ana/pictemp/sqw.png",

)
