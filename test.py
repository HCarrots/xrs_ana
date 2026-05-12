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

plt.figure(figsize=(10, 6))
plt.plot(processor.eloss, processor.signals, label='Raw Signal (Mask)')


plt.title('Background Fitting (Mask Region Only)')
plt.xlabel('Energy Loss')
plt.ylabel('Signal')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7, color='gray')
plt.show()

plt.figure(figsize=(10, 6))

I0 = 1e13
sample_ro = 8.79
M = 164.93
sample_thickness = 0.1
pz = processor.xrs_energy_correction(0,8.79,sample_thickness,whichq=1)

plt.plot(processor.eloss,processor.signals[:,1])
plt.show()
plt.plot(pz,processor.signals[:,1])
plt.show()

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
    edge_pz=1.7,
    make_plots=True,
)

plt.plot(processor.eloss,processor.valence,label = "valence")
plt.plot(processor.eloss,processor.valasymmetry,label = "valasymmetry")
plt.show()

processor.get_all_valprof(1,0,make_plots= False,wait_for_input =False,interp_kind ="linear",
        q_epsilon = 1e-12,
        return_components = False,)
fit = processor.remv_alence_prof(
    whichq=[0, 2, 3],     
    eoffset=0.0,
    fit_shift=True,
    initial_scale=1.0,
    initial_shift=0.0,
    fit_start=None,
    make_plots=True,
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
    emin=5,
    emax=250,
    show_error=True,
)