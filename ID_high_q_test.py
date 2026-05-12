import os, sys
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-xrsana")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
import warnings
from scipy import interpolate, optimize
from scipy.optimize import curve_fit
from scipy.constants import Avogadro
from xrsana import math_functions, xrs_ComptonProfiles, xrs_utilities
from xrsana import xrs_extraction

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


HFCP_PATH = os.path.join(PROJECT_ROOT, "xrsana", "resources", "data", "ComptonProfiles.dat")


data_extract = xrs_extraction.edge_extraction(data,['Ho'],[1.0],{'Ho':['N4']})
data_extract.truncate(5)
data_extract.areanorm(whichq = 9.842, emin=10, emax=None)
data_extract.energycorrect(whichq = 9.842, alpha=90, densities=[8.79], samthickness=1.0)
data_extract.removeelastic(whichq=9.842, range1=[10, 20], range2=[30, 40], guess=None, stoploop=True, overwrite=True)
pz,val=data_extract.extractval_test(whichq=9.842, mirror=False, linrange1=[155,165], linrange2=None, element='Ho', edge='N4')

plt.plot(pz, val, label='Valence momentum')
plt.show()
"""
columns = xrs_extraction._as_columns(whichq=9.842)
plt.plot(data_extract.pzscale, data_extract.valencepz[:, columns], label='Valence momentum')
plt.xlabel('pz (a.u.)')
plt.ylabel('J(pz)')
plt.legend()
plt.grid(True)
plt.show()
"""
data_extract.getallvalprof(whichq=9.842, smoothgval=0.0, stoploop=False)
data_extract.remvalenceprof_test(whichq=9.842, eoffset=0.0, element='Ho', edge='N4')
data_extract.averageqs(whichq=9.842, errorweighing=True)
data_extract.save_average_Sqw("average_Sqw.dat", emin=5, emax=None, normrange=None)