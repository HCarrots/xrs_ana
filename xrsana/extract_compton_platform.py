#!/usr/bin/env python
"""
Extract Compton platform data from ID33 experimental files using the xrstools framework.

This script reads two experimental data files (ID33_HL.dat and ID33_VD.dat)
from the sandbox/data directory and uses the XRStools data processing framework
to extract the Compton platform data for each analyzer chamber.

Workflow:
  1. Read the experimental data (energy loss vs intensity)
  2. Create a data container compatible with XRStools edge_extraction
  3. Compute Hartree-Fock (HF) Compton profiles for the sample
  4. Fit and remove the elastic peak (Pearson7) + linear background
  5. The remaining signal is the extracted Compton platform data
  6. Save results and generate comparison plots

Usage:
    python extract_compton_platform.py
"""

import numpy as np
import os
import sys

import matplotlib

import matplotlib.pyplot as plt

from xrsana import xrs_extraction, xrs_ComptonProfiles

'''
# =====================================================================
# Configuration Parameters (modify these for your experiment)
# =====================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

HL_FILE = os.path.join(DATA_DIR, "ID33_HL.dat")
VD_FILE = os.path.join(DATA_DIR, "ID33_VD.dat")

# Sample composition
#   Silicon:   FORMULAS=['Si'],   STOICH_WEIGHTS=[1.0], EDGES={'Si':['K']}
#   Diamond:   FORMULAS=['C'],    STOICH_WEIGHTS=[1.0], EDGES={'C':['K']}
#   Water:     FORMULAS=['H2O'],  STOICH_WEIGHTS=[1.0], EDGES={'O':['K']}
#   LiCl:      FORMULAS=['LiCl'], STOICH_WEIGHTS=[1.0], EDGES={'Li':['K']}
#   Silicon L23: FORMULAS=['Si'], STOICH_WEIGHTS=[1.0], EDGES={'Si':['L23']}
FORMULAS = ["B"]
STOICH_WEIGHTS = [1.0]
EDGES = {"B": ["K"]}

# Incident energy (keV) - typical value for ESRF ID33
E0 = 9.7

# Scattering angles (degrees) for each analyzer chamber
# HL = Horizontal Left, VD = Vertical Down
TTH_HL = 57.58
TTH_VD = -18.0

# Pre-normalization energy range (eV)
# Experimental signals are scaled so their integral in this range
# matches the HF Compton profile integral
PRENORM_RANGE = [5.0, 1000.0]

# Fit ranges for background removal (eV)
# range1: pre-edge region (Compton tail before the edge)
# range2: post-edge region (after the edge)
# For C K-edge: the edge is at ~284 eV
RANGE1_HL = [50.0, 250.0]
RANGE2_HL = [300.0, 600.0]
RANGE1_VD = [50.0, 250.0]
RANGE2_VD = [300.0, 600.0]

# Energy range for saving final results (eV)
EMIN_SAVE = 0.0
EMAX_SAVE = None  # None means save the full range

# Normalization range for the saved S(q,w)
NORM_RANGE = (
    None  # e.g. [520.0, 600.0] for area-normalization, None for no normalization
)

# Use logarithmic scale for plots
USE_LOG_PLOT = False

'''
# =====================================================================
# Data Container
# =====================================================================


class ExpDataContainer:
    """
    Data container compatible with XRStools edge_extraction.

    Provides all attributes expected by edge_extraction and HF_dataset:
      eloss   : 1D array of energy loss values (eV)
      signals : 2D array of signals, shape (n_eloss, n_analyzers)
      errors  : 2D array of Poisson errors, shape (n_eloss, n_analyzers)
      E0      : incident energy (keV)
      tth     : list of scattering angles (degrees) per analyzer
      cenom   : list of center-of-mass values (keV) per analyzer
    """

    def __init__(self, eloss, signals, errors, E0, tth):
        self.eloss = np.asarray(eloss, dtype=np.float64)
        # Ensure signals and errors are 2D arrays for compatibility.
        signals = np.asarray(signals, dtype=np.float64)
        errors = np.asarray(errors, dtype=np.float64)
        self.signals = signals[:, None] if signals.ndim == 1 else signals
        self.errors = errors[:, None] if errors.ndim == 1 else errors
        self.E0 = float(E0)
        self.tth = list(tth)
        self.cenom = [float(E0)] * len(tth)



def read_experimental_data(filename):
    """
    Read an experimental .dat file with tab-separated columns:
      Column 1: Energy loss (eV)
      Column 2: Intensity (a.u.)

    Returns:
        eloss     : 1D array of energy loss values (eV)
        intensity : 1D array of intensity values
    """
    data = np.loadtxt(filename, skiprows=1, delimiter="\t")
    return data[:, 0], data[:, 1]


def get_hf_core(
    eloss,
    E0,
    tth,
    formulas,
    stoich_weights,
    edges,
    element=None,
    edge=None,
    columns=None,
    HFcore_shift=0.0,
    normalize=False,
    norm_range=None,
    norm_area=1.0,
):
    """
    Return only the Hartree-Fock core profile on the supplied eloss grid.

    This is the lightweight path for notebooks when the HF core is needed
    without fitting or subtracting any experimental background. Set
    normalize=True to scale the integral over norm_range to norm_area.
    """
    tth = [tth] if np.isscalar(tth) else list(tth)
    dummy = np.zeros((len(eloss), len(tth)))
    exp_data = ExpDataContainer(eloss, dummy, np.ones_like(dummy), E0, tth)
    hf_data = xrs_extraction.HF_dataset(exp_data, formulas, stoich_weights, edges)

    if element is None:
        element = next(iter(edges))
    if edge is None:
        edge = edges[element][0]
    if columns is None:
        columns = range(len(tth))

    core = np.mean(hf_data.C_edges[element][edge][:, list(columns)], axis=1)
    if HFcore_shift:
        core = np.interp(hf_data.eloss, hf_data.eloss + HFcore_shift, core)
    if normalize:
        if norm_range is None:
            inds = np.arange(len(hf_data.eloss))
        else:
            inds = np.where(
                np.logical_and(
                    hf_data.eloss >= norm_range[0],
                    hf_data.eloss <= norm_range[1],
                )
            )[0]
        area = np.trapezoid(core[inds], hf_data.eloss[inds])
        if len(inds) < 2 or not np.isfinite(area) or np.isclose(area, 0.0):
            raise ValueError("Cannot normalize HF core: selected integral is zero or invalid.")
        core = core * (norm_area / area)
    return core


def extract_compton_platform(
    eloss,
    intensity,
    E0,
    tth,
    formulas,
    stoich_weights,
    edges,
    prenorm_range,
    range1,
    range2,
):
    """
    Extract Compton platform data from experimental data using XRStools.

    Steps:
      1. Create an ExpDataContainer with the experimental data
      2. Initialize edge_extraction (computes HF Compton profiles)
      3. Average over the analyzer
      4. Remove the elastic peak and linear background (Pearson7 fit
         guided by the HF core Compton profile)
      5. Save the extracted Compton platform data and HF profiles

    Returns:
        extraction_obj : the edge_extraction object containing all results
    """
    n_pts = len(eloss)

    signals = intensity.reshape(-1, 1)
    errors = np.sqrt(np.maximum(np.abs(intensity), 1.0)).reshape(-1, 1)

    exp_data = ExpDataContainer(eloss, signals, errors, E0, [tth])

    element = list(edges.keys())[0]
    edge = list(edges.values())[0][0]

    print(
        f"    Computing HF Compton profiles for {formulas} "
        f"(element={element}, edge={edge})..."
    )

    # Pre-normalize manually to avoid issues with edge_extraction's prenormrange
    print("    Pre-normalizing experimental data...")
    hf = xrs_ComptonProfiles.HFProfile(
        formulas,
        stoich_weights,
        "/home/hushiqi/work/xrstools/XRStools/resources/data/ComptonProfiles.dat",
    )
    hf.get_elossProfiles(E0, [tth])
    hf_J_interp = np.interp(eloss, hf.eloss, hf.J_total[:, 0])
    HFnorm = np.trapezoid(hf_J_interp, eloss)
    prenorm_inds = np.where(
        np.logical_and(eloss >= prenorm_range[0], eloss <= prenorm_range[1])
    )[0]
    EXPnorm = np.trapezoid(intensity[prenorm_inds], eloss[prenorm_inds])
    norm_factor = HFnorm / EXPnorm
    signals = intensity * norm_factor
    errors = errors * norm_factor
    print(f"      Normalization factor: {norm_factor:.8f}")

    # Create edge_extraction object with pre-normalized data
    exp_data = ExpDataContainer(eloss, signals, errors, E0, [tth])
    extraction_obj = xrs_extraction.edge_extraction(
        exp_data, formulas, stoich_weights, edges, prenormrange=None
    )

    print("    Averaging analyzer signal...")
    extraction_obj.analyzerAverage([0], errorweighing=False)
    hf.get_elossProfiles(E0, [tth])
    hf_J_interp = np.interp(eloss, hf.eloss, hf.J_total[:, 0])
    HFnorm = np.trapezoid(hf_J_interp, eloss)
    prenorm_inds = np.where(
        np.logical_and(eloss >= prenorm_range[0], eloss <= prenorm_range[1])
    )[0]
    EXPnorm = np.trapezoid(intensity[prenorm_inds], eloss[prenorm_inds])
    norm_factor = HFnorm / EXPnorm
    signals = intensity * norm_factor
    errors = errors * norm_factor
    print(f"      Normalization factor: {norm_factor:.8f}")

    # Create edge_extraction object with pre-normalized data
    exp_data = ExpDataContainer(eloss, signals, errors, E0, [tth])
    extraction_obj = xrs_extraction.edge_extraction(
        exp_data, formulas, stoich_weights, edges, prenormrange=None
    )

    print("    Averaging analyzer signal...")
    # For a single analyzer, analyzerAverage expects a column index
    extraction_obj.analyzerAverage([0], errorweighing=False)

    print(f"    Fitting background (range1={range1}, range2={range2})...")
    try:
        extraction_obj.removeCorePearsonAv(
            element,
            edge,
            range1,
            range2,
            weights=[2, 1],
            HFcore_shift=0.0,
            show_plots=False,
        )
    except Exception as e:
        print(f"    WARNING: removeCorePearsonAv failed ({e})")
        print("    Falling back to removeCorePearsonAv_new...")
        extraction_obj.removeCorePearsonAv_new(
            element, edge, range1, range2, HFcore_shift=0.0, reg_lam=10
        )
   
    hf_data = np.column_stack(
        [
            extraction_obj.eloss,
            extraction_obj.HF_dataset.J_total[:, 0],
            extraction_obj.HF_dataset.C_total[:, 0],
            extraction_obj.HF_dataset.V_total[:, 0],
            extraction_obj.HF_dataset.q_vals[:, 0],
        ]
    )
    return extraction_obj
