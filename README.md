# xrs_ana
This project is derived from the [xrstools](https://github.com/christophsahle/XRStools) projec in ESRF. It is compatible with the ID33 beamline of the High Energy Photon Source (HEPS) and has been migrated and adapted for Python 3.15.
## Development Plan
1. Fix the functions related to xrs_prediction and support two usage modes: command-line interface (CLI) and Jupyter Notebook (ipynb) environment.
2. Performing function fix in the future.

## Browse Reduced Crystal Data

Reduced HEPS ID33 text exports can be browsed in a local web UI:

```bash
xrsana-browse-data /home/hushiqi/work/xrs_ana/ex_space/analysis/data/Ho
```

During development, the module can also be run directly:

```bash
python -m xrsana.data_browser /home/hushiqi/work/xrs_ana/ex_space/analysis/data/Ho
```

The browser renders one labeled spectrum image per crystal and shows the crystal
number, momentum transfer magnitude (`q_ave`), and scattering angle (`2theta`).
It also supports filtering by crystal name and q range.
