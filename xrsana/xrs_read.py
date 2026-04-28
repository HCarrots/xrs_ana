import numpy as np
import os
import re
from datetime import datetime

from xrsana import xrs_utilities, xrs_scans
import h5py
# These values are used in read_Lerix class but may be useful elsewhere? LJRH
TINY = 1.e-7
MAX_FILESIZE = 100*1024*1024  # 100 Mb limit
MIN_FILESIZE = 1024 # 1 kB minimum to avoid empty files
COMMENTCHARS = '#;%*!$'
NAME_MATCH = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)*$").match
VALID_SNAME_CHARS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
VALID_NAME_CHARS = '.%s' % VALID_SNAME_CHARS
RESERVED_WORDS = ('and', 'as', 'assert', 'break', 'class', 'continue',
                'def', 'del', 'elif', 'else', 'eval', 'except', 'exec',
                'execfile', 'finally', 'for', 'from', 'global', 'if',
                'import', 'in', 'is', 'lambda', 'not', 'or', 'pass',
                'print', 'raise', 'return', 'try', 'while', 'with',
                'group', 'end', 'endwhile', 'endif', 'endfor', 'endtry',
                'enddef', 'True', 'False', 'None')

def print_citation_message():
	"""Prints plea for citing the XRStools article when using this software.

	"""
	print ('                                                                                ')
	print (' ############################# Welcome to XRStools #############################')
	print (' # If you are using this software, please cite the following work:             #')
	print (' # Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari: #')
	print (' # "Planning, performing, and analyzing X-ray Raman scattering experiments."   #')
	print (' # Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.                 #')
	print (' ###############################################################################')
	print ('                                                                                ')

class read_heps_id33:
    """
    Reader for the reduced HEPS ID33 XRS text export.

    Expected layout for a scan named ``Ho_scan1``::

        Ho_scan1/
            Ho_scan1_info.txt
            Ho_scan1_rois.txt
            Ho_scan1_data.txt

    The HEPS export already contains the reduced spectrum, so this class wraps it
    into the same attributes used by the LERIX reader: ``eloss``, ``signals``,
    ``errors``, ``E0``, ``cenom``, ``tth`` and ``scans``.
    """
    def __init__(
        self,
        exp_dir,
        scan_name=None,
        data_file=None,
        roi_file=None,
        info_file=None,
        intensity_columns=None,
        energy_column='Energy Transfer (eV)',
        use_roi_selection=True,
        q_column='q_ave',
        e0_column='center',
        q_range=None,
        modules=None,
        remove_bad_fits=True,
    ):
        import pandas as pd

        self.path = os.path.abspath(exp_dir)
        self.scan_name = scan_name or os.path.basename(os.path.normpath(self.path))
        self.info_file = info_file or os.path.join(self.path, self.scan_name + '_info.txt')
        self.roi_file = roi_file or os.path.join(self.path, self.scan_name + '_rois.txt')
        self.data_file = data_file or os.path.join(self.path, self.scan_name + '_data.txt')

        self.header_attrs = self._read_info(self.info_file)
        self.rois = pd.read_csv(self.roi_file, sep='\t')
        self.dataframe = pd.read_csv(self.data_file, sep='\t')

        self.energy_column = energy_column
        if intensity_columns is None:
            intensity_columns = [
                col for col in self.dataframe.columns
                if col not in ('Unnamed: 0', energy_column)
            ]
        self.intensity_columns = intensity_columns

        if q_range is None:
            q_range = self.header_attrs.get('q_range')
        if modules is None:
            modules = self.header_attrs.get('module_add_list')

        selected_rois = self._select_rois(
            self.rois,
            use_roi_selection=use_roi_selection,
            q_range=q_range,
            modules=modules,
            remove_bad_fits=remove_bad_fits,
        )
        if selected_rois.empty:
            selected_rois = self.rois
        self.selected_rois = selected_rois

        self.scans = {}
        self.elastic_scans = []
        self.nixs_scans = [self.scan_name + '_data.txt']
        self.wide_scans = []
        self.is_checked = []
        self.resolution = {}
        self.data = {}

        self.key = {
            'Analyzer%02d' % (idx + 1): idx
            for idx in range(len(self.intensity_columns))
        }
        self.cenom_dict = {
            key: {
                'average': {
                    'fwhm': float(np.nanmean(selected_rois['width'])) if 'width' in selected_rois else np.nan,
                    'e0': float(np.nanmean(selected_rois[e0_column])) if e0_column in selected_rois else np.nan,
                }
            }
            for key in self.key
        }

        self.eloss = self.dataframe[energy_column].to_numpy(dtype='float64')
        self.signals = self.dataframe[self.intensity_columns].to_numpy(dtype='float64')
        if self.signals.ndim == 1:
            self.signals = self.signals[:, np.newaxis]
        self.errors = np.sqrt(np.absolute(self.signals))

        e0_eV = float(np.nanmean(selected_rois[e0_column])) if e0_column in selected_rois else 0.0
        self.E0 = e0_eV / 1e3
        self.cenom = [self.E0] * len(self.intensity_columns)
        self.energy = self.E0 + self.eloss / 1e3

        q_values = selected_rois[q_column].to_numpy(dtype='float64') if q_column in selected_rois else np.array([])
        q_average = float(np.nanmean(q_values)) if q_values.size else np.nan
        self.q = q_values
        self.q_average = q_average
        self.tth = self._make_tth(q_average, len(self.intensity_columns))

        scan_attrs = dict(self.header_attrs)
        scan_attrs.update({
            'beamline': scan_attrs.get('beamline', 'HEPS ID33'),
            'selected_rois': selected_rois['crystal'].tolist() if 'crystal' in selected_rois else [],
            'q_average': q_average,
        })

        onescan = xrs_scans.scan(
            [],
            1,
            self.energy,
            np.ones_like(self.eloss),
            self.signals,
            scan_attrs,
            self.dataframe,
            'nixs',
        )
        onescan.eloss = self.eloss
        onescan.signals = self.signals
        onescan.errors = self.errors
        self.scans[self.scan_name] = onescan

    def _read_info(self, filename):
        attrs = {}
        if not os.path.isfile(filename):
            return attrs
        import ast

        with open(filename, 'r') as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or '=' not in line:
                    continue
                key, value = [part.strip() for part in line.split('=', 1)]
                try:
                    attrs[key] = ast.literal_eval(value)
                except Exception:
                    try:
                        attrs[key] = float(value)
                    except ValueError:
                        attrs[key] = value
        return attrs

    def _select_rois(self, rois, use_roi_selection=True, q_range=None, modules=None, remove_bad_fits=True):
        selected = rois.copy()
        if not use_roi_selection:
            return selected

        if modules and 'crystal' in selected:
            modules = [modules] if isinstance(modules, str) else list(modules)
            prefixes = tuple(module + '-' for module in modules)
            selected = selected[selected['crystal'].astype(str).str.startswith(prefixes)]

        remove_list = set(self.header_attrs.get('ROIremoveList', []))
        if remove_bad_fits:
            remove_list.update(self.header_attrs.get('BadFitList', []))
        if remove_list and 'crystal' in selected:
            selected = selected[~selected['crystal'].isin(remove_list)]

        if q_range and 'q_ave' in selected:
            q_min, q_max = q_range
            selected = selected[(selected['q_ave'] >= q_min) & (selected['q_ave'] <= q_max)]

        return selected

    def _make_tth(self, q_average, number_of_columns):
        if not np.isfinite(q_average) or self.E0 <= 0.0:
            return [np.nan] * number_of_columns

        k = 0.5067731 * self.E0
        cos_tth = 1.0 - (q_average ** 2) / (2.0 * k ** 2)
        tth = np.degrees(np.arccos(np.clip(cos_tth, -1.0, 1.0)))
        return [float(tth)] * number_of_columns

    def load_elastics(self, *args, **kwargs):
        print('HEPS ID33 reduced data: elastic calibration is already stored in the ROI table.')
        return self

    def load_nixs(self, *args, **kwargs):
        print('HEPS ID33 reduced data loaded: {}'.format(self.data_file))
        return self

    def load_wides(self, *args, **kwargs):
        print('HEPS ID33 reduced data has no separate wide scan in this export.')
        return self

    def join_nixs_wide(self, *args, **kwargs):
        print('HEPS ID33 reduced data has no separate wide scan to join.')
        return self

    def update_cenom(self, analyzers='all'):
        print("{} {}".format("E0 was found to be (keV): ", self.E0))
        return self.E0

    def scan_info(self, file):
        return (1, self.scan_name, 'nixs', os.path.basename(file))

    def save_H5(self, H5name='HEPS_ID33_data.H5'):
        H5path = H5name if os.path.isabs(H5name) else os.path.join(self.path, H5name)
        if not os.path.isdir(os.path.dirname(H5path)):
            print('H5 path directory does not exist!')
            return

        with h5py.File(H5path, 'a') as H5file:
            group_name = os.path.basename(H5path)
            if group_name in H5file:
                del H5file[group_name]
            g = H5file.create_group(group_name)
            H5_xrs = g.create_group('XRS')
            h5group = H5_xrs.create_group(self.scan_name)
            h5group.create_dataset('energy', data=self.energy)
            h5group.create_dataset('signals', data=self.signals)
            h5group.create_dataset('eloss', data=self.eloss)
            h5group.create_dataset('errors', data=self.errors)

            g.create_dataset('energy', data=self.energy)
            g.create_dataset('signals', data=self.signals)
            g.create_dataset('eloss', data=self.eloss)
            g.create_dataset('errors', data=self.errors)
            g.create_dataset('tth', data=self.tth)
