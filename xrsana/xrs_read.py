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

    The reader also accepts ``*_all_data.txt`` files, and it discovers the scan
    prefix from the standard ``*_info.txt``, ``*_rois.txt`` and data-file triplet.
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
        auto_load=True,
    ):
        self.path = os.path.abspath(exp_dir)
        self.requested_scan_name = scan_name
        self.user_info_file = info_file
        self.user_roi_file = roi_file
        self.user_data_file = data_file
        self.user_intensity_columns = intensity_columns
        self.energy_column = energy_column
        self.use_roi_selection = use_roi_selection
        self.q_column = q_column
        self.e0_column = e0_column
        self.q_range = q_range
        self.modules = modules
        self.remove_bad_fits = remove_bad_fits

        self.nixs_name = scan_name
        self.wide_name, self.wide_scans = 'wide', []
        self.elastic_name, self.elastic_scans = 'elastic', []
        self.scan_files = self.sort_dir()

        self.scans = {}
        self.nixs_scans = [group['data_file'] for group in self.scan_files]
        self.scan_name = (
            self.scan_files[0]['scan_name']
            if len(self.scan_files) == 1
            else os.path.basename(os.path.normpath(self.path))
        )
        self.sample_name = self.scan_name
        self.scannumbers = []
        self.is_checked = []
        self.resolution = {}
        self.header_attrs = {}
        self.rois = None
        self.selected_rois = None
        self.dataframe = None
        self.info_file = None
        self.roi_file = None
        self.data_file = None
        self.intensity_columns = []
        self.data = {}
        self.key = {}
        self.cenom_dict = {}
        self.eloss = np.array([])
        self.signals = np.array([])
        self.errors = np.array([])
        self.E0 = np.nan
        self.cenom = []
        self.energy = np.array([])
        self.q = np.array([])
        self.q_average = np.nan
        self.tth = []

        if not self.isValidDir(self.path):
            return
        if auto_load:
            self.load_nixs()

    def sort_dir(self, path=None):
        """Find standard HEPS ID33 reduced data triplets."""
        if path is None:
            path = self.path
        if self.user_data_file:
            data_file = os.path.abspath(self.user_data_file)
            return [self._make_scan_group(data_file, self.user_info_file, self.user_roi_file)]

        groups = self._find_scan_groups(path)
        if not groups:
            for entry in sorted(os.listdir(path)) if os.path.isdir(path) else []:
                child = os.path.join(path, entry)
                if os.path.isdir(child):
                    groups.extend(self._find_scan_groups(child))

        if self.requested_scan_name:
            requested = str(self.requested_scan_name).lower()
            groups = [
                group for group in groups
                if group['scan_name'].lower() == requested
                or os.path.basename(group['scan_dir']).lower() == requested
            ]
        return sorted(groups, key=lambda group: group['data_file'])

    def _find_scan_groups(self, path):
        if not os.path.isdir(path):
            return []
        groups = []
        for filename in os.listdir(path):
            if filename.endswith('_all_data.txt') or filename.endswith('_data.txt'):
                groups.append(self._make_scan_group(os.path.join(path, filename)))
        return [
            group for group in groups
            if os.path.isfile(group['info_file']) and os.path.isfile(group['roi_file'])
        ]

    def _make_scan_group(self, data_file, info_file=None, roi_file=None):
        data_file = os.path.abspath(data_file)
        scan_dir = os.path.dirname(data_file)
        basename = os.path.basename(data_file)
        if basename.endswith('_all_data.txt'):
            scan_name = basename[:-len('_all_data.txt')]
        elif basename.endswith('_data.txt'):
            scan_name = basename[:-len('_data.txt')]
        else:
            scan_name = os.path.splitext(basename)[0]
        return {
            'scan_dir': scan_dir,
            'scan_name': scan_name,
            'data_file': data_file,
            'info_file': os.path.abspath(info_file) if info_file else os.path.join(scan_dir, scan_name + '_info.txt'),
            'roi_file': os.path.abspath(roi_file) if roi_file else os.path.join(scan_dir, scan_name + '_rois.txt'),
        }

    def isValidDir(self, dir):
        if not os.path.isdir(dir) and not self.user_data_file:
            print('Check the directory you have supplied')
            return False
        if not self.nixs_scans:
            print('No HEPS ID33 reduced data triplet found. Expected *_info.txt, *_rois.txt and *_data.txt or *_all_data.txt.')
            return False
        return True

    def scan_info(self, file):
        group = self._group_for_data_file(file)
        try:
            scan_number = self.nixs_scans.index(group['data_file']) + 1
        except ValueError:
            scan_number = 1
        return (scan_number, group['scan_name'], 'nixs', os.path.basename(group['data_file']))

    def _group_for_data_file(self, file):
        abs_file = os.path.abspath(file)
        for group in self.scan_files:
            if os.path.abspath(group['data_file']) == abs_file:
                return group
        return self._make_scan_group(abs_file)

    def readscan_heps_id33(self, file):
        import pandas as pd

        group = self._group_for_data_file(file)
        scan_info = self.scan_info(group['data_file'])
        header_attrs = self._read_info(group['info_file'])
        rois = pd.read_csv(group['roi_file'], sep='\t')
        dataframe = pd.read_csv(group['data_file'], sep='\t')

        q_range = self.q_range if self.q_range is not None else header_attrs.get('q_range')
        modules = self.modules if self.modules is not None else header_attrs.get('module_add_list')
        selected_rois = self._select_rois(
            rois,
            header_attrs=header_attrs,
            use_roi_selection=self.use_roi_selection,
            q_range=q_range,
            modules=modules,
            remove_bad_fits=self.remove_bad_fits,
        )
        if selected_rois.empty:
            selected_rois = rois.copy()

        energy_column = self._resolve_energy_column(dataframe)
        intensity_columns, selected_rois = self._resolve_intensity_columns(dataframe, selected_rois, energy_column)
        eloss = dataframe[energy_column].to_numpy(dtype='float64')
        signals = dataframe[intensity_columns].to_numpy(dtype='float64')
        if signals.ndim == 1:
            signals = signals[:, np.newaxis]
        errors = np.sqrt(np.absolute(signals))

        e0_values = selected_rois[self.e0_column].to_numpy(dtype='float64') if self.e0_column in selected_rois else np.array([])
        e0_eV = float(np.nanmean(e0_values)) if e0_values.size else 0.0
        E0 = e0_eV / 1e3
        energy = E0 + eloss / 1e3
        q_values = self._q_values_for_columns(selected_rois, intensity_columns)
        q_average = float(np.nanmean(q_values)) if q_values.size else np.nan
        tth = self._make_tth(q_values, len(intensity_columns), E0)

        scan_attrs = dict(header_attrs)
        scan_attrs.update({
            'beamline': scan_attrs.get('beamline', 'HEPS ID33'),
            'selected_rois': selected_rois['crystal'].tolist() if 'crystal' in selected_rois else [],
            'q_average': q_average,
        })

        onescan = xrs_scans.scan(
            [],
            scan_info[0],
            energy,
            np.ones_like(eloss),
            signals,
            scan_attrs,
            dataframe,
            'nixs',
        )
        onescan.eloss = eloss
        onescan.signals = signals
        onescan.errors = errors
        onescan.cenom = (e0_values / 1e3).tolist() if e0_values.size else [E0] * len(intensity_columns)
        onescan.tth = tth
        onescan.q = q_values

        self.scans[scan_info[1]] = onescan
        self.header_attrs = header_attrs
        self.rois = rois
        self.selected_rois = selected_rois
        self.dataframe = dataframe
        self.info_file = group['info_file']
        self.roi_file = group['roi_file']
        self.data_file = group['data_file']
        self.intensity_columns = intensity_columns
        self.scan_name = scan_info[1]
        self.E0 = E0
        self.cenom = onescan.cenom
        self.q = q_values
        self.q_average = q_average
        self.tth = tth
        self._update_analyzer_metadata(selected_rois, intensity_columns)
        return onescan

    def _resolve_energy_column(self, dataframe):
        if self.energy_column in dataframe.columns:
            return self.energy_column
        columns = [col for col in dataframe.columns if not str(col).startswith('Unnamed:')]
        if not columns:
            raise ValueError('Could not find an energy column in {}'.format(self.data_file))
        return columns[0]

    def _resolve_intensity_columns(self, dataframe, selected_rois, energy_column):
        if self.user_intensity_columns is not None:
            intensity_columns = list(self.user_intensity_columns)
            return intensity_columns, selected_rois

        excluded = set([energy_column])
        excluded.update([col for col in dataframe.columns if str(col).startswith('Unnamed:')])
        candidate_columns = [col for col in dataframe.columns if col not in excluded]
        if 'crystal' not in selected_rois:
            return candidate_columns, selected_rois

        roi_crystals = selected_rois['crystal'].astype(str).tolist()
        crystal_columns = [crystal for crystal in roi_crystals if crystal in dataframe.columns]
        if not crystal_columns:
            return candidate_columns, selected_rois

        selected_rois = selected_rois.set_index('crystal').loc[crystal_columns].reset_index()
        return crystal_columns, selected_rois

    def _q_values_for_columns(self, selected_rois, intensity_columns):
        if self.q_column not in selected_rois:
            return np.array([])
        q_values = selected_rois[self.q_column].to_numpy(dtype='float64')
        if len(q_values) == len(intensity_columns):
            return q_values
        q_average = float(np.nanmean(q_values)) if q_values.size else np.nan
        return np.full(len(intensity_columns), q_average, dtype='float64')

    def _update_analyzer_metadata(self, selected_rois, intensity_columns):
        self.key = {
            'Analyzer%02d' % (idx + 1): idx
            for idx in range(len(intensity_columns))
        }
        self.cenom_dict = {}
        for key, idx in self.key.items():
            if idx < len(selected_rois):
                row = selected_rois.iloc[idx]
                fwhm = float(row['width']) if 'width' in selected_rois else np.nan
                e0 = float(row[self.e0_column]) if self.e0_column in selected_rois else np.nan
            else:
                fwhm = float(np.nanmean(selected_rois['width'])) if 'width' in selected_rois else np.nan
                e0 = float(np.nanmean(selected_rois[self.e0_column])) if self.e0_column in selected_rois else np.nan
            self.cenom_dict[key] = {'average': {'fwhm': fwhm, 'e0': e0}}
            self.data[key] = {}

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

    def _select_rois(self, rois, header_attrs=None, use_roi_selection=True, q_range=None, modules=None, remove_bad_fits=True):
        header_attrs = header_attrs or self.header_attrs
        selected = rois.copy()
        if not use_roi_selection:
            return selected

        if modules and 'crystal' in selected:
            modules = [modules] if isinstance(modules, str) else list(modules)
            prefixes = tuple(module + '-' for module in modules)
            selected = selected[selected['crystal'].astype(str).str.startswith(prefixes)]

        remove_list = set(header_attrs.get('ROIremoveList', []))
        if remove_bad_fits:
            remove_list.update(header_attrs.get('BadFitList', []))
        if remove_list and 'crystal' in selected:
            selected = selected[~selected['crystal'].isin(remove_list)]

        if q_range and self.q_column in selected:
            q_min, q_max = q_range
            selected = selected[(selected[self.q_column] >= q_min) & (selected[self.q_column] <= q_max)]

        return selected

    def _make_tth(self, q_values, number_of_columns, E0=None):
        E0 = self.E0 if E0 is None else E0
        q_values = np.asarray(q_values, dtype='float64')
        if q_values.size == 0:
            q_values = np.full(number_of_columns, np.nan, dtype='float64')
        elif q_values.size != number_of_columns:
            q_values = np.full(number_of_columns, float(np.nanmean(q_values)), dtype='float64')
        if E0 <= 0.0:
            return [np.nan] * number_of_columns

        k = 0.5067731 * E0
        cos_tth = 1.0 - (q_values ** 2) / (2.0 * k ** 2)
        tth = np.degrees(np.arccos(np.clip(cos_tth, -1.0, 1.0)))
        return [float(value) for value in tth]

    def interp_scan_data(self, x, y, common_x):
        y = np.asarray(y)
        if y.ndim == 1:
            return np.interp(common_x, x, y)
        out = np.zeros((len(common_x), y.shape[1]))
        for col in range(y.shape[1]):
            out[:, col] = np.interp(common_x, x, y[:, col])
        return out

    def average_scans(self, chosen_scans):
        scans = [self.scans[self.scan_info(scan)[1]] for scan in chosen_scans]
        lengths = [len(scan.energy) for scan in scans]
        if len(set(lengths)) == 1:
            energy = np.array([scan.energy for scan in scans]).mean(axis=0)
            signals = np.array([scan.signals for scan in scans]).mean(axis=0)
            eloss = np.array([scan.eloss for scan in scans]).mean(axis=0)
            errors = np.array([scan.errors for scan in scans]).mean(axis=0)
            return energy, signals, eloss, errors

        print('Scans have different lengths; interpolating onto a common energy grid before averaging.')
        start = max([scan.energy[0] for scan in scans])
        stop = min([scan.energy[-1] for scan in scans])
        common_energy = np.linspace(start, stop, min(lengths))
        signals = np.array([self.interp_scan_data(scan.energy, scan.signals, common_energy) for scan in scans]).mean(axis=0)
        eloss = np.array([self.interp_scan_data(scan.energy, scan.eloss, common_energy) for scan in scans]).mean(axis=0)
        errors = np.array([self.interp_scan_data(scan.energy, scan.errors, common_energy) for scan in scans]).mean(axis=0)
        return common_energy, signals, eloss, errors

    def load_elastics(self, *args, **kwargs):
        print('HEPS ID33 reduced data: elastic calibration is already stored in the ROI table.')
        return self

    def load_nixs(self, exp_dir=None, scans='all', analyzers='all'):
        if scans == 'all':
            chosen_scans = self.nixs_scans
        elif isinstance(scans, list):
            chosen_scans = [self.nixs_scans[idx] for idx in scans]
        else:
            print("scans must be list of scan indices or all")
            return self

        for file in chosen_scans:
            print("{} {}".format("Reading HEPS ID33 NIXS scan: ", file))
            self.readscan_heps_id33(file)

        if chosen_scans:
            self.energy, self.signals, self.eloss, self.errors = self.average_scans(chosen_scans)
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

    def save_H5(self, H5name='HEPS_ID33_data.H5'):
        if not self.scans and self.nixs_scans:
            self.load_nixs()
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
