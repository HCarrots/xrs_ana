import numpy as np
import traceback , os , re,sys
from datetime import datetime

from xrsana import xrs_utilities,xrs_scans
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

print_citation_message()

class read_heps_id33:
    def __init__(self):
        pass

class read_lerix:
    def __init__(self,exp_dir,elastic_name='elastic',nixs_name='nixs',wide_name='wide',energycolumn=25,monitorcolumn=4,ancolumns=range(6,25)):
        self.scans         = {} # was a dictionary before
        self.path          = os.path.abspath(os.path.split(exp_dir)[0])
        self.monicolumn    = monitorcolumn
        self.encolumn      = energycolumn
        self.ancolumns     = ancolumns
        self.nixs_name, self.nixs_scans       = nixs_name, []
        self.wide_name, self.wide_scans       = wide_name, []
        self.elastic_name, self.elastic_scans = elastic_name, []
        #split scans into NIXS and elastic and begin instance of XRStools scan class for each scan
        if not self.isValidDir(self.path):
            sys.exit(0)
        for file in self.sort_dir():
                scan_info = self.scan_info(file)
                if scan_info[2]=='elastic':
                    self.elastic_scans.append(file)
                if scan_info[2]=='nixs':
                    self.nixs_scans.append(file)
                if scan_info[2]=='wide':
                    self.wide_scans.append(file)
                else:
                    continue
        self.header_attrs  = {} #dictionary of useful information from scan files, inc. e0, comments, scan_time+date
        self.key           = {'Analyzer01':0, 'Analyzer02':1, 'Analyzer03':2,'Analyzer04':3,'Analyzer05':4,
                            'Analyzer06':5,'Analyzer07':6,'Analyzer08':7,'Analyzer09':8,'Analyzer10':9,'Analyzer11':10,'Analyzer12':11
                            ,'Analyzer13':12,'Analyzer14':13,'Analyzer15':14,'Analyzer16':15,'Analyzer17':16,'Analyzer18':17,'Analyzer19':18}
        self.scannumbers   = []
        self.scan_name     = []
        self.sample_name   = []
        self.eloss_avg     = [] #elastic eloss average
        self.signals_avg   = [] #elastic signals average used to plot analyzer resolutions at the end
        self.energy        = []
        self.signals       = []
        self.errors        = []
        self.is_checked    = [] #inserted later to save a list of the chosen analyzers after using .plot_data() save function
        self.tth           = []
        self.resolution    = {}
        self.E0            = []
        self.cenom         = []
        self.cenom_dict    = {}
        self.data          = {}
        #make cenom_dict which holds e0 and fwhm for all analyser
        for key in self.key.keys():
            self.cenom_dict[key] = {}
            self.data[key]       = {}
    ################################################################################
    # Get Ascii Info - parse a file and return the key details
    ################################################################################
    def getfloats(self, txt, allow_times=True):
        """
        function goes through a line and returns the line as a list of strings
        """
        words = [w.strip() for w in txt.replace(',', ' ').split()]
        # mktime = time.mktime
        for i, w in enumerate(words):
            val = None
            try:
                val = float(w)
            except ValueError:
                pass
        words[i] = val
        return(words)

    def colname(self, txt):
        """Function to replace bad characters with '_''s making a line of strings
        easier to handle."""
        return self.fixName(txt.strip().lower()).replace('.', '_')

    def isValidName(self, filename):
        """Function checks that a filename isn't in the list of reserved pythonic
        words. Returns corrected name or False"""
        if filename in RESERVED_WORDS:
            return False
        tnam = filename[:].lower()
        return NAME_MATCH(tnam) is not None

    def fixName(self, filename, allow_dot=True):
        if self.isValidName(filename):
            return filename
        if self.isValidName('_%s' % filename):
            return '_%s' % filename
        chars = []
        valid_chars = VALID_SNAME_CHARS
        if allow_dot:
            valid_chars = VALID_NAME_CHARS
        for s in filename:
            if s not in valid_chars:
                s = '_'
            chars.append(s)
        filename = ''.join(chars)
        # last check (name may begin with a number or .)
        if not self.isValidName(filename):
            filename = '_%s' % filename
        return filename

    def strip_headers(self, headers):
        #reorganise the headers and remove superfluous lines and commentchars
        header = []
        for hline in headers:
            hline = hline.strip().replace('\t', ' ')
            if len(hline) < 1:
                continue
            if hline[0] in COMMENTCHARS:
                hline = hline[1:].lstrip() #assumes reading l2r
            if len(hline) <1:
                continue
            header.append(hline)
        return(header)

    def separate_infile(self, text):
        """Function parses an 20ID ASCII file in reverse and separates it into
        headers, footers and data"""
        _labelline = None
        ncol = None
        dat, footers, headers = [], [], []
        try:
            text.reverse()
        except:
            text[::-1]
        section = 'FOOTER'
        for line in text:
            line = line.strip()
            if len(line) < 1: #remove any blank lines
                continue
            if section == 'FOOTER' and not None in self.getfloats(line):
                section = 'DATA'
            elif section == 'DATA' and None in self.getfloats(line):
                section = 'HEADER'
                _labelline = line
                if _labelline[0] in COMMENTCHARS:
                    _labelline = _labelline[1:].strip()
            if section == 'FOOTER': #reading footers but not using them currently
                footers.append(line)
            elif section == 'HEADER':
                headers.append(line)
            elif section == 'DATA':
                rowdat  = self.getfloats(line)
                if ncol is None:
                    ncol = len(rowdat)
                if ncol == len(rowdat):
                    dat.append(rowdat)
        return(headers, dat, footers)

    def pull_id20attrs(self, header):
        """Function takes headers of 20ID ASCII and parses it for key information
        for header_attrs - N.B. could be shortened by looping through a list of
        important key words rather than doing one by one."""
        bounds, steps, int_times = [], [], []
        header_attrs = {}
        line = -2
        #iterate through the header and pull out useful information and send it to header_attrs Dictionary
        for hhline in map(str.lower,header):
            line = line + 1 #counting to return the user comments which are on the next line
            try:
                if str(header[comment_line].strip()) == 'Scan config:':
                    header_attrs['User Comments'] = ""
                    pass
                else:
                    header_attrs['User Comments'] = str(header[comment_line].strip())
            except:
                pass
            if hhline.startswith('beamline'):
                words = hhline.split('beamline',1)
                header_attrs['beamline'] = str(words[1].strip())
            elif hhline.startswith('e0'):
                if ':' in hhline:
                    words = hhline.split(':',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
                elif '=' in hhline:
                    words = hhline.split('=',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
            elif hhline.startswith('user comment'):
                comment_line = line
            elif "scan time" in hhline:
                #search for scan date and time see: https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                try:
                    words = hhline.split('scan time',1)
                    header_attrs['scan_time'] = datetime.strptime(words[1].strip(), '%H hrs %M min %S sec.').time()
                    header_attrs['scan_date'] = datetime.strptime(words[0].split('panel',1)[1].strip().strip(';'), '%m/%d/%Y  %I:%M:%S %p').date()
                except:
                    continue
            elif "scan bounds" in hhline:
                words = hhline.split('scan bounds',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        bounds.append(float(i))
                    except:
                        pass
                header_attrs['scan_bounds'] = bounds
            elif "scan step(s)" in hhline:
                words = hhline.split('scan step(s)',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        steps.append(float(i))
                    except:
                        pass
                header_attrs['scan_steps'] = steps
            elif "integration times" in hhline:
                words = hhline.split('integration times',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        int_times.append(float(i))
                    except:
                        pass
                header_attrs['int_times'] = int_times
        return(header_attrs)

    def get_col_headers(self, header):
        col_headers = []
        for i in self.colname(header[0]).split('___'): #need three _ to work
            if not i:
                continue
            col_headers.append(i.strip('_'))
        return(col_headers)

    def scan_info(self, f):
        """get the scan number, name, type and file extention from the title of
        the scan assuming typical format e.g. elastic.0001, nixs.0001
        returns:
        [0] -> scan number (e.g. 0001)
        [1] -> scan name (e.g. nixs0001)
        [2] -> scan_type (e.g. nixs)
        [3] -> file name (e.g. nixs.0001)"""
        f = os.path.basename(f) #this allows both directories and files to be passed to get scan_info
        fn,fext = os.path.splitext(f)
        if str.lower(fn)==str.lower(self.nixs_name):
            scan_type = 'nixs'
        elif str.lower(fn)==str.lower(self.elastic_name):
            scan_type = 'elastic'
        elif str.lower(fn)==str.lower(self.wide_name):
            scan_type = 'wide'
        else:
            print(""">> LERIX >> WARNING  \n You have probably called the scan_info
            function without specifying a correct \n <class>.nixs/wide/elastic_name if you
            are calling scan_info manually - you can change this by setting:\n\n
            <class>.nixs_name = '<nixs_name>'""")
            sys.exit()
        scan_number = fext.lstrip('.')
        scan_number = int(scan_number)
        scan_name = scan_type + '%04d' %scan_number
        return(scan_number, scan_name, scan_type, f)

    # def sort_dir(self, path=None):
    #     """Returns a list of directory contents after filtering out scans without
    #     the correct format or size e.g. 'elastic.0001, nixs.0001 '"""
    #     if not path: # allows sort_dir() to be called without a path.
    #         path = self.path
    #     # regular expression search for *.[0-9][0-9][0-9][0-9] in path if is a file (therfore not dir)
    #     res = [f for f in os.listdir(path) if (re.search(r".[0-9]{4}",f)) and (os.path.isfile(os.path.join(path,file)))]
    #     for scan in res:
    #         if not (MIN_FILESIZE>os.stat(os.path.join(path,file)).st_size<MAX_FILESIZE):
    #             print(file, 'is outside of the accepted file limits 1KB -> 100 MB. Skipping.')
    #             res.remove(file)
    #     sorted_dir = sorted([file for file in res if os.path.splitext(file)[0] in [graphite.elastic_name,graphite.nixs_name,graphite.wide_name]]) #case sensitive
    #     return(sorted_dir)

    def sort_dir(self, path=None):
        """Returns a list of directory contents after filtering out scans without
        the correct format or size e.g. 'elastic.0001, nixs.0001 '"""
        if not path: # allows sort_dir() to be called without a path.
            path = self.path
        # regular expression search for *.[0-9][0-9][0-9][0-9] in path if is a file (therfore not dir)
        res = [f for f in os.listdir(path) if (re.search(r".[0-9]{4}",f)) and (os.path.isfile(os.path.join(path,f)))]
        for scan in res:
            if not (MIN_FILESIZE<os.path.getsize(os.path.join(path,scan))<MAX_FILESIZE):
                print(scan, 'is outside of the accepted file limits 1KB -> 100 MB. Skipping.')
                res.remove(scan)
        scan_names = [self.elastic_name.lower(), self.nixs_name.lower(), self.wide_name.lower()]
        sorted_dir = sorted([f for f in res if os.path.splitext(f)[0].lower() in scan_names])
        return(sorted_dir)

    def isValidDir(self,dir):
        """Show that the scan directory is valid, that the directory holds a scan
        with the correct elastic name, nixs name and then let the user know if it
        has not found a wide scan. Returns True if valid directory."""
        if not os.path.isdir(dir):
            print('Check the directory you have supplied')
            return False
        elif not os.path.isfile(os.path.join(dir, self.elastic_name+'.0001')):
            print("The directory you supplied does not have a elastic.0001 file!!! \n If your elastic scan has a different name (*case sensitive), please specify as: 'elastic_name'")
            return False
        elif not os.path.isfile(os.path.join(dir, self.nixs_name+'.0001')):
            print("The directory you supplied does not have a NIXS.0001 file!!! \n If your raman scan has a different name (*case sensitive), please specify as: 'NIXS_name'")
            return False
        elif not os.path.isfile(os.path.join(dir, self.wide_name+'.0001')):
            print("No wide scans found. Continuing...")
            return True
        else:
            return True

    ################################################################################
    # Read Scan
    ################################################################################
    def update_cenom(self, analyzers="all"):
        """Internal Function to get the centre of mass of the elastic peak and
        the E0 for each elastic scan using XRStools"""
        self.cenom = []
        if not self.cenom_dict: #check that the cenom-dict is populated first.
            print('Cenom Dictionary is empty. Please load elastics first!')
            return(0)
        if analyzers == "all":
            print("Running 'Update_Cenom' Script for All analysers")
            analyzers = sorted(self.key.keys())
        elif type(analyzers) is list:
            print("Running 'Update_Cenom' Script for analysers:", analyzers)
            tmp = []
            for i in analyzers:
                i = i - 1
                tmp.append(sorted(self.key.keys())[i])
            analyzers = tmp
        else:
            print("list of analysers for E0 calculation must be a list type.")
            return(False)
        for analyzer in analyzers:
            avg_fwhm, avg_cenom = [],[]
            for scan in self.elastic_scans:
                scan = self.scan_info(scan)[1] # change from elastic.0001 to elastic0001
                avg_cenom.append(self.cenom_dict[analyzer][scan]['e0'])
                avg_fwhm.append(self.cenom_dict[analyzer][scan]['fwhm'])
            # annoying bit of code because nanmean() of a list of just nan's returns a warning, this step avoids ugly (yet harmless) warnings.
            if np.all(np.isnan(avg_cenom)):
                avg_cenom = np.nan
            else:
                avg_cenom = np.nanmean(np.array(avg_cenom)) # was just nanmean
            if np.all(np.isnan(avg_fwhm)):
                avg_fwhm = np.nan
            else:
                avg_fwhm = np.nanmean(np.array(avg_fwhm)) # was just nanmean
            self.cenom_dict[analyzer].update({'average': {'fwhm': avg_fwhm, 'e0': avg_cenom}})
            self.cenom.append(avg_cenom / 1e3) #divide by thousand to go from eV to keV
        self.E0 = np.nanmean(np.array(self.cenom)) # was just nanmean
        print("{} {}".format("E0 was found to be (keV): ", self.E0))
        print("{} {}".format("Average FWHM for the elastics is (eV): ", avg_fwhm))
        for i in range(len(analyzers)):
            if np.isnan(self.cenom[i]):
                print(analyzers[i], 'Elastic peak is less than 100 counts, setting to average e0')
                self.cenom[i] = self.E0
            else:
                continue
        #self.resolution['Resolution'] = round(np.mean(resolution),3)

    def readscan_20ID(self, file):
        """Read an ID20-type ASCII file and return header attributes and data as
        a dictionary. Takes a file path.
        header_attrs -> int_times, scan_steps, scan_bounds, e0, comments, beamline,
                        scan_time, scan_date
        data         -> dictionary of np.array (float64) with callable column names

        New XRStools uses scan class:
        # create an instance of "scan" class for every scan
        onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
        # assign one dictionary entry to each scan
        self.scans[scanname] = onescan
        edfmats: 2D images from pixel counting detector (not required here)
        number: scan number
        energy: energy axis
        monitor: i0 axis
        counters: interpreted as the detector channels
        motors: motor positions (not required)
        data: full data matrix
        scantype: elastic, nixs or long
        """
        import pandas as pd
        scan_info = self.scan_info(file)
        qixs_list = []
        f = open(file, "r") #read starts here
        text = f.read()
        text = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
        headers, dat, footers = self.separate_infile(text)
        try:
            dat = [map(list,zip(*dat))[i][::-1] for i in range(len(dat[1]))] # this function does the inverse ([::-1]) transposition of the dat object, doesn't seem to work in windows
        except:
            dat = [list(map(list,zip(*dat)))[i][::-1] for i in range(len(dat[1]))]
        names = self.get_col_headers(self.strip_headers(headers)) #returns a list of names in the order found in the data file.
        data = pd.DataFrame(np.array(dat).T,columns = np.array(names).T, dtype='float64') #returns a pandas array with the data arranged into labelled columns
        data = np.delete(data.values,np.where(np.diff(data.values[:,self.encolumn]) < 0),axis=0) # Energy values must be strictly increasing!
        tmp_monitor = data[:,self.monicolumn]
        tmp_energy  = data[:,self.encolumn]
        tmp_signals = data[:,self.ancolumns]
        tmp_errors  = np.sqrt(np.absolute(tmp_signals))
        scan_attrs  = self.pull_id20attrs(self.strip_headers(headers)) #get scan_attrs
        if scan_info[2]=='elastic':
            for analyzer in sorted(self.key.keys()): #The analyzer channels in the scan ASCII
                # check counts are high-enough, using XIA filters avoids broadening FWHM
                if 200 <= np.max(tmp_signals[:,self.key[analyzer]]) <= 4000:
                    try:
                        fwhm, cenom = xrs_utilities.fwhm(tmp_energy,tmp_signals[:,self.key[analyzer]])
                    except:
                        print('Elastic scan is wrong shape, skipping')
                        fwhm, cenom = np.nan, np.nan
                else:
                    fwhm, cenom = np.nan, np.nan
                self.cenom_dict[analyzer].update({scan_info[1]: {'fwhm': fwhm, 'e0': cenom}})
        elif scan_info[2]=='nixs' or scan_info[2]=='wide':
            #create empty array with shape energy.v.signals
            eloss = np.zeros(tmp_signals.shape)
            self.tth = list(range(9,180,9)) #assign tth to self
            try:
                e_zero = self.E0 * 1e3 # convert e0 back to eV to perform subtraction
                tmp_eloss = np.subtract(tmp_energy,e_zero)
            except:
                print('>> No elastic <class>.E0 found! Make sure elastic scans have been loaded first')
                e_zero = scan_attrs['e0'] * 1e3 # convert e0 back to eV to perform subtraction
                tmp_eloss = np.subtract(tmp_energy,e_zero)
        # create an instance of "scan" class for every scan
        edfmats, motors = [],scan_attrs #no 2D pixel detector at LERIX
        number          = scan_info[0]
        energy          = tmp_energy
        monitor         = tmp_monitor
        counters        = tmp_signals
        data            = data
        scantype        = scan_info[2]
        onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
        self.scans[scan_info[1]] = onescan
        if scan_info[2]=='nixs' or scan_info[2]=='wide':
            self.scans[scan_info[1]].eloss = tmp_eloss
            self.scans[scan_info[1]].signals = np.divide(tmp_signals.T,monitor).T #transpose seems to be necessary, but don't know why?
            self.scans[scan_info[1]].errors = tmp_errors
        f.close()

    ################################################################################
    # Begin the reading
    ################################################################################
    def _average_scan_attributes(self, chosen_scans, attributes):
        if not chosen_scans:
            raise ValueError('No scans selected for averaging.')

        scan_objs = [self.scans[self.scan_info(scan)[1]] for scan in chosen_scans]
        energies = [np.asarray(scan.energy) for scan in scan_objs]

        same_grid = all(
            energy.shape == energies[0].shape and np.allclose(energy, energies[0])
            for energy in energies
        )
        if same_grid:
            return [
                np.array([np.asarray(getattr(scan, attr)) for scan in scan_objs]).mean(axis=0)
                for attr in attributes
            ]

        common_min = max(energy[0] for energy in energies)
        common_max = min(energy[-1] for energy in energies)
        base_energy = energies[np.argmin([len(energy) for energy in energies])]
        common_energy = base_energy[(base_energy >= common_min) & (base_energy <= common_max)]

        if common_energy.size < 2:
            common_energy = np.linspace(common_min, common_max, min(len(energy) for energy in energies))

        averaged = []
        for attr in attributes:
            interpolated = []
            for scan, energy in zip(scan_objs, energies):
                values = np.asarray(getattr(scan, attr))
                if attr == 'energy':
                    interpolated.append(common_energy)
                elif values.ndim == 1:
                    interpolated.append(np.interp(common_energy, energy, values))
                else:
                    interpolated.append(np.column_stack([
                        np.interp(common_energy, energy, values[:, col])
                        for col in range(values.shape[1])
                    ]))
            averaged.append(np.array(interpolated).mean(axis=0))
        return averaged

    def load_elastics(self,exp_dir=None,scans='all',analyzers='all'):
        """Function to load scan data from a typical APS 20ID Non-Resonant inelastic
        X-ray scattering experiment. With data in the form of elastic.0001, allign.0001
        and NIXS.0001. Function reteurns the averaged energy loss, signals, errors, E0
        and 2theta angles for the scans in the chosen directory."""
        if exp_dir is None:
            exp_dir = self.path
        if scans == 'all':
            chosen_scans = self.elastic_scans
        elif isinstance(scans,list):
            chosen_scans = [self.elastic_scans[i] for i in scans]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading elastic scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)
        self.update_cenom(analyzers)

    def load_nixs(self,exp_dir=None,scans='all',analyzers='all'):
        """Blah Blah"""
        if exp_dir is None:
            exp_dir = self.path
        if scans == 'all':
            chosen_scans = self.nixs_scans
        elif isinstance(scans,list):
            chosen_scans = [self.nixs_scans[i] for i in scans]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading NIXS scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)
        #average the data over the chosen scans
        self.energy, self.signals, self.eloss, self.errors = self._average_scan_attributes(
            chosen_scans,
            ['energy', 'signals', 'eloss', 'errors']
        )

    def load_wides(self,exp_dir=None,scans='all',analyzers='all',join=False):
        """Blah Blah"""
        if exp_dir is None:
            exp_dir = self.path
        if scans == 'all':
            chosen_scans = self.wide_scans
        elif isinstance(scans,list):
            scans[:] = [x - 1 for x in scans] #scan 1 will be the 0th item in the list
            chosen_scans = [self.wide_scans[i] for i in scans]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading Wide scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)
        if join:
            if not np.any(self.eloss):
                try:
                    self.load_nixs(exp_dir=exp_dir,scans=scans,analyzers=analyzers)
                except:
                    print('Must load elastics & nixs first!! >> tried to load NIXS but failed.')
                    return(0)
                self.join_nixs_wide()
            else:
                self.join_nixs_wide()

    def join_nixs_wide(self,scaling='auto'):
        """
        Joins together the loaded nixs_scans and wide_scans so that the compton (below edge) AND
        the nixs scans become one. This is important for compton subtraction.
        inputs:
            scaling == scales wide to nixs data can be 'auto', 'none' or a floating value.
        (i) Check that eloss is scrictly increasing.
        (1) Find indices (idx) of values in wide eloss array closest to the min and max of the nixs eloss array.
        (2) Remove the wide eloss values between these (1) indices. (e.g. cut ready for paste)
        (3) create new array with nixs_eloss values inserted into the wide eloss_values.
        (4) Use the wide_eloss indices (idx) to remove values between the indices for the wide_signals data (analyser wise).
        (5) Insert the nixs_signals into the created hole. (analyser-wise).
        energy,signals,errors,eloss
        TO DO: works but is ugly and can be shortened with a sensible loop.
        """
        # Average wide scan data.
        wide_energy, wide_signals, wide_eloss, wide_errors = self._average_scan_attributes(
            self.wide_scans,
            ['energy', 'signals', 'eloss', 'errors']
        )
        #check to see all values of wide_eloss are strictly increasing
        if not np.all(wide_eloss[1:] >= wide_eloss[:-1], axis=0) and np.all(self.eloss[1:] >= self.eloss[:-1], axis=0):
            print("Eloss in the wide and nixs scan/s must be strictly increasing! stopping.")
            return(0)
        self_eloss = np.asarray(self.eloss)
        wide_eloss = np.asarray(wide_eloss)

        nixs_min = np.min(self_eloss)
        nixs_max = np.max(self_eloss)

        wide_min = wide_eloss[np.argmin(np.abs(wide_eloss - nixs_min))]
        wide_max = wide_eloss[np.argmin(np.abs(wide_eloss - nixs_max))]

        idx = (wide_eloss >= wide_min) & (wide_eloss <= wide_max)
        #idx = (wide_eloss>=min(wide_eloss, key=lambda x:abs(x-min(self.eloss))))*(wide_eloss<=min(wide_eloss, key=lambda x:abs(x-max(self.eloss)))) # idx is boolean with len=wide_eloss, True when in the range of nixs_eloss
        # (A) Generate the common eloss scale.
        wide_eloss = np.delete(wide_eloss,np.where(idx)[0],axis=0)  # deletes wide_eloss values where False.
        try:
            self.eloss = np.insert(wide_eloss,np.where(idx)[0][0],self.eloss,axis=0)  # inserts nixs_eloss at the 0th value of idx
        except:
            self.eloss = np.append(wide_eloss,self.eloss,axis=0)
        # (B) Generate the combined signals with scaling options auto, none, float
        wide_signals_crossover = wide_signals[np.where(idx)[0][0]] # take wide_signal value at crossover point
        wide_signals = np.delete(wide_signals,np.where(idx)[0],axis=0)
        if scaling == 'auto':
            scaling = np.divide(wide_signals_crossover,self.signals[0])
            try:
                self.signals = np.insert(np.divide(wide_signals,scaling),np.where(idx)[0][0],self.signals,axis=0)
            except: # if wide scan doesn't go past the nixs scan, just append nixs onto end of wide.
                self.signals = np.append(np.divide(wide_signals,scaling),self.signals,axis=0)
        elif isinstance(scaling,float):
            scaling = float(scaling)
            try:
                self.signals = np.insert(wide_signals,np.where(idx)[0][0],scaling*self.signals,axis=0)
            except: # if wide scan doesn't go past the nixs scan, just append nixs onto end of wide.
                self.signals = np.append(wide_signals,scaling*self.signals,axis=0)
        elif scaling == 'none':
            try:
                self.signals = np.insert(wide_signals,np.where(idx)[0][0],self.signals,axis=0)
            except: # if wide scan doesn't go past the nixs scan, just append nixs onto end of wide.
                self.signals = np.append(wide_signals,self.signals,axis=0)
        # (C) Generate the common energy scale.
        wide_energy = np.delete(wide_energy,np.where(idx)[0],axis=0)  # deletes wide_eloss values where False.
        try:
            self.energy = np.insert(wide_energy,np.where(idx)[0][0],self.energy,axis=0)  # inserts nixs_eloss at the 0th value of idx
        except:
            self.energy = np.append(wide_energy,self.energy,axis=0)
        # (D) Generate the common errors scale.
        wide_errors = np.delete(wide_errors,np.where(idx)[0],axis=0)  # deletes wide_eloss values where False.
        try:
            self.errors = np.insert(wide_errors,np.where(idx)[0][0],self.errors,axis=0)  # inserts nixs_eloss at the 0th value of idx
        except:
            self.errors = np.append(wide_errors,self.errors,axis=0)
        print("Successfully joined the wide and nixs scans.")

    def plot_data(self,analyzer=False):
        """<classObj>.plot_data() Function that can be called to plot the eloss
        data for each channel and build an average by clicking a button.
        Does not require matplotlib >2.1"""
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Button, Cursor
        channels = []
        for analyzer in self.resolution:
            if analyzer.startswith('Analyzer'):
                if self.resolution[analyzer] < 1.0:
                    channels.append(int(analyzer.lstrip('Analyzer'))-1)
        data = np.average(self.signals[:,channels],axis=1)
        fig, ax = plt.subplots()
        ax.plot(self.eloss, data, lw=2)
        ax.set_xlabel('Energy Loss (eV)')
        ax.set_ylabel('S(q,w) [1/eV]')
        ax.set_title('Plotting Raman Analysers')
        plt.subplots_adjust(left=0.3)
        checkbuttonaxis = plt.axes([0.02, 0.15, 0.2, 0.8])
        anlabels, anvals = list(self.key), (False,)*len(list(self.key))
        anstates = dict(zip(anlabels,anvals))
        analyzers = CheckButtons(checkbuttonaxis, anlabels, anvals)
        buttonaxis = plt.axes([0.01, 0.01, 0.3, 0.09])
        bStatus  = Button(buttonaxis,'Save Averaged Analyzers')

        def onclick(label):
            """Tell the user what they have clicked - also good for de-bugging
            """
            anstates[label] = not anstates[label]
            print('un'*(not anstates[label]) + 'checked %s' %label)
            func()

        def savebutton(val):
            import pandas as pd
            if not self.is_checked:
                print('please select your chosen analysers first!')
            else:
                print('selected analysers (python counting):  ', self.is_checked)
                save_signals = np.average(self.signals[:,self.is_checked],axis=1)
                save_errors = np.average(self.errors[:,self.is_checked],axis=1)
                df = pd.DataFrame(list(zip(self.eloss,save_signals,save_errors)), columns=['eloss','signals','errors'])
                print(df)
                try:
                    filename = getattr(self, 'average_save_path', os.path.join(self.path, 'Result.csv'))
                    df.to_csv(filename,sep=',',na_rep='nan')
                    print('Saved as: ',filename)
                except:
                    print("{} {}".format(">> Warning >>", "file save was unsuccessful"))

        def func():
            ax.clear()
            self.is_checked = []
            for ii in anlabels:
                if anstates[ii]:
                    self.is_checked.append(self.key[ii])
            ax.plot(self.eloss, np.average(self.signals[:,self.is_checked],axis=1),lw=2)
            cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
            ax.autoscale(True)
            ax.set_xlabel('Energy Loss (eV)')
            ax.set_title('Plotting Raman Analysers')
            plt.draw()

        bStatus.on_clicked(savebutton)
        analyzers.on_clicked(onclick)
        cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
        plt.show()

    def save_H5(self,H5name='20ID_APS_data.H5'):
        #if the user asks, call function to write all info to H5 file
        H5path = os.path.join(self.path,H5name)
        if not os.path.isdir(os.path.dirname(H5path)):
            print('H5 path directory does not exist!')
        if os.path.isfile(H5path):
            H5file = h5py.File(H5path, "a")
        else:
            H5file = h5py.File(H5path, "w")
        if H5name in H5file:
            del H5file[H5name]
        g = H5file.create_group(H5name) #H5 subgroup with the name of the sample
        H5_ela = g.create_group('elastic') #H5 subgroup for elastics
        H5_xrs = g.create_group('XRS')     #H5 subgroup for NIXS
        all_scans = self.elastic_scans + self.nixs_scans + self.wide_scans
        for file in all_scans:
            scan_info = self.scan_info(file)
            if scan_info[2] == 'elastic':
                h5group = H5_ela.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("cenoms",data=self.scans[scan_info[1]].cenom)
            elif scan_info[2]=='nixs':
                h5group = H5_xrs.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("eloss",data=self.scans[scan_info[1]].eloss)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
            elif scan_info[2]=='wide':
                h5group = H5_xrs.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("eloss",data=self.scans[scan_info[1]].eloss)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
        g.create_dataset("energy",data=self.energy)
        g.create_dataset("signals",data=self.signals)
        g.create_dataset("eloss",data=self.eloss)
        g.create_dataset("errors",data=self.errors)
        g.create_dataset("tth",data=self.tth)
        #g.create_dataset("Mean Resolutions", data=np.array(self.resolution.items()))
        #Never forget to close an open H5 file!!!
        H5file.close()
