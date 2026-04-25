#!/usr/bin/python
# Filename: xrs_scans.py
#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"


import numpy as np
from . import xrs_utilities, math_functions
import h5py

from itertools import groupby
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import Rbf
from matplotlib import pylab as plt
from scipy import ndimage


class Scan:
    """ **Scan**

    Class for manipulating scan data from the Hydra and Fourc spectrometers. 

    Attributes:
        edf_mats      (np.array): Array containing all 2D detector images that belong to the scan.
        number             (int): Scan number.
        scan_type       (string): Keyword, used later to group scans (add similar scans, etc.).
        energy        (np.array): Array containing the energy axis (1D).
        monitor       (np.array): Array containing the monitor signal (1D).
        counters    (dictionary): Counters with assiciated scan data.
        motors            (list): Motor positions for the scan.
        eloss         (np.array): Array of the energy loss scale.
        signals       (np.array): Array of signals extracted from the ROIs.
        errors        (np.array): Array of Poisson errors.
        cenom             (list): Center of mass for each ROI (used if scan is an elastic line scan).
        signals_pw        (list): Pixelwise (PW) data, one array of PW data per ROI.
        errors_pw         (list): Pixelwise (PW) Poisson errors, one array of PW errors per ROI.
        cenom_pw          (list): Center of mass for each pixel.
        signals_pw_interp (list): Interpolated signals for each pixel.
        Ein              (float): Incident energy, used if energy2 is scanned.
        raw_signals       (dict): Dictionary of raw data (summed, line-by-line, pixel-by-pixel).
        raw_errors        (dict): Dictionary of raw erros (summed, line-by-line, pixel-by-pixel).

    """
    normalizationDict={} 
    def __init__( self ):
        self.edfmats    = np.array([])
        self.scan_number= None
        self.scan_type  = None
        self.energy     = np.array([])
        self.monitor    = np.array([])
        self.counters   = {}
        self.motors     = []
        self.eloss      = np.array([])
        self.signals    = np.array([])
        self.errors     = np.array([])
        self.cenom      = []
        self.signals_pw = []
        self.errors_pw  = []
        self.cenom_pw   = []
        self.signals_pw_interp = []
        self.Ein        = None
        self.scan_motor = None

        self.raw_signals = {}
        self.raw_errors  = {}

        self.__signals_normalized__ = False

    def assign( self, detector_arrays, scan_number, energy_scale, monitor_signal, counters, \
                motor_positions, scan_data, scan_type='generic' ):
        """ **assign**

        Method to group together existing data from a scan (for backward compatibility).

        Args:
            detector_arrays (np.array): Array of all 2D images that belong to the scan.
            scan_number         (int): Scan number.
            energy_scale   (np.array): Array of the energy axis.
            monitor_signal (np.array): Array of the monitor signal.
            counters     (dictionary): Counters with associated scan data.
            motor_positions    (list): Motor positions for the scan.
            scan_data     (np.array): Matrix with scan data.
            scantype            (str): Keyword, used later to group scans (add similar scans, etc.).

        """
        self.edfmats    = np.array(detector_arrays)
        self.scan_number= scan_number
        self.energy     = np.array(energy_scale)
        self.monitor    = np.array(monitor_signal)
        self.counters   = counters
        self.motors     = motor_positions
        self.scan_data  = scan_data
        self.scan_type  = scan_type

    def save_hdf5( self, fname ):
        """ **save_hdf5**
        Save a scan in an HDF5 file.
        Note:
            HDF5 files are strange for overwriting files.

        Args:
            fname (str): Path and filename for the HDF5 file.
        """
        if isinstance(fname, h5py.Group):
            f=fname
        else:
            f = h5py.File(fname, "w")

        h5_md = f.require_group("motorDict")
        for mn, mv in self.motors.items():
            h5_md[mn] = mv
            
        h5_md = f.require_group("counters")
        for mn, mv in self.counters.items():
            h5_md[mn] = mv
            
        for attr in ['edfmats', 'scan_number', 'energy', 'monitor', 'scan_type','signals','errors']:
            f[attr] = eval( 'self.' + attr )
            
        for key in self.used_masks.keys():
            hgroup = f.require_group(key)
            pos, mask  = self.used_masks[key]
            hgroup["mask"] = mask
            hgroup["mask_pos"] = pos
            if hasattr( self, "insets") :
                hgroup["insets"]  = self.insets  [key]
        if not isinstance(fname, h5py.Group):
            f.close()

    def load_hdf5( self, fname ):
        """ **load_hdf5**
        Load a scan from an HDF5 file.
        Args:
            fname (str): Filename of the HDF5 file.
        """
        if isinstance(fname, h5py.Group):
            f=fname
        else:
            f = h5py.File(fname, "r")

        for attr in ['edfmats', 'scan_number', 'energy', 'monitor', 'counters', 'motors', 'scan_type']:
           eval( 'self.' + attr +' = f[attr]')

        self.used_masks = {}
        for key in f:
            print(" has key ", key)
            if str(key)[:3] == "ROI" :
                mygroup = f[key]
                self.used_masks[key] = [mygroup["mask_pos"][:].tolist(),  np.array(mygroup["mask"][:])] 


           
        f.close()

        

    def get_raw_signals( self, roi_obj, method='sum', scaling=None, rot_angles=None, storeInsets=False ):
        """ **get_raw_signals**

        Applies given ROIs to detector images. 

        Applies the provided ROIs to detector images in the specified mannar:
        summing, line-by-line, or pixel-by-pixel. Depending on the choice, the
        resulting data array is 2D (sum), 3D (line-by-line), or 4D (pixel-by-pixel).
        The scanned direction is always the first dimension of the resulting data 
        matrix.

        Args:
            roi_obj (instance): Instance of the 'XRStools.xrs_rois.roi_object' class defining the ROIs.
            method    (string): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', 'pixel', or 'column'. Default is 'sum'.
            scaling (np.array): Array of float-type scaling factors (factor for each ROI).
        Returns:
            None if there are no detector images to apply the ROIs to.

        """


        self.used_masks={}

        if storeInsets:
            self.insets={}
            for key, (pos, M) in roi_obj.red_rois.items():
                S     = M.shape
                self.insets[key] =  np.zeros(   [  len(self.edfmats), S[0] , S[1] ],  self.edfmats.dtype   ) 


        for ii in range(len(self.edfmats)):
            for key, (pos, M) in roi_obj.red_rois.items():
                S     = M.shape
                inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
                self.used_masks[key] = (   pos, M )
                if storeInsets:
                    self.insets[key][ii] = self.edfmats[ii, inset[0], inset[1]] * (M/M.max())

                
        # sum
        if method == 'sum':
            print('selected method is \'sum\': summing up pixels from each ROI.')
            signals = {} # dict (one entry per ROI, with vector (one entry per energy point))
            errors  = {} # sqrt of the sum of counts
            for key, (pos, M) in roi_obj.red_rois.items():
                signals[key] = np.zeros((len(self.energy)))
                errors[key]  = np.zeros((len(self.energy)))

            for ii in range(len(self.edfmats)):
                ind = 0
                for key, (pos, M) in roi_obj.red_rois.items():
                    S     = M.shape
                    inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
                    
                    signals[key][ii] = np.sum( self.edfmats[ii, inset[0], inset[1]] * (M/M.max()))
                    errors[key][ii]  = np.sqrt(signals[key][ii])
                    signals[key][ii] /= self.monitor[ii]
                    errors[key][ii]  /= self.monitor[ii]
                    ind += 1

        # row
        elif method == 'row':
            print('selected method is \'row\': summing over non-dispersive direction for each ROI.')
            signals = {} # dict (one entry per ROI, with 2D matrix (energy vs row))
            errors  = {} # sqrt of the sum of counts
            rot_angles_dict = {} # put possible rotation angles into dict
            counter = 0
            for key, (pos, M) in sorted(roi_obj.red_rois.items()):
                signals[key] = np.zeros((len(self.energy), M.shape[0]))
                errors[key]  = np.zeros((len(self.energy), M.shape[0]))
                if rot_angles is not None:
                    rot_angles_dict[key] = rot_angles[counter]
                    counter += 1

            for ii in range(len(self.edfmats)):
                ind = 0
                for key, (pos, M) in roi_obj.red_rois.items():
                    S     = M.shape
                    inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
                    
                    
                    # rotate raw_signals and raw_errors if method is 'line' and angles are provided
                    if rot_angles:
                        if len(roi_obj.red_rois) is not len(rot_angles):
                            print('Only %d rotation angles provided for %d ROIs. Will end here.'%(len(rot_angles), len(self.raw_signals)))
                            return

                        # rotate images before summation
                        orig_slice    = self.edfmats[ii, inset[0], inset[1]] * (M/M.max())
                        slice_for_sum = ndimage.interpolation.rotate( orig_slice, rot_angles_dict[key],\
                                        reshape=False, order=0, mode='constant' )
                    else:
                        slice_for_sum = self.edfmats[ii, inset[0], inset[1]] * (M/M.max())
                    signals[key][ii,:] = np.sum( slice_for_sum , axis=1)
                    errors[key][ii,:]  = np.sqrt(signals[key][ii,:])
                    signals[key][ii,:] /= self.monitor[ii]
                    errors[key][ii,:]  /= self.monitor[ii]
                    ind += 1

        # column
        elif method == 'column':
            print('selected method is \'column\': summing over dispersive direction for each ROI.')
            signals = {} # dict (one entry per ROI, with 2D matrix (energy vs row))
            errors  = {} # sqrt of the sum of counts
            for key, (pos, M) in roi_obj.red_rois.items():
                signals[key] = np.zeros((len(self.energy), M.shape[1]))
                errors[key]  = np.zeros((len(self.energy), M.shape[1]))

            for ii in range(len(self.edfmats)):
                ind = 0
                for key, (pos, M) in roi_obj.red_rois.items():
                    S     = M.shape
                    inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))

                    
                    signals[key][ii,:] = np.sum( self.edfmats[ii, inset[0], inset[1]] * (M/M.max()), axis=0)
                    errors[key][ii,:]  = np.sqrt( signals[key][ii,:] )
                    signals[key][ii,:] /= self.monitor[ii]
                    errors[key][ii,:]  /= self.monitor[ii]
                    ind += 1

        # pixel
        elif method == 'pixel' or method == 'pixel2':
            print('selected method is \'pixel\': returning ROI pixel-by-pixel.')
            signals = {} # dict (one entry per ROI, with 3D matrix (energy vs pixel_0 vs pixel_1))
            errors  = {} # sqrt of the sum of counts
            for key, (pos, M) in roi_obj.red_rois.items():
                signals[key] = np.zeros((len(self.energy), M.shape[0], M.shape[1]))
                errors[key]  = np.zeros((len(self.energy), M.shape[0], M.shape[1]))

            for ii in range(len(self.edfmats)):
                ind = 0
                for key, (pos, M) in roi_obj.red_rois.items():
                    S     = M.shape
                    inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))

                    
                    signals[key][ii,:,:] = self.edfmats[ii, inset[0], inset[1]] * (M/M.max())
                    errors[key][ii,:,:]  = np.sqrt( signals[key][ii,:,:]  )
                    signals[key][ii,:,:] /= self.monitor[ii]
                    errors[key][ii,:,:]  /= self.monitor[ii]
                    ind += 1

        # unknown method
        else:
            print( 'Unknown integration method. Use either \'sum\', \'row\', or \'pixel\'.' )
            return

        # set normalization 
        self.__signals_normalized__ = True

        # assign
        for key,ii in zip(sorted(roi_obj.red_rois), range(len(roi_obj.red_rois))):
            if np.any(scaling):
                self.raw_signals[key] = signals[key] * scaling[ii]
                self.raw_errors[key]  = errors[key]  * scaling[ii]
            else:
                self.raw_signals[key] = signals[key]
                self.raw_errors[key]  = errors[key]

        # clear image stack after pixel-by-pixel integration
        print( 'Clearing detector image stack of scan No. %s.' % self.scan_number )
        self.edfmats = np.array([])

    def get_signals( self, method='sum', cenom_dict=None, comp_factor=None, scaling=None, PIXEL_SIZE=0.055 ):
        """ **get_signals**

        Turns pixel-, column- or sum-wise raw-data into data.

        Takes the raw-data after application of the ROIs and applies the chosen compensation
        scheme.

        Args:
            method        (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.
            cenom_dict   (dict): Dictionary with one entry per ROI holding information about
                the center of mass of the according elastic line.
            comp_factor (float): Factor used in the RIXS-style line-by-line compensation.
            scaling  (np.array): Array of float-type scaling factors (factor for each ROI).

        """

        if method == 'sum':
            self.signals = np.zeros(( len(self.energy), len(self.raw_signals) ))
            self.errors  = np.zeros(( len(self.energy), len(self.raw_signals) ))
            for key,ii in zip(sorted(self.raw_signals), range(len(self.raw_signals))):
                if not self.__signals_normalized__:
                    self.signals[:,ii] = self.raw_signals[key]/self.monitor
                    self.errors[:,ii]  = self.raw_errors[key]/self.monitor
                else:
                    self.signals[:,ii] = self.raw_signals[key]
                    self.errors[:,ii]  = self.raw_errors[key]

        elif method == 'pixel':
            # assign sign of compensation
            if self.scan_motor == 'energy':
                comp_direction = 1.0
            elif self.scan_motor == 'anal energy':
                comp_direction = -1.0
            else:
                print('Unknown energy motor, will break here.')
                return
            self.signals = np.zeros(( len(self.energy), len(self.raw_signals) ))
            self.errors  = np.zeros(( len(self.energy), len(self.raw_signals) ))
            for key,ii in zip(sorted(self.raw_signals), range(len(self.raw_signals))):
                S = cenom_dict[key].shape
                master_cenom =  cenom_dict[key][int(S[0]/2.),int(S[1]/2.)]
                for dim1 in range(self.raw_signals[key].shape[1]):
                    for dim2 in range(self.raw_signals[key].shape[2]):
                        x    = self.energy + ( comp_direction*cenom_dict[key][dim1,dim2] - comp_direction*master_cenom )
                        y    = self.raw_signals[key][:,dim1,dim2]
                        rbfi = Rbf( x, y, function='linear' )
                        self.signals[:,ii] += rbfi(self.energy)
                        y    = self.raw_errors[key][:,dim1,dim2]
                        rbfi = Rbf( x, y, function='linear' )
                        self.errors[:,ii] += rbfi(self.energy)**2
                if not self.__signals_normalized__:
                    self.errors[:,ii] = np.sqrt( self.errors[:,ii] )/self.monitor
                    self.signals[:,ii] /= self.monitor
                else:
                    self.errors[:,ii] = np.sqrt( self.errors[:,ii] )

        elif method == 'pixel2':
            # assign sign of compensation
            if self.scan_motor == 'energy':
                comp_direction = 1.0

            elif self.scan_motor == 'anal energy':
                comp_direction = -1.0
            else:
                print('Unknown energy motor, will break here.')
                return
            self.signals = np.zeros(( len(self.energy), len(self.raw_signals) ))
            self.errors  = np.zeros(( len(self.energy), len(self.raw_signals) ))
            energy = self.energy #* 1e3 # energy in eV
            ## meanmon = np.mean(self.monitor)
            for key,ii in zip(sorted(self.raw_signals), range(len(self.raw_signals))):
                S = cenom_dict[key].shape
                #the_hist = np.histogram(cenom_dict[key][cenom_dict[key]>0.0 ], bins=10)
                #master_cenom =  np.average(the_hist[1][1:], weights= the_hist[0])#cenom_dict[key][int(S[0]/2.),int(S[1]/2.)]
                master_cenom = np.mean(cenom_dict[key][cenom_dict[key]>0.0 ])
                y  = self.raw_signals[key]
                yn = y# (y.T/self.monitor).T * meanmon
                yc = np.zeros(( len(energy),S[0]*S[1] ))
                for dim1 in range(S[0]):
                    for dim2 in range(S[1]):
                        sort  = np.argsort(energy)
                        if cenom_dict[key][dim1,dim2] > 0.0:
                            yc[:,dim1*dim2] = np.interp( energy[sort] -(comp_direction*(cenom_dict[key][dim1,dim2] - \
                                                    master_cenom)), energy[sort], yn[sort,dim1,dim2], \
                                                    left=float('nan'), right=float('nan') )
                        else:
                            yc[:,dim1*dim2] = yn[sort,dim1,dim2]
                self.signals[:,ii] = xrs_utilities.nonzeroavg(yc)
                self.errors[:,ii]  = np.sqrt(self.signals[:,ii])
                # and normalize to I0
                if not self.__signals_normalized__:
                    self.signals[:,ii] /= self.monitor
                    self.errors[:,ii] /= self.monitor

        elif method == 'row':

            # assign sign of compensation !!! test this !!!
            if self.scan_motor == 'energy':
                comp_direction = 1.0
            elif self.scan_motor == 'anal energy':
                comp_direction = -1.0
            else:
                print('Unknown energy motor, will break here.')
                return
            self.signals = np.zeros(( len(self.energy), len(self.raw_signals) ))
            self.errors  = np.zeros(( len(self.energy), len(self.raw_signals) ))
            energy = self.energy * 1e3 # energy in eV
            ### meanmon = np.mean(self.monitor)
            for key,ii in zip(sorted(self.raw_signals), range(len(self.raw_signals))):
                y  = self.raw_signals[key]
                yn = (y.T/self.monitor).T #### * meanmon CHECK THIS REMOVAL
                meanii = len(range(yn.shape[1]))//2
                yc = np.zeros_like(y)
                for jj in range(yn.shape[1]):
                    sort  = np.argsort(energy)
                    yc[:,jj] = np.interp( energy[sort] + (jj-meanii)*comp_direction*comp_factor*PIXEL_SIZE, energy[sort], yn[sort,jj], 
                                        left=float('nan'), right=float('nan') )
                self.signals[:,ii] = xrs_utilities.nonzeroavg(yc)
                self.errors[:,ii]  = np.sqrt(self.signals[:,ii])
                # and normalize to I0
                if not self.__signals_normalized__:
                    self.signals[:,ii] /= self.monitor
                    self.errors[:,ii] /= self.monitor

        else:
            print( 'Unknown integration method. Use either \'sum\', \'row\', or \'pixel\'.' )
            return

        # set normalization
        self.__signals_normalized__ = True

    def apply_rois( self, roi_obj, scaling=None ):
        """ **apply_rois**

        Sums up intensities in each ROI.

        Note:
            Old function, keeping for backward compatibility.

        Args:
            roi_obj (instance): Instance of the 'XRStools.xrs_rois.roi_object' class defining the ROIs.
            scaling (np.array): Array of float-type scaling factors (factor for each ROI).

        Returns:
            None if there are no detector images to apply the ROIs to.

        """

        # make sure there is data loaded
        if self.edfmats.size == 0:
            print( 'Please assign detector images first.' )
            return

        # apply ROIs
        signals = np.zeros((len(self.energy), len(roi_obj.red_rois)))
        for ii in range(len(self.edfmats)):
            ind = 0
            for key, (pos, M) in roi_obj.red_rois.items():
                S     = M.shape
                inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
                signals[ii,ind] = np.sum( self.edfmats[ii, inset[0], inset[1]] * (M/M.max()))
                ind += 1

        # assign
        self.signals = signals
        self.errors  = np.sqrt(signals)

        # apply scaling (if applicable)
        if np.any(scaling):
            # make sure, there is one scaling factor for each ROI
            assert len(scaling) == signals.shape[1]
            for ii in range(signals.shape[1]):
                self.signals[:,ii] *= scaling[ii]
                self.errors[:,ii]  *= scaling[ii]

    def apply_rois_pw( self,roi_obj, scaling=None ):
        """ **apply_rois_pw**

        Pixel-wise reading of the ROIs' pixels into a list of arrays.

        I.e. each n-pixel ROI will have n Spectra, saved in a 2D array.

        Args:
        roi_obj (instance): Instance of the 'XRStools.xrs_rois.roi_object' class defining the ROIs.
        scaling (list) or (np.array): Array or list of float-type scaling factors (one factor for each ROI).

        """
        data   = [] # list of 2D arrays (energy vs. intensity for each pixel inside a single ROI) 
        errors = []
        for n in range(len(roi_obj.indices)): # each ROI
            roidata = np.zeros((len(self.edfmats),len(roi_obj.indices[n]))) # 2D np array energy vs pixels in current roi
            for m in range(len(self.edfmats)): # each energy point along the scan
                for l in range(len(roi_obj.indices[n])): # each pixel on the detector
                    roidata[m,l] = self.edfmats[m,roi_obj.indices[n][l][0],roi_obj.indices[n][l][1]]
            data.append(roidata) # list which contains one array (energy point, pixel) per ROI
            errors.append(np.sqrt(roidata))

        self.signals_pw = data
        self.errors_pw  = errors

        # apply scaling (if applicable)
        if np.any(scaling):
            # make sure, there is one scaling factor for each ROI
            assert len(scaling) == len(roi_obj.indices)
            for ii in range(len(roi_obj.indices)):
                self.signals_pw[ii] *= scaling[ii]
                self.errors_pw[ii]  *= scaling[ii]

    def append_scan( self, scan, method='sum', where='right' ):
        """ **append_scan**

        Appends scan to the current scan, either at higher energies 
        (where='right') or at lower energies (where='left').

        Args:
            scan       (obj): Object of the Scan class.
            method     (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.
            where      (str): Keyword specifying if the scan should be appended
                at lower (where='left') or highger energies (where='righ'). Default is 
                'right'.        
 
        """

        if where == 'right':

            self.energy  = np.append( self.energy , scan.energy )
            self.monitor = np.append( self.monitor, scan.monitor )
            for key in self.raw_signals:
                self.raw_signals[key] = np.append( self.raw_signals[key], scan.raw_signals[key], axis=0 )
                self.raw_errors[key]  = np.append( self.raw_errors[key], scan.raw_errors[key], axis=0 )
        if where == 'left':

            self.energy  = np.append( scan.energy,  self.energy )
            self.monitor = np.append( scan.monitor, self.monitor )
            for key in self.raw_signals:
                self.raw_signals[key] = np.append( scan.raw_signals[key], self.raw_signals[key], axis=0 )
                self.raw_errors[key]  = np.append( scan.raw_errors[key], self.raw_errors[key], axis=0 )

    def insert_scan( self, scan, method='sum', where=None ):
        """ **insert_scan**

        Inserts another scan into the current instance.

        Args:
            scan       (obj): Object of the Scan class.
            method     (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.
            where     (list): Optional tuple of energy values (high and low (in keV))
                for where to insert the scan. By default (None), the lowest and highest
                energy values of the given scan will be used. 

        """

        # find where given scan should be inserted
        if not where:
            low_inds  = np.where( self.energy < scan.get_E_start() )[0]
            high_inds = np.where( self.energy > scan.get_E_end() )[0]
        else:
            low_inds  = np.where( self.energy < where[0] )[0]
            high_inds = np.where( self.energy > where[1] )[0]       

        # 'sum'
        if method == 'sum':

            self.energy  = np.append( self.energy[low_inds], np.append( scan.energy, self.energy[high_inds] ) )
            self.monitor = np.append( self.monitor[low_inds], np.append( scan.monitor, self.monitor[high_inds] ) )
            for key in self.raw_signals:
                self.raw_signals[key] = np.append( self.raw_signals[key][low_inds], 
                                        np.append( scan.raw_signals[key], self.raw_signals[key][high_inds], axis=0 ), axis=0 )
                self.raw_errors[key] = np.append( self.raw_errors[key][low_inds], 
                                        np.append( scan.raw_errors[key], self.raw_errors[key][high_inds], axis=0 ), axis=0 )

        # 'pixel'
        elif method in [ 'pixel' , 'pixel2']:

            self.energy  = np.append( self.energy[low_inds], np.append( scan.energy, self.energy[high_inds] ) )
            self.monitor  = np.append( self.monitor[low_inds], np.append( scan.monitor, self.monitor[high_inds] ) )
            for key in self.raw_signals:
                self.raw_signals[key] = np.append( self.raw_signals[key][low_inds,:,:], 
                                        np.append( scan.raw_signals[key], self.raw_signals[key][high_inds,:,:], axis=0 ), axis=0 )
                self.raw_errors[key] = np.append( self.raw_errors[key][low_inds,:,:], 
                                        np.append( scan.raw_errors[key], self.raw_errors[key][high_inds,:,:], axis=0 ), axis=0 )

        # 'row'
        elif method == 'row':
            self.energy  = np.append( self.energy[low_inds], np.append( scan.energy, self.energy[high_inds] ) )
            self.monitor  = np.append( self.monitor[low_inds], np.append( scan.monitor, self.monitor[high_inds] ) )
            for key in self.raw_signals:
                self.raw_signals[key] = np.append( self.raw_signals[key][low_inds,:], 
                                        np.append( scan.raw_signals[key], self.raw_signals[high_inds,:], axis=0 ), axis=0 )
                self.raw_errors[key] = np.append( self.raw_errors[key][low_inds,:], 
                                        np.append( scan.raw_errors[key], self.raw_errors[high_inds,:], axis=0 ), axis=0 )

        else:
            print( 'Unknown method. Use either \'sum\', \'row\', or \'pixel\'.' )
            return

    def add_scan(self, scan, method='sum', interp=False):
        """ **add_scan**

        Adds signals from a different scan.

        Args:
            scan       (obj): Object of the Scan class.
            method     (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.
            interp (boolean): Boolean specifying if norm-conserving linear wavelet
                interpolation should be used, False by default. 

        """
        # @@@@@@@@@@@@@@@@@@   DA VERIFICARE TUTTO L'ACCOUNTING
        # @@@@@@@@@@@@@@@@@@   changed by christoph: 13/07/2018
        # @@@@@@@@@@@@@@@@@@        please double check

        if not interp:
            
            # sum monitors
            av_monitor = np.sum((self.monitor, scan.monitor), axis=0)

            for key in self.raw_signals:
                if method in [ 'sum', 'row', 'column']:

                    # recover raw counts
                    signals1 = self.raw_signals[key]*self.monitor
                    signals2 = scan.raw_signals[key]*scan.monitor

                    # sum signals
                    av_signals = np.sum((signals1, signals2), axis=0)

                    # errors of summed signals
                    av_errors  = np.sqrt(av_signals)

                    # assign and renormalize
                    self.raw_signals[key] = av_signals/av_monitor
                    self.raw_errors[key]  = av_errors/av_monitor

                if method in [ 'pixel', 'pixel2']:

                    # recover raw counts
                    signals1 = self.raw_signals[key]*self.monitor[:,None,None]
                    signals2 = scan.raw_signals[key]*scan.monitor[:,None,None]

                    # sum signals
                    av_signals = np.sum((signals1, signals2), axis=0)
              
                    # errors of summed signals
                    av_errors  = np.sqrt(av_signals)

                    # assign and renormalize
                    self.raw_signals[key] = av_signals/av_monitor[:,None,None]
                    self.raw_errors[key]  = av_errors/av_monitor[:,None,None]

        if interp:
            # find longest scan dimension
            dim0 = np.amax([len(self.energy), len(scan.energy)])

            # sum the monitor signal (always 1D)
            mm = np.zeros((dim0, 2))*np.nan
            mm[0:len(self.monitor),0] = self.monitor
            mm[0:len(scan.monitor),1] = scan.monitor
            av_monitor = np.nansum( mm , axis=1)
            
            # add also unfinished scans
            for key in self.raw_signals:

                # construct a matrix for nanmean depending on the shape of the raw_signals
                if method in ['sum']:
                    # recover raw counts
                    signals1 = self.raw_signals[key]*self.monitor
                    signals2 = scan.raw_signals[key]*scan.monitor
                    # signals
                    yy = np.zeros((dim0, 2))*np.nan
                    yy[0:len(signals1),0] = signals1
                    yy[0:len(signals2),1] = signals2
                    av_signals = np.nansum( yy , axis=1)
                    # errors of summed signals
                    av_errors  = np.sqrt(av_signals)
                    #
                    self.raw_signals[key] = av_signals/av_monitor
                    self.raw_errors[key]  = av_errors/av_monitor

                if method in ['pixel', 'pixel2']:
                    # recover raw counts
                    signals1 = self.raw_signals[key]*self.monitor[:,None,None]
                    signals2 = scan.raw_signals[key]*scan.monitor[:,None,None]
                    # signals
                    y1 = np.zeros((dim0, signals1.shape[1], signals1.shape[2]))*np.nan
                    y2 = np.zeros((dim0, signals1.shape[1], signals1.shape[2]))*np.nan
                    y1[0:len(signals1),:,:] = signals1
                    y2[0:len(signals2),:,:] = signals2
                    av_signals = np.nansum( (y1,y2) , axis=0 )
                    # errors of summed signals
                    av_errors  = np.sqrt( av_signals )
                    #
                    self.raw_signals[key] = av_signals/av_monitor[:,None,None]
                    self.raw_errors[key]  = av_errors/av_monitor[:,None,None]

                if method in ['row']:
                    # recover raw counts
                    signals1 = self.raw_signals[key]*self.monitor[:,None]
                    signals2 = scan.raw_signals[key]*scan.monitor[:,None]
                    # signals
                    y1 = np.zeros((dim0, signals1.shape[1]))*np.nan
                    y2 = np.zeros((dim0, signals1.shape[1]))*np.nan
                    y1[0:len(signals1),:,:] = signals1
                    y2[0:len(signals2),:,:] = signals2
                    av_signals = np.nansum( (y1,y2) , axis=0 )
                    # errors of summed signals
                    av_errors  = np.sqrt( av_signals )
                    #
                    self.raw_signals[key] = av_signals/av_monitor[:,None]
                    self.raw_errors[key]  = av_errors/av_monitor[:,None]

        if interp=='Rbf':
            #rbfi = Rbf( scan.energy, scan.monitor, function='linear' )
            #self.monitor += rbgi( self.energy )
            for key in self.raw_signals:
                # remove zeros in errors
                self.raw_errors[key][self.raw_errors[key] == 0 ] = 1.0
                scan.raw_errors[key][scan.raw_errors[key] == 0 ] = 1.0

                # recover raw counts
                signals1 = self.raw_signals[key]*self.monitor
                signals2 = scan.raw_signals[key]*scan.monitor

                # sum interpolated signals
                rbfi = Rbf( scan.energy, signals2, function='linear' )
                av_signals = signals1 + rbfi( self.energy )

                # sum interpolated monitors
                rbfi = Rbf( scan.energy, scan.monitor, function='linear' )
                av_monitor = self.monitor + rbfi( self.energy )

                # errors of summed signals
                av_errors = np.sqrt(av_signals)

                # assign and renormalize
                self.raw_signals[key] = av_signals/av_monitor 
                self.raw_errors[key]  = av_errors/av_monitor  
                
    def get_type(self):
        """ **get_type**

        Returns the type of the scan.

        """
        return self.scan_type

    def get_scan_number(self):
        """ **get_scan_number**

        Returns the number of the scan.

        """
        return self.scan_number

    def get_shape(self):
        """ **get_shape**

        Returns the shape of the matrix holding the signals.

        """
        if not np.any(self.signals):
            print( 'please apply the ROIs first.' )
            return
        else:
            return self.signals.shape

    def get_num_of_rois(self):
        """ **get_num_of_rois**

        Returns the number of ROIs applied to the scan.

        """
        if not self.signals.any():
            print( 'please apply the ROIs first.' )
            return
        else:
            return self.signals.shape[1]

    def get_E_start(self):
        """ **get_E_start**

        Returs the first energy value.

        """
        return self.energy[0]

    def get_E_end(self):
        """ **get_E_end**

        Returs the last energy value.

        """
        return self.energy[-1]

    def get_resolution( self, keV2eV=True ):
        """ **get_resolution**

        Returns the ROI-wise resolution based on the 
        xrs_utilities fwhm method.
        """
        resolutions = []
        for ii in range(self.signals.shape[1]):
            x = self.energy
            y = self.signals[:,ii]
            if keV2eV:
                resolutions.append( xrs_utilities.fwhm(x,y)[0]*1.0e3 )
            else:
                resolutions.append( xrs_utilities.fwhm(x,y)[0] )
        return np.array(resolutions), np.mean(resolutions), np.std(resolutions)


class scan:
    """
    Container class, holding information of single scans performed with 2D detectors. 
    """
    def __init__(self,detector_arrays,scannumber,energy_scale,monitor_signal,counters,motor_positions,scan_data,scantype='generic'):
        # rawdata
        self.edfmats  = np.array(detector_arrays)     # numpy array of all 2D images that belong to the scan
        self.number   = scannumber                    # scan number 
        self.scantype = scantype                      # keyword, later used to group scans (add similar scans, etc.)
        self.energy   = np.array(energy_scale)        # numpy array of the energy axis
        self.monitor  = np.array(monitor_signal)      # numpy array of the monitor signal
        # some things maybe not imediately necessary 
        self.counters = counters                      # names of all counters for this scan (maybe unnecessary)
        self.motors   = motor_positions               # all motor positions for this scan
        self.scan_data = np.array(scan_data)          # all tabular data for this scan
        # data (to be filled after defining rois)
        self.eloss    = []     # numpy array of the energy loss scale
        self.signals  = []     # numpy array of signals extracted from the ROIs
        self.errors   = []     # numpy array of all Poisson errors
        self.cenom    = []     # list with center of masses (used if scan is an elastic line scan)
        # pixel-wise things
        self.signals_pw = []   
        self.cenom_pw   = []
        self.signals_pw_interp = []

    def applyrois(self,indices,scaling=None):
        """
        Sums up intensities found in the ROIs of each detector image
        and stores it into the self.signals attribute.
        roi_object = instance of the 'rois' class redining the ROIs
        scaling    = numpy array of numerical scaling factors (has to be one for each ROIs)
        """
        data = np.zeros((len(self.edfmats),len(indices)))
        for n in range(len(indices)): # each roi
            for m in range(len(self.edfmats)): # each energy point along the scan
                for l in range(len(indices[n])): # each pixel on the detector
                    data[m,n] += self.edfmats[m,indices[n][l][0],indices[n][l][1]]
        self.signals = np.array(data)
        self.errors  = np.sqrt(data)
        if np.any(scaling):
            assert len(scaling) == len(indices) # make sure, there is one scaling factor for each roi
            for ii in range(len(indices)):
                self.signals[:,ii] *= scaling[ii]
                self.errors[:,ii]  *= scaling[ii]

    def applyrois_pw(self,indices,scaling=None):
        """
        Pixel-wise reading of the ROI's pixels into a list of arrays. I.e. each n-pixel ROI will have n Spectra, saved in a 2D array.
        Parameters
        ----------
        indices : list
            List of indices (attribute of the xrs_rois class).
        scaling : list of flaots, optional
            Python list of scaling factors (one per ROI defined) to be applied to all pixels of that ROI.

        """
        data = [] # list of 2D arrays (energy vs. intensity for each pixel inside a single ROI) 
        for n in range(len(indices)): # each ROI
            roidata = np.zeros((len(self.edfmats),len(indices[n]))) # 2D np array energy vs pixels in current roi
            for m in range(len(self.edfmats)): # each energy point along the scan
                for l in range(len(indices[n])): # each pixel on the detector
                    roidata[m,l] = self.edfmats[m,indices[n][l][0],indices[n][l][1]]
            data.append(roidata) # list which contains one array (energy point, pixel) per ROI
        self.signals_pw = data

    def get_eloss_pw(self):
        """
        Finds the center of mass for each pixel in each ROI, sets the energy loss scale in 
        and interpolates the signals to a common energy loss scale. Finds the resolution (FWHM) for
        each pixel.

        """
        if not self.scantype == 'elastic': # return, if scantype is not elastic
            print( 'This method is meant for elastic line scans only!' )
            return
        if not self.signals_pw: # return, if there is no data
            print( 'Please use the applyrois_pw function first!' )
            return
        else:
            cenom_pw      = []
            resolution_pw = []
            for roiind in range(len(self.signals_pw)): # each ROI
                oneroi_cenom      = []
                oneroi_resolution = []
                for pixelind in range(len(self.signals_pw[roiind])): # each pixel in the ROI
                    oneroi_cenom.append(xrs_utilities.find_center_of_mass(self.energy,self.signals_pw[roiind][:,pixelind]))
                    try: 
                        FWHM,x0 = xrs_utilities.fwhm((self.energy - oneroi_cenom[roiind][pixelind])*1e3,self.signals_pw[roiind][:,pixelind])
                        oneroi_resolution.append(FWHM)
                    except:
                        oneroi_resolution.append(0.0)
                cenom_pw.append(oneroi_cenom)
                resolution_pw.append(oneroi_resolution)
            # define mean of first ROI as 'master' energy loss scale for all Pixels
            self.eloss = (self.energy - np.mean(cenom_pw[0]))*1e3 # average energy loss in eV from first ROI
            # !!!! mayb it is better to just save the CENOM values and to the interpolation later with the stiched 
            # !!!! whole spectrum... 
            # define eloss-scale for each ROI and interpolate onto the 'master' eloss-scale
            #for roiind in range(len(self.signals_pw)): # each ROI
            #   oneroi_signals = np.zeros_like(self.signals_pw[roiind])
            #   for pixelind in range(len(self.signals_pw[roiind])):
            #       x = (self.energy-self.cenom_pw[roiind][pixelind])*1e3
            #       y =
            #       f = interp1d(x,y, bounds_error=False,fill_value=0.0)
            #       self.signals[:,n] = f(self.eloss)
            #self.signals_pw_interp = []

    def get_type(self):
        return self.scantype

    def get_scannumber(self):
        return self.number

    def get_shape(self):
        if not np.any(self.signals):
            print( 'please apply the ROIs first.' )
            return
        else:
            return np.shape(self.signals)

    def get_numofrois(self):
        if not self.signals.any():
            print( 'please apply the ROIs first.' )
            return
        else:
            return np.shape(self.signals)[1]

class scangroup:
    """
    Container class holding information from a group of scans.
    """
    def __init__(self,energy,signals,errors,grouptype='generic'):
        self.energy    = energy     # common energy scale
        self.eloss     = []         # common energy loss scale
        self.signals   = signals    # summed up signals
        self.errors    = errors     # Poisson errors
        self.grouptype = grouptype  # keyword common to all scans
        self.signals_orig = signals # keep a copy of uninterpolated data

    def get_type(self):
        return self.grouptype

    def get_meanenergy(self):
        return np.mean(self.energy)

    def get_estart(self):
        return self.energy[0]

    def get_eend(self):
        return self.energy[-1]

    def get_meanegridspacing(self):
        return np.mean(np.diff(self.energy))

    def get_maxediff(self):
        return (self.energy[-1]-self.energy[0])

class Scan_group:
    """
    Container class holding information from a group of scans.
    """
    def __init__(self, energy, signals, errors, group_type='generic'):
        self.energy    = energy     # common energy scale
        self.eloss     = []         # common energy loss scale
        self.signals   = signals    # summed up signals
        self.errors    = errors     # Poisson errors
        self.grouptype = group_type  # keyword common to all scans
        self.raw_signals = {}
        self.raw_errors  = {}

    def get_type(self):
        return self.grouptype

    def get_meanenergy(self):
        return np.mean(self.energy)

    def get_estart(self):
        return self.energy[0]

    def get_eend(self):
        return self.energy[-1]

    def get_meanegridspacing(self):
        return np.mean(np.diff(self.energy))

    def get_maxediff(self):
        return (self.energy[-1]-self.energy[0])


class offDiaDataSet:
    """ **offDiaDataSet**
    Class to hold information from an off-diagonal dataset.
    """
    def __init__(self):
        self.RCmonitor    = np.array([])
        self.signalMatrix = np.array([])
        self.errorMatrix  = np.array([])
        self.motorMatrix  = np.array([])
        self.I0Matrix     = np.array([])
        self.energy       = np.array([])
        self.eloss        = np.array([])
        self.ROI_number   = 0
        self.G_vector     = np.array([])
        self.q0           = np.array([])
        self.qh           = np.array([])
        self.k0           = np.array([])
        self.kh           = np.array([])
        self.kprime       = np.array([])
        self.alignedSignalMatrix = np.array([])
        self.alignedErrorMatrix  = np.array([])
        self.alignedRCmonitor    = np.array([])
        self.masterRCmotor= np.array([])

    def filterDetErrors(self,threshold=3000000):
        inds = np.where(self.signalMatrix >= threshold)
        for ii in range(len(inds[0])):
            self.signalMatrix[inds[0][ii], inds[1][ii]] = 0.0
            self.signalMatrix[inds[0][ii], inds[1][ii]] = np.interp(inds[0][ii], [inds[0][ii]-1,inds[0][ii]+1] , [self.signalMatrix[inds[0][ii], inds[1][ii]-1],self.signalMatrix[inds[0][ii], inds[1][ii]+1]])

    def normalizeSignals(self):
        self.signalMatrix /= self.I0Matrix
        self.errorMatrix  /= self.I0Matrix
        self.RCmonitor    /= self.I0Matrix
        self.signalMatrix *= np.mean(self.I0Matrix)
        self.errorMatrix  *= np.mean(self.I0Matrix)
        self.RCmonitor    *= np.mean(self.I0Matrix)

    def alignRCmonitor(self):
        # check if data exists
        if not np.any(self.RCmonitor):
            print('Please load some data first.')
            return
        #x = np.arange(len(self.motorMatrix.T[:,0]))

        #RCmonitor   = self.RCmonitor#[:,diagonal_inds[0]:-diagonal_inds[1]]
        #motorMatrix = self.motorMatrix#[:,diagonal_inds[0]:-diagonal_inds[1]]
        #signalMatrix= self.signalMatrix#[:,diagonal_inds[0]:-diagonal_inds[1]]

        RCposition = []
        RCmax = []
        for ii in range(len(self.RCmonitor)):
            x = self.motorMatrix[ii,:]
            y = self.RCmonitor[ii,:]
            try:
                guess = [x[np.where(y == np.amax(y))[0]][0], 0.01, 1.0, np.amax(y), 1.]
                popt, pcov = optimize.curve_fit(math_functions.pearson7_forcurvefit, x, y,p0=guess) 
                RCposition.append(popt[0])
            except: 
                RCposition.append(xrs_utilities.find_center_of_mass(x,y))
            RCmax.append(np.amax(y))

        RCmax = np.array(RCmax)
        RCmax /= np.mean(RCmax)

        # possibly correct for deviations from Bragg's law
        #RCfit = np.polyval(np.polyfit(self.energy,self.RCposition),self.energy)
        #for ii in range(len(RCfit)):

        master_phi = self.motorMatrix[10,:]-RCposition[10]
        signalMatrix = np.zeros((len(self.energy),len(master_phi)))
        errorMatrix  = np.zeros((len(self.energy),len(master_phi)))
        RCmonitor    = np.zeros((len(self.energy),len(master_phi)))
        for ii in range(len(self.energy)):
            signalMatrix[ii,:] = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.signalMatrix[ii,:])*RCmax[ii]
            errorMatrix[ii,:]  = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.errorMatrix[ii,:])*RCmax[ii]
            RCmonitor[ii,:]    = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.RCmonitor[ii,:])*RCmax[ii]

        self.alignedSignalMatrix = signalMatrix
        self.alignedErrorMatrix  = errorMatrix
        self.masterRCmotor       = master_phi
        self.alignedRCmonitor    = RCmonitor

    def alignRCmonitor2(self):
        # check if data exists
        if not np.any(self.RCmonitor):
            print('Please load some data first.')
            return
        #x = np.arange(len(self.motorMatrix.T[:,0]))

        #RCmonitor   = self.RCmonitor#[:,diagonal_inds[0]:-diagonal_inds[1]]
        #motorMatrix = self.motorMatrix#[:,diagonal_inds[0]:-diagonal_inds[1]]
        #signalMatrix= self.signalMatrix#[:,diagonal_inds[0]:-diagonal_inds[1]]

        RCposition = []
        for ii in range(len(self.RCmonitor)):
            x = self.motorMatrix[ii,:]
            y = self.RCmonitor[ii,:]
            #try:
            #    guess = [x[np.where(y == np.amax(y))[0]][0], 0.01, 1.0, np.amax(y), 1.]
            #    popt, pcov = optimize.curve_fit(math_functions.pearson7_forcurvefit, x, y,p0=guess) 
            #    RCposition.append(popt[0])
            #except: 
            RCposition.append(xrs_utilities.find_center_of_mass(x,y))

        # possibly correct for deviations from Bragg's law
        #RCfit = np.polyval(np.polyfit(self.energy,self.RCposition),self.energy)
        #for ii in range(len(RCfit)):

        master_phi = self.motorMatrix[10,:]-RCposition[10]
        signalMatrix = np.zeros((len(self.energy),len(master_phi)))
        errorMatrix  = np.zeros((len(self.energy),len(master_phi)))
        RCmonitor    = np.zeros((len(self.energy),len(master_phi)))
        for ii in range(len(self.energy)):
            signalMatrix[ii,:] = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.signalMatrix[ii,:])
            errorMatrix[ii,:]  = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.errorMatrix[ii,:])
            RCmonitor[ii,:]    = np.interp(master_phi,self.motorMatrix[ii,:]-RCposition[ii],self.RCmonitor[ii,:])

        self.alignedSignalMatrix = signalMatrix
        self.alignedErrorMatrix  = errorMatrix
        self.masterRCmotor       = master_phi
        self.alignedRCmonitor    = RCmonitor

    def alignRCmonitorCC(self,repeat=2):
        """ **alignRCmonitorCC**
        Use cross-correlation to align data matrix according to the Rockin-Curve monitor.
        """
        # check if data exists
        if not np.any(self.RCmonitor):
            print('Please load some data first.')
            return

        signalMatrix = np.zeros((len(self.energy),len(self.motorMatrix[0,:])))
        errorMatrix  = np.zeros((len(self.energy),len(self.motorMatrix[0,:])))
        RCmonitor    = np.zeros((len(self.energy),len(self.motorMatrix[0,:])))

        # first iteration
        for ii in range(len(self.RCmonitor)):
            x0 = self.RCmonitor[0,:]
            x  = self.RCmonitor[ii,:]
            y  = np.correlate(x0,x,mode='same')
            ind = np.where(y == np.amax(y))[0]
            if ii == 0:
                master_phi = self.motorMatrix[ii,:] - self.motorMatrix[ii,ind]

            signalMatrix[ii,:] = np.interp(master_phi,self.motorMatrix[ii,:]-self.motorMatrix[ii,ind],self.signalMatrix[ii,:])
            errorMatrix[ii,:]  = np.interp(master_phi,self.motorMatrix[ii,:]-self.motorMatrix[ii,ind],self.errorMatrix[ii,:])
            RCmonitor[ii,:]    = np.interp(master_phi,self.motorMatrix[ii,:]-self.motorMatrix[ii,ind],self.RCmonitor[ii,:])

        # further iterations
        if repeat:
            for jj in range(repeat):
                for ii in range(len(RCmonitor)):
                    x0 = RCmonitor[0,:]
                    x  = RCmonitor[ii,:]
                    y  = np.correlate(x0,x,mode='same')
                    ind = np.where(y == np.amax(y))[0]
                    signalMatrix[ii,:] = np.interp(master_phi,master_phi-master_phi[ind],self.signalMatrix[ii,:])
                    errorMatrix[ii,:]  = np.interp(master_phi,master_phi-master_phi[ind],self.errorMatrix[ii,:])
                    RCmonitor[ii,:]    = np.interp(master_phi,master_phi-master_phi[ind],self.RCmonitor[ii,:])

        self.alignedSignalMatrix = signalMatrix
        self.alignedErrorMatrix  = errorMatrix
        self.masterRCmotor       = master_phi
        self.alignedRCmonitor    = RCmonitor

    def deglitchSignalMatrix(self,startpoint,stoppoint,threshold):
        signalMatrix = self.alignedSignalMatrix
        for ii in range(signalMatrix.shape[1]):
            ind = np.where(signalMatrix[startpoint:stoppoint,ii] >= threshold)[0]
            if np.any(ind):
                signalMatrix[startpoint+ind, ii] = (signalMatrix[startpoint+ind-1,ii] + signalMatrix[startpoint+ind+1,ii] )/2.0
        self.alignedSignalMatrix = signalMatrix

    def interpolateMatrix(self,master_matrix,master_energy,master_RCmotor):
        from scipy import interpolate
        signalMatrix = self.alignedSignalMatrix
        y=self.energy
        x=self.masterRCmotor
        xx, yy = np.meshgrid(x, y)
        z = signalMatrix
        f = interpolate.interp2d(x, y, z, kind='cubic')
        ynew = master_energy
        xnew = master_RCmotor
        znew = f(xnew, ynew)
        self.interpSignalMatrix = znew
        self.interpEnergy  = ynew
        self.interpRCmotor = xnew

    def removeElastic(self,fitrange=[-6.0,2.0]):
        self.energy = np.array(self.energy)
        cenom = []
        for ii in range(self.signalMatrix.shape[1]):
            inds = np.where(np.logical_and(np.array(self.energy) >= fitrange[0], np.array(self.energy)<= fitrange[1]))[0]
            x = self.energy[inds]
            y = self.alignedSignalMatrix[inds,ii]
            guess = [x[np.where(y==np.amax(y))[0]], 0.002, 100.0,  np.amax(y) ,0.0]
            fitfct = lambda a: np.sum( (y - math_functions.pearson7(x,a) )**2.0 )
            params = optimize.minimize(fitfct, guess, method='SLSQP').x
            cenom.append(params[0])
            back = math_functions.pearson7(self.energy,params )
            plt.ion()
            plt.cla()
            plt.plot(self.energy,self.alignedSignalMatrix[:,ii],self.energy,back,self.energy,self.alignedSignalMatrix[:,ii]-back)
            plt.waitforbuttonpress()
            self.alignedSignalMatrix[:,ii] = self.alignedSignalMatrix[:,ii]-back
        self.eloss = (self.energy - np.mean(cenom))*1.0e3

    def removeElastic2(self, fitrange1, fitrange2, guess=None):
        """ **removeElastic2**
        Subtract Pearson7 plus linear.
        """
        self.alignedSignalMatrixB = np.zeros_like(self.alignedSignalMatrix)
        self.energy = np.array(self.energy)
        cenom = []
        for ii in range(self.signalMatrix.shape[1]):
            region1 = np.where(np.logical_and(np.array(self.energy) >= fitrange1[0], np.array(self.energy) <= fitrange1[1]))[0]
            region2 = np.where(np.logical_and(np.array(self.energy) >= fitrange2[0], np.array(self.energy) <= fitrange2[1]))[0]
            inds    = np.append(region1,region2)
            x = self.energy[inds]
            y = self.alignedSignalMatrix[inds,ii]
            if not guess:
                guess = [x[np.where(y==np.amax(y))[0]], 0.001, 1.0,  np.amax(y) ,0.0, -0.1, 0.01]
            fitfct = lambda a: np.sum( (y - math_functions.pearson7_linear(x,a) )**2.0 )
            params = optimize.minimize(fitfct, guess, method='COBYLA',tol=1e-20).x
            cenom.append(params[0])
            back = math_functions.pearson7_linear(self.energy,params )
            plt.ion()
            plt.cla()
            plt.plot(self.energy,self.alignedSignalMatrix[:,ii],self.energy,back,self.energy,self.alignedSignalMatrix[:,ii]-back)
            plt.waitforbuttonpress()
            self.alignedSignalMatrixB[:,ii] = self.alignedSignalMatrix[:,ii]-back
        self.eloss = (self.energy - np.mean(cenom))*1.0e3

    def windowSignalMatrix(self,estart,estop):
        inds = np.where(np.logical_and(self.eloss>=estart, self.eloss<= estop))[0]
        self.alignedSignalMatrix = self.alignedSignalMatrix[inds[0]:inds[-1],:]

    def removeLinearBack(self,fitrange1,fitrange2):
        self.energy = np.array(self.energy)
        for ii in range(self.signalMatrix.shape[1]):
            region1 = np.where(np.logical_and(np.array(self.energy) >= fitrange1[0], np.array(self.energy) <= fitrange1[1]))
            region2 = np.where(np.logical_and(np.array(self.energy) >= fitrange2[0], np.array(self.energy) <= fitrange2[1]))
            region  = np.append(region1,region2)
            x = np.array(self.energy[region])
            y = self.alignedSignalMatrix[region,ii]
            back = np.polyval(np.polyfit(x,y,1),self.energy)
            self.alignedSignalMatrix[:,ii] -= back

    def removeConstBack(self,fitrange1,fitrange2):
        self.energy = np.array(self.energy)
        for ii in range(self.signalMatrix.shape[1]):
            region1 = np.where(np.logical_and(np.array(self.energy) >= fitrange1[0], np.array(self.energy) <= fitrange1[1]))
            region2 = np.where(np.logical_and(np.array(self.energy) >= fitrange2[0], np.array(self.energy) <= fitrange2[1]))
            region  = np.append(region1,region2)
            x = np.array(self.energy[region])
            y = self.alignedSignalMatrix[region,ii]
            back = np.polyval(np.polyfit(x,y,0),self.energy)
            self.alignedSignalMatrix[:,ii] -= back

    def removePearsonBack(self,fitrange1,fitrange2):
        self.energy = np.array(self.energy)
        for ii in range(self.signalMatrix.shape[1]):
            region1 = np.where(np.logical_and(np.array(self.energy) >= fitrange1[0], np.array(self.energy) <= fitrange1[1]))
            region2 = np.where(np.logical_and(np.array(self.energy) >= fitrange2[0], np.array(self.energy) <= fitrange2[1]))
            region  = np.append(region1,region2)
            x = np.array(self.energy[region])
            y = self.alignedSignalMatrix[region,ii]
            guess = [x[np.where(y==np.amax(y))[0]], 0.002, 100.0,  np.amax(y) ,0.0]
            fitfct = lambda a: np.sum( (y - math_functions.pearson7(x,a) )**2.0 )
            params = optimize.minimize(fitfct, guess, method='SLSQP').x
            back = math_functions.pearson7(self.energy,params)
            self.alignedSignalMatrix[:,ii] -= back

    def replaceSignalByConstant(self,fitrange):
        self.energy = np.array(self.energy)
        for ii in range(self.signalMatrix.shape[1]):
            inds = np.where(np.logical_and(np.array(self.energy) >= fitrange[0], np.array(self.energy) <= fitrange[1]))[0]
            x = np.array(self.energy[inds])
            y = self.alignedSignalMatrix[inds,ii]
            back = np.polyval(np.polyfit(x,y,0),self.energy)
            self.alignedSignalMatrix[:,ii] = back

def findgroups(scans):
    """
    this groups together instances of the scan class based on their  "scantype" attribute and returns ordered scans
    """ 
    allscannames = []
    for scan in scans:
        print( scan )
        allscannames.append(scan)
    allscannames.sort() # 
    allscans = []
    for scan in allscannames:
         allscans.append(scans[scan])
    allscans = sorted(allscans,key=lambda x:x.get_type())
    rawgroups = []
    results = groupby(allscans,lambda x:x.get_type())
    print( 'The following scangroups were found:' )
    for key,thegroup in results:
        print( key )
        ls = -1
        thegroup = list(thegroup)
        for t in thegroup:
            if ls!=-1 and len(t.monitor)!=ls:
                print( " Scan Number :" +str( t.scan_number  )+" is added to group of key " +key + " but has lenght "+  str(len(t.monitor)) + "  versus " + str( ls)   )
        rawgroups.append(list(thegroup))
    return rawgroups

def makegroup(groupofscans,grouptype=None):
    """
    takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup        class (turns several instances of the "scan" class into an instance of the "scangroup" class)
    """
    if not grouptype:
        grouptype = groupofscans[0].get_type() # the type of the sum of scans will be the same as the first from the list   
    theenergy   = groupofscans[0].energy # all scans are splined onto energy grid of the first scan
    thesignals  = np.zeros(groupofscans[0].get_shape())
    theerrors   = np.zeros(groupofscans[0].get_shape())
    themonitors = np.zeros(np.shape(groupofscans[0].monitor))
    for scan in groupofscans:
        f  = interpolate.interp1d(scan.energy,scan.monitor, bounds_error=False, fill_value=0.0)
        moni = f(theenergy)
        themonitors += moni
        for n in range(thesignals.shape[-1]):
            f = interpolate.interp1d(scan.energy,scan.signals[:,n], bounds_error=False, fill_value=0.0)
            signal = f(theenergy)
            thesignals[:,n] += signal*moni
    for n in range(thesignals.shape[-1]):
        theerrors[:,n]  = np.sqrt(thesignals[:,n])
    # and normalize 
    for n in range(thesignals.shape[-1]):
        thesignals[:,n] = thesignals[:,n]/themonitors
        theerrors[:,n]  = theerrors[:,n]/themonitors

    group = scangroup(theenergy,thesignals,theerrors,grouptype)
    return group

def makegroup_nointerp(groupofscans,grouptype=None,absCounts=False):
    """
    takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup class (turns several instances of the "scan" class into an instance of the "scangroup" class), same as makegroup but withouth interpolation to account for encoder differences... may need to add some "linspace" function in case the energy scale is not monotoneous... 
    """
    if not grouptype:
        grouptype = groupofscans[0].get_type() # the type of the sum of scans will be the same as the first from the list   
    theenergy   = groupofscans[0].energy
    thesignals  = np.zeros(groupofscans[0].get_shape())
    theerrors   = np.zeros(groupofscans[0].get_shape())
    themonitors = np.zeros(np.shape(groupofscans[0].monitor))

    
    for scan in groupofscans:
        themonitors += scan.monitor

        for n in range(thesignals.shape[-1]):
            thesignals[:,n] += scan.signals[:,n]* scan.monitor
    for n in range(thesignals.shape[-1]):
        theerrors[:,n]  = np.sqrt(thesignals[:,n])
        
    # and normalize
    for n in range(thesignals.shape[-1]):
        thesignals[:,n] = thesignals[:,n]/themonitors
        theerrors[:,n]  = theerrors[:,n]/themonitors
         

    group = scangroup(theenergy,thesignals,theerrors,grouptype)
    return group

def make_scan_group_sum( group_of_scans, group_type=None, abs_counts=False):
    """ **make_scan_group_sum**

    Sums together a list of scans with equal scan_type and returns an instance 
    of the Scan_group class.

    Args:
        group_of_scans (list): List containing instances of the Scan class to be summed up.
        group_type      (str): Keyword defining the type of scans, if None (default) the 
            group_type will be the type of the first scan in the group_of_scans list.
        abs_counts  (boolean): Boolean defining if results should be returned in absolute
            count units (ct/s).
        time_counter    (str): Counter name for the counting time mnemonic.

    Returns:
        group (obj): Instance of the Scan_group container class.

    """
    # if not provided, type of returned group will be same as first in group_of_scans
    if not group_type:
        group_type = group_of_scans[0].get_type()   

    # create the arrays for energy, signals, and errors
    theenergy   = group_of_scans[0].energy
    thesignals  = np.zeros( group_of_scans[0].get_shape() )
    theerrors   = np.zeros( group_of_scans[0].get_shape() )
    themonitors = np.zeros( np.shape(group_of_scans[0].monitor) )

    
    # sum up the different scans in the list
    for scan in group_of_scans:
        themonitors += scan.monitor

        for n in range(thesignals.shape[-1]):
            thesignals[:,n] += scan.signals[:,n]*scan.monitor
            
    for n in range(thesignals.shape[-1]):
        theerrors[:,n]  = np.sqrt(thesignals[:,n])

    # and normalize by the sum of the monitor signal
    for n in range(thesignals.shape[-1]):
        thesignals[:,n] = thesignals[:,n]/themonitors
        theerrors[:,n]  = theerrors[:,n]/themonitors


            
    group = Scan_group( theenergy, thesignals, theerrors, group_type )

    return group

def make_scan_group_pixel( group_of_scans, group_type=None, abs_counts=False ):
    """ **make_scan_group_pixel**

    Pixel-by-pixel summation of a list of scans with equal scan_type and returns an instance 
    of the Scan_group class.

    Args:
        group_of_scans (list): List containing instances of the Scan class to be summed up.
        group_type      (str): Keyword defining the type of scans, if None (default) the 
            group_type will be the type of the first scan in the group_of_scans list.
        abs_Counts  (boolean): Boolean defining if results should be returned in absolute
            count units (ct/s).
        time_counter    (str): Counter name for the counting time mnemonic.

    Returns:
        group (obj): Instance of the Scan_group container class.

    """
    # if not provided, type of returned group will be same as first in group_of_scans
    if not group_type:
        group_type = group_of_scans[0].get_type()   

    # create the arrays/dicts for energy, signals, and errors
    energy   = group_of_scans[0].energy
    monitors = np.zeros(np.shape(groupofscans[0].monitor))
    raw_signals = {}
    raw_errors  = {}
    for key in group_of_scans[0].raw_signals:
        raw_signals[key] = np.zeros_like(group_of_scans[0].raw_signals[key])
        raw_errors[key]  = np.zeros_like(group_of_scans[0].raw_errors[key])

    # sum up the different scans in the list (pixel-by-pixel)
    for scan in group_of_scans:
        monitors += scan.monitor
        for key in scan.raw_signals:
            raw_signals[key] += scan.raw_signals[key]*scan.monitor

    for key in raw_signals:
        raw_errors[key] = np.sqrt(raw_signals[key])

    # and normalize by the sum of the monitor signal
    for key in raw_signals:
        raw_signals[key] /= monitors
        raw_errors[key]  /= monitors

    # create the group
    group = Scan_group( energy, np.array([]), np.array([]), group_type )
    group.raw_signals = raw_signals
    group.raw_errors  = raw_errors

    return group



def append2Scan_right(group1,group2,inds=None,grouptype='spectrum'):
    """
    append two instancees of the scangroup class, return instance of scangroup
    append group2[inds] to the right (higher energies) of group1
    if inds is not None, only append rows indicated by inds to the first group 
    """
    # assert isinstance(group1,scangroup) and isinstance(group2,scangroup)
    energy  = group1.energy
    signals = group1.signals
    errors  = group1.errors
    gtype   = group1.grouptype
    if not inds:
        energy  = np.append(energy,group2.energy)
        signals = np.append(signals,np.squeeze(group2.signals),0)
        errors  = np.append(errors,np.squeeze(group2.errors),0)
        return scangroup(energy,signals,errors,grouptype=gtype)
    else:
        energy  = np.append(energy,group2.energy[inds])
        signals = np.append(signals,np.squeeze(group2.signals[inds,:]),0)
        errors  = np.append(errors,np.squeeze(group2.errors[inds,:]),0)
        return scangroup(energy,signals,errors,grouptype=gtype)



def append2Scan_left(group1,group2,inds=None,grouptype='spectrum'):
    """
    append two instancees of the scangroup class, return instance of scangroup
    append group1[inds] to the left (lower energies) of group2
    if inds is not None, only append rows indicated by inds to the first group 
    """
    assert isinstance(group1,scangroup) and isinstance(group2,scangroup)
    if not inds:
        energy  = group1.energy
        signals = group1.signals
        errors  = group1.errors
        gtype   = group1.grouptype
        energy  = np.append(energy,group2.energy)
        signals = np.append(signals,np.squeeze(group2.signals),0)
        errors  = np.append(errors,np.squeeze(group2.errors),0)
        return scangroup(energy,signals,errors,grouptype=gtype)
    else:
        energy  = group1.energy[inds]
        signals = group1.signals[inds,:]
        errors  = group1.errors[inds,:]
        gtype   = group1.grouptype
        energy  = np.append(energy,group2.energy)
        signals = np.append(signals,np.squeeze(group2.signals),0)
        errors  = np.append(errors,np.squeeze(group2.errors),0)
        return scangroup(energy,signals,errors,grouptype=gtype)

def insertScan(group1,group2,grouptype='spectrum'):
    """
    inserts group2 into group1
    NOTE! there is a numpy insert function, maybe it would be better to use that one!
    """
    # find indices for below and above group2
    lowinds  = np.where(group1.energy<group2.get_estart())
    highinds = np.where(group1.energy>group2.get_eend())

    energy  = np.append(group1.energy[lowinds],np.append(group2.energy,group1.energy[highinds]))
    signals = np.append(np.squeeze(group1.signals[lowinds,:]),np.append(group2.signals,np.squeeze(group1.signals[highinds,:]),0),0)
    errors  = np.append(np.squeeze(group1.errors[lowinds,:]),np.append(group2.errors,np.squeeze(group1.errors[highinds,:]),0),0)

    return scangroup(energy,signals,errors,grouptype)

def catScansLong(groups,include_elastic):
    """
    takes a longscan and inserts other backgroundscans (scans that have 'long' in their name) and other scans and inserts them into the long scan. 
    """
    # the long scan
    spectrum = groups['long']

    # groups that don't have 'long' in the grouptype    
    allgroups  = []
    for group in groups:
        if not 'long' in group:
            allgroups.append(groups[group])
    allgroups.sort(key = lambda x:x.get_estart())

    # groups that have 'long' in the grouptype  
    longgroups = []
    for group in groups:
        if 'long' in group and group != 'long':
            longgroups.append(groups[group])
    longgroups.sort(key = lambda x:x.get_estart())

    # if there are other longscans: insert those first into the long scan
    for group in longgroups:
        spectrum = insertScan(spectrum,group)

    # insert other scans into the long scan
    for group in allgroups:
        spectrum = insertScan(spectrum,group)

    # cut off the elastic line if present in the groups
    if 'elastic' in groups and not include_elastic:
        inds = np.where(spectrum.energy > groups['elastic'].get_eend())[0]
        return spectrum.energy[inds], spectrum.signals[inds,:], spectrum.errors[inds,:]
    else:
        return spectrum.energy, spectrum.signals, spectrum.errors

def catScans(groups,include_elastic):
    """
    concatenate all scans in groups, return the appended energy, signals, and errors
    """
    # sort the groups by their start-energy
    allgroups  = []
    for group in groups:
        allgroups.append(groups[group])
    allgroups.sort(key = lambda x:x.get_estart())

    # assign first group to the spectrum (which is an instance of the scangroup class as well)
    spectrum = scangroup(allgroups[0].energy,allgroups[0].signals,allgroups[0].errors,grouptype='spectrum')

    # go through all other groups and append them to the right of the spectrum
    if len(allgroups)>1: # check if there are groups to append
        for group in allgroups[1:]:
            spectrum = append2Scan_right(spectrum,group)

    # cut off the elastic line if present in the groups
    if 'elastic' in groups and not include_elastic:
        inds = np.where(spectrum.energy > groups['elastic'].get_eend())[0]
        return spectrum.energy[inds], spectrum.signals[inds,:], spectrum.errors[inds,:]
    else:
        return spectrum.energy, spectrum.signals, spectrum.errors

def catScans_pixel( groups, include_elastic ):
    """ **catScans_pixel**

    Stitch together all scans in groups in a pixel-by-pixel fashion for 
    the case of no available long (overview) scan.

    Args:
        groups (list): List of scan-groups to be stitched together.
        include_elastic (boolean): Boolean switch if the elastic line
            should be included in the final spectrum or not.

    Returns:
        energy  (np.array): Array of energy loss values.
        raw_signals (dict): Dictionary (one entry per ROI) of 
            pixel-by-pixel intensities.
        raw_errors  (dict): Dictionary (one entry per ROI) of 
            pixel-by-pixel Poisson errors.

    """
    # sort the groups by their start-energy
    all_groups  = []
    for group in groups:
        all_groups.append(groups[group])
    all_groups.sort(key = lambda x:x.get_estart())

    # create a scan-group for the stitched spectrum
    spectrum = Scan_group( all_groups[0].energy, np.array([]), np.array([]), group_type='spectrum' )
    spectrum.raw_signals = all_groups[0].raw_signals
    spectrum.raw_errors  = all_groups[0].raw_errors

    # go through all other groups and append them to the right of the spectrum
    if len(all_groups)>1: # check if there are groups to append
        for group in all_groups[1:]:
            spectrum = append2Scan_right_pixel( spectrum, group )

    # cut off the elastic line if present in the groups
    if 'elastic' in groups and not include_elastic:
        inds = np.where(spectrum.energy > groups['elastic'].get_eend())[0]
        return spectrum.energy[inds], spectrum.signals[inds,:], spectrum.errors[inds,:]
    else:
        return spectrum.energy, spectrum.signals, spectrum.errors


def appendScans(groups,include_elastic):
    """
    try including different background scans... 
    append groups of scans ordered by their first energy value.
    long scans are inserted into gaps that at greater than two times the grid of the finer scans
    """
    # find all types of groups  
    grouptypes = [key for key in groups.keys()]

    if 'long' in grouptypes:
        print( " going to refine " )
        return catScansLong(groups,include_elastic)
    else: 
        return catScans(groups,include_elastic)

def appendScans_pixel( groups, include_elastic ):
    """ **appendScans_pixel**

    Decides if there is a long scan available and selects
    accordingly which stitching method to use.

    Args:
        groups (list): List of scan-groups to be stitched together.
        include_elastic (boolean): Boolean switch if the elastic line
            should be included in the final spectrum or not.

    Returns:
        energy  (np.array): Array of energy loss values.
        raw_signals (dict): Dictionary (one entry per ROI) of 
            pixel-by-pixel intensities.
        raw_errors  (dict): Dictionary (one entry per ROI) of 
            pixel-by-pixel Poisson errors.
        
    """
    # find all types of groups  
    grouptypes = [key for key in groups.keys()]

    if 'long' in grouptypes:
        print( " going to refine " )
        return catScansLong_pixel( groups, include_elastic )
    else: 
        return catScans_pixel( groups, include_elastic )


def stitch_groups_to_spectrum(groups, method='sum', include_elastic=False ):
    """ **stitch_groups_to_spectrum**

    Takes a dictionary of instances of the Scan class and stitches them
    together to produce a spectrum. Long scans and scans that have 'long'
    in ther scan_type attribute are treated specially.

    Args:
        groups (list): List of instances of the Scan class.
        method  (str): Keyword describing the kind of integration scheme
            to be used (possible values are 'sum', 'pixel', or 'row'), 
            default is 'sum'.
        include_elastic (boolean): Boolean flag deciding if the elastic
            should be included in the final spectrum.

    Returns:
        Scan class instance with the stitched spectrum.

    """

    # create the final spectrum
    spectrum = Scan()

    # find all groups that are not long scans, sort by acending energy 
    all_groups  = []
    for group in groups:
        if not 'long' in group:
            all_groups.append(groups[group])
    all_groups.sort( key = lambda x:x.get_E_start() )

    # groups that have 'long' in the grouptype  
    long_groups = []
    for group in groups:
        # if there is a long scan
        if group == 'long':
            spectrum.energy  = groups[group].energy
            spectrum.monitor = groups[group].monitor
            #if method is 'sum':
            #   spectrum.signals = groups[group].signals
            #   spectrum.errors  = groups[group].errors
            #elif method is 'pixel' or 'row':
            spectrum.raw_signals = groups[group].raw_signals
            spectrum.raw_errors  = groups[group].raw_errors
            #else:
            #   print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
            #   return
        # if there is a scan that has 'long' in its name
        if 'long' in group and group != 'long':
            long_groups.append(groups[group])
            long_groups.sort( key = lambda x:x.get_E_start() )

    # if long exists, insert backround groups, then other scans
    if np.any(spectrum.energy):
        for group in long_groups:
            spectrum.insert_scan( group, method=method )
        for group in all_groups:
            spectrum.insert_scan( group, method=method )

    # if no long scan exists, just append all others
    else:
        spectrum.energy  = all_groups[0].energy
        spectrum.monitor = all_groups[0].monitor
        #if method is 'sum':
        spectrum.raw_signals = all_groups[0].raw_signals
        spectrum.raw_errors  = all_groups[0].raw_errors
        #elif method is 'pixel' or 'row':
        #   spectrum.raw_signals = all_groups[0].raw_signals
        #   spectrum.raw_errors  = all_groups[0].raw_errors
        #else:
        #   print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
        #   return
        for group in all_groups[1::]:
            spectrum.append_scan(group, method=method)

    # cut elastic line if applicable
    if 'elastic' in groups and not include_elastic:
        
        inds = np.where(spectrum.energy > groups['elastic'].get_E_end())[0]
        spectrum.energy  = spectrum.energy[inds]
        spectrum.monitor = spectrum.monitor[inds]
        for key in spectrum.raw_signals:
            if method == 'sum':
                spectrum.raw_signals[key] = spectrum.raw_signals[key][inds]
                spectrum.raw_errors[key]  = spectrum.raw_errors[key][inds]
            if method in [ 'pixel' , 'pixel2']:
                spectrum.raw_signals[key] = spectrum.raw_signals[key][inds,:,:]
                spectrum.raw_errors[key]  = spectrum.raw_errors[key][inds,:,:]
            if method == 'row':
                spectrum.raw_signals[key] = spectrum.raw_signals[key][inds,:]
                spectrum.raw_errors[key]  = spectrum.raw_errors[key][inds,:]    

    return spectrum

def get_XES_spectrum( groups ):
    """ **get_XES_spectrum**

    Constructs a XES spectrum from the given scan groups (sums of 
    separate emission scans). 

    Args:
        groups (dict): Dictionary of groups with partial XES scans. 

    """
    # sort the groups by their end-energy (is the smallest in XES/energy2-scans)
    allgroups  = []
    for group in groups:
        if groups[group].energy[0] > groups[group].energy[-1]:
            groups[group].energy = np.flipud( groups[group].energy )
            for ii in range(groups[group].signals.shape[1]):
                groups[group].signals[:,ii] = np.flipud( groups[group].signals[:,ii] )
        allgroups.append(groups[group])
    allgroups.sort(key = lambda x:x.get_E_end())

    # find lowest energy, highest energy, smallest increment, define grid
    eStart = np.amin([   np.amin(allgroups[ii].energy) for ii in range(len(allgroups))  ])
    eStop  = np.amax([   np.amax(allgroups[ii].energy) for ii in range(len(allgroups))  ])
    eStep  = np.amin(np.abs(  [np.diff(group.energy)[0] for group in allgroups]  ))

    energy  = np.arange(eStart,eStop,eStep)
    signals = np.zeros((len(energy),group.signals.shape[1],len(allgroups)))
    errors  = np.zeros((len(energy),group.signals.shape[1],len(allgroups)))

    # interpolate all groups onto grid, put into a matrix
    for group,ii in zip(allgroups,range(len(allgroups))):
        for jj in range(group.signals.shape[1]):
            interp_signals   = np.interp(energy, group.energy, group.signals[:,jj], left=0.0, right=0.0)
            signals[:,jj,ii] = interp_signals
            interp_errors    = np.interp(energy, group.energy, group.errors[:,jj], left=0.0, right=0.0)
            errors[:,jj,ii]  = interp_errors


    # sum up and weigh by number of available non-zero points
    sum_signals = np.zeros( (len(energy), group.signals.shape[1]))
    sum_errors  = np.zeros( (len(energy), group.signals.shape[1]))
    for jj in range(group.signals.shape[1]):
        for ii in range(len(energy)):
            nParts = len(np.where(signals[ii,jj]>0.0)[0])
            sum_signals[ii,jj] = np.sum(signals[ii,jj])/nParts
            sum_errors[ii,jj]  = np.sqrt(np.sum(errors[ii,jj]**2.0))/nParts

    spectrum = scangroup(energy, sum_signals, sum_errors, grouptype='spectrum')
    return spectrum.energy, spectrum.signals, spectrum.errors




def catXESScans(groups):
    """
    Concatenate all scans in groups, return the appended energy, signals, and errors.
    This needs to be a bit smarter to also work for scans that are scanned from small to large energy...
    """
    # sort the groups by their end-energy (is the smallest in XES/energy2-scans)
    allgroups  = []
    for group in groups:
        if groups[group].energy[0] > groups[group].energy[-1]:
            groups[group].energy = np.flipud( groups[group].energy )
            for ii in range(groups[group].signals.shape[1]):
                groups[group].signals[:,ii] = np.flipud( groups[group].signals[:,ii] )
        allgroups.append(groups[group])
    allgroups.sort(key = lambda x:x.get_eend())

    # find lowest energy, highest energy, smallest increment, define grid
    eStart = np.amin([   np.amin(allgroups[ii].energy) for ii in range(len(allgroups))  ])
    eStop  = np.amax([   np.amax(allgroups[ii].energy) for ii in range(len(allgroups))  ])
    eStep  = np.amin(np.abs(  [np.diff(group.energy)[0] for group in allgroups]  ))

    energy  = np.arange(eStart,eStop,eStep)
    signals = np.zeros((len(energy),group.signals.shape[1],len(allgroups)))
    errors  = np.zeros((len(energy),group.signals.shape[1],len(allgroups)))

    # interpolate all groups onto grid, put into a matrix
    for group,ii in zip(allgroups,range(len(allgroups))):
        for jj in range(group.signals.shape[1]):
            interp_signals   = np.interp(energy, group.energy, group.signals[:,jj], left=0.0, right=0.0)
            signals[:,jj,ii] = interp_signals
            interp_errors    = np.interp(energy, group.energy, group.errors[:,jj], left=0.0, right=0.0)
            errors[:,jj,ii]  = interp_errors


    # sum up and weigh by number of available non-zero points
    sum_signals = np.zeros( (len(energy), group.signals.shape[1]))
    sum_errors  = np.zeros( (len(energy), group.signals.shape[1]))
    for jj in range(group.signals.shape[1]):
        for ii in range(len(energy)):
            nParts = len(np.where(signals[ii,jj]>0.0)[0])
            sum_signals[ii,jj] = np.sum(signals[ii,jj])/nParts
            sum_errors[ii,jj]  = np.sqrt(np.sum(errors[ii,jj]**2.0))/nParts

    spectrum = scangroup(energy,sum_signals,sum_errors,grouptype='spectrum')
    return spectrum.energy, spectrum.signals, spectrum.errors

def appendXESScans(groups):
    """
    try including different background scans... 
    append groups of scans ordered by their first energy value.
    long scans are inserted into gaps that at greater than two times the grid of the finer scans
    """
    # find all types of groups  
    grouptypes = [key for key in groups.keys()]
    return catXESScans(groups)

def create_sum_image(scans,scannumbers):
    """
    Returns a summed image from all scans with numbers 'scannumbers'.
    scans       = dictionary of objects from the scan-class
    scannumbers = single scannumber, or list of scannumbers from which an image should be constructed
    """
    # make 'scannumbers' iterable (even if it is just an integer)
    numbers = []
    if not isinstance(scannumbers,list):
        numbers.append(scannumbers)
    else:
        numbers = scannumbers

    key = 'Scan%03d' % numbers[0]
    image = np.zeros_like(scans[key].edfmats[0,:,:])
    for number in numbers:
        key = 'Scan%03d' % number
        for ii in range(scans[key].edfmats.shape[0]):
            image += scans[key].edfmats[ii,:,:]

    return image


def create_diff_image(scans,scannumbers,energy_keV):
    """
    Returns a summed image from all scans with numbers 'scannumbers'.
    scans       = dictionary of objects from the scan-class
    scannumbers = single scannumber, or list of scannumbers from which an image should be constructed
    """
    # make 'scannumbers' iterable (even if it is just an integer)
    numbers = []
    if not isinstance(scannumbers,list):
        numbers.append(scannumbers)
    else:
        numbers = scannumbers

    key = 'Scan%03d' % numbers[0]
    below_image = np.zeros_like(scans[key].edfmats[0,:,:])
    above_image = np.zeros_like(scans[key].edfmats[0,:,:])

    # find indices below and above 'energy'
    below_inds = scans[key].energy < energy_keV
    above_inds = scans[key].energy > energy_keV
    for number in numbers:
        key = 'Scan%03d' % number
        for ii in below_inds:
            below_image += scans[key].edfmats[ii,:,:]
        for ii in above_inds:
            above_image += scans[key].edfmats[ii,:,:]

    return (above_image - below_image)

def findRCscans(scans):
    """ **findRCscans**
    Returns a list of scans with name RC.
    """
    RCscans = []
    for key in scans:
        if scans[key].get_type() == 'RC':
            RCscans.append(scans[key])
    return RCscans

def sum_scans_to_group( group, method='sum', interp=False ):
    """ **sum_scans_to_group**

    Sums up all scans in the list group to form a scan-group.

    Args:
        group (list): List of scans to be added up. 
        method (str): Keyword describing which data analysis method
            should be used. Possible values are 'sum', 'pixel', 'row'.
            Default is 'sum'.
        interp (boolean): Flag if interpolation onto the energy grid 
            of the first scan should be applied (default is False).

    Returns:
        An instance of the Scan class containing the summed up results.

    """

    # initialize result
    summed_group = Scan()
    summed_group.energy  = group[0].energy
    summed_group.monitor = group[0].monitor

    # sum up
    summed_group.raw_signals = group[0].raw_signals
    summed_group.raw_errors  = group[0].raw_errors
    for scan in group[1::]:
        summed_group.add_scan( scan, method=method, interp=interp )

    return summed_group
