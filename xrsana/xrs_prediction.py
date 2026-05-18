#!/usr/bin/python
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

import argparse
import ast
import io
import json
import os
import shelve
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

os.environ.setdefault('MPLCONFIGDIR', os.path.join('/tmp', 'xrsana-matplotlib'))

import h5py
import numpy as np
import pylab

from . import xrs_utilities
from scipy import constants, optimize

installation_dir = os.path.dirname(os.path.abspath(__file__))

def cla():
    pylab.cla()

class detector:
    """
    Class to describe detector related things. All default values are meant
    for the ESRF MAXIPIX detector.
    """

    def __init__(self, energy=9.68, thickness=500, material='Si', pixel_size=[256,768]):
        self.energy     = np.array(energy)    # analyzer energy [keV]
        self.thickness  = np.array(thickness) # thickness of the active material [microns]
        self.material   = material            # detector active material
        self.efficiency = []

    def set_energy(self,energy):
        self.energy = energy
    def get_energy(self):
        if np.any(self.energy):
            return self.energy
        else:
            print( 'No energy set, please set an enegy first!')
            return
    def set_thickness(self,thickness):
        self.thickness = thickness
    def get_thickness(self):
        if np.any(self.thickness):
            return self.thickness
        else:
            print( 'No thickness set, please set a thickness first!')
    def set_material(self,material):
        if isinstance(material,str):
            self.material = material.upper()
        else:
            print( 'material must be passed as a string. Default is \'Si\'!')
    def get_material(self):
        if any(self.material):
            return self.material
        else:
            print( 'No material set, please set the material first!')
    def set_size(self,size):
        if np.shape(np.array(size)) == (2,):
            self.size = np.array(size)
        else:
            print( 'size must be passed as a 2x0 numpy array or list of two entries.')
    def get_size(self):
        if np.any(self.size):
            return self.size
        else:
            print( 'No size has been set')
    def get_efficiency(self,energy=None):
        """
        calculates the detector efficiency at the given energy (simply given by 
        the absorption of the detector active material).
        """
        if not energy:
            energy = self.energy
        thickness = self.thickness*1e-4 # conversion to cm (needed for mpr routine)
        material  = self.material
        murho,rhov,mv = xrs_utilities.mpr(energy,material)
        return 1.0 - np.exp(-thickness*murho)

class analyzer:
    """
    Class to describe things related to the analyzer crystal used. Default values are for a Si(660) crystal.
    """
    
    def __init__(self,material='Si', hkl=[6,6,0], mask_d=60.0, bend_r=1.0, energy_resolution = 0.5, diced=False, thickness=500.0, database_dir=installation_dir):
        self.material     = material               # analyzer material
        self.hkl          = np.array(hkl)          # [hkl] indices of reflection used (shape (3,) numpy array)
        self.mask_d       = mask_d                 # analyzer mask diameter in [mm]
        self.bend_r       = bend_r                 # bending radius of the crystal [mm]
        self.diced        = diced                  # boolean (True or False) if a diced crystal is used or not (defalt is False)
        self.thickness    = thickness              # thickness of the analyzer crystal
        self.energy_resolution = energy_resolution # energy resolution [eV]
        self.energy_of_refl_calculation = None     # energy(dspace(hkl,material)) !!! check this again !!! may be obsolete or misleading
        self.database_dir     = database_dir       # path to a folder, where once calculated reflectivities are stored to spead up second time usage
        # output
        self.solid_angle      = []                 # solid angle of the analyzers
        self.efficiency       = []                 # factor, to be calculated
        self.reflectivity     = []                 # analyzer reflectivity (to be calculated)
        self.deviation_meV    = []                 # x-axis for reflectivity [meV]
        self.deviation_arcsec = []                 # x-axis for reflectivity [arc seconds]

    def set_material(self,material):
        if isinstance(material,str):
            self.material = material.upper()
        else:
            print( 'material must be passed as a string. Default is \'Si\'!')
    def get_material(self):
        if self.material:
            return self.material
        else:
            print( 'No material set, please set the material first!')
    def set_hkl(self,hkl):
        if np.shape(np.array(hkl)) == (3,):
            self.hkl = np.array(hkl)

    def get_hkl(self):
        if np.any(self.hkl):
            return self.hkl
        else:
            print( 'No hkl set, please set a hkl first!')
    def set_mask_d(self,mask_d):
        if isinstance(mask_d,int) or isinstance(mask_d,float):
            self.mask_d = np.array(mask_d)
        else:
            print( 'mask_d (analyzer mask diameter in mm) must be passed as either integer or float!')
    def get_mask_d(self):
        if np.any(self.mask_d):
            return self.mask_d
        else:
            print( 'No mask_d has been set')
    def set_bend_r(self,bend_r):
        if isinstance(bend_r,int) or isinstance(bend_r,float):
            self.bend_r = np.array(bend_r)
        else:
            print( 'bend_r (analyzer bending radius in m) must be passed as integer or float!')
    def get_bend_r(self):
        if np.any(self.bend_r):
            return self.bend_r
        else:
            print( 'No bend_r has been set')
    def set_diced(self,diced):
        if diced == True or diced == False:
            self.diced = diced
        else:
            print( 'diced must be either True or False')
    def get_diced(self):
        if self.diced:
            return self.diced
        else:
            print( 'diced has not been set')
    def set_thickness(self,thickness):
        self.thickness = thickness

    def get_thickness(self):
        if np.any(self.thickness):
            return self.thickness
        else:
            print( 'No thickness set, please set a thickness first!')
            return
    def get_energy_resolution(self):
        return self.energy_resolution

    def get_energy_resolution_keV(self):
        return self.energy_resolution*1.0e-3

    def get_solid_angle(self):
        det_area = 2.0*np.pi*(self.mask_d/2)**2.0
        sample_det_distance = self.bend_r*1.0e3/2
        return det_area/(4.0*np.pi*sample_det_distance**2.0)

    def get_reflectivity(self,energy,dev = np.arange(-50.0,150.0,1.0), alpha = 0.0):
        """
        Calculates the reflectivity curve for a given analyzer crystal. Checks 
        in the directory self.database_dir, if desired reflectivity curve has
        been calculated before.
        IN:
        energy = energy at which the reflectivity is to be calculated in [keV]
        dev    = deviation parameter for which the curve is to be calculated
        alpha  = deviation angle from exact Bragg angle [deg]
        """
        hkl      = self.get_hkl()
        material = self.get_material()
        bend_r   = self.get_bend_r()
        # print energy, hkl, material, bend_r, type(dev), type(alpha)
        # try opening reflectivity from file
        filename = os.path.join(
            self.database_dir,
            material + '_'  + 'hkl' + str(hkl) + '_'+ str(energy) + 'keV' + '.dat',
        )
        try:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
            database = shelve.open(filename)
            try:
                self.reflectivity     = database['reflectivity']
                self.deviation_meV    = database['deviation_meV']
                self.deviation_arcsec = database['deviation_arcsec']
                self.energy_of_refl_calculation = database['energy_of_refl_calculation']
            finally:
                database.close()
        except (FileNotFoundError, OSError, KeyError):
            # if no file exists, calculate reflectivity from scratch
            #print ('>>>>>>>>>>>>>>>', energy, hkl, material, bend_r, dev, alpha)
            reflectivity, e_scale, dev, e0 = xrs_utilities.taupgen(energy,hkl,material, bend_r,dev,alpha)
            self.reflectivity     = reflectivity
            self.deviation_meV    = e_scale
            self.deviation_arcsec = dev
            self.energy_of_refl_calculation = e0
            # save reflectivity for next time use
            # filename = self.database_dir +  material + '_'  + 'hkl' + str(hkl) + '_'+ str(energy) + 'keV' + '.dat'
            filename =   material + '_'  + 'hkl' + str(hkl) + '_'+ str(energy) + 'keV' + '.dat'
            # database =   shelve.open(filename)
            # database['reflectivity']     = reflectivity
            # database['deviation_meV']    = e_scale
            # database['deviation_arcsec'] = dev
            # database['energy_of_refl_calculation'] = e0

    def plot_reflectivity(self,mode='energy'):
        """
        Generates and opens a plot of the calculated reflectivity curve.
        mode = keyword for which x-axis is to be used, can be 'energy' or 'angle'
        """
        cla()
        dev_energy   = self.deviation_meV
        dev_arcsec   = self.deviation_arcsec
        reflectivity = self.reflectivity
        hkl          = self.get_hkl()
        E0           = self.energy_of_refl_calculation
        if mode == 'energy':
            pylab.plot(dev_energy,reflectivity)
            pylab.xlabel(['deviation from Bragg angle [meV]'])
            pylab.ylabel(['reflectivity [arb. units]'])
            titlestring = 'Takagi-Taupin curve for the ' + str(hkl) + ' reflection at %.2f' %E0 + ' keV.'
            pylab.title(titlestring)
            pylab.show()
        elif mode == 'angle':
            pylab.plot(dev_arcsec,reflectivity)
            pylab.xlabel(['deviation from Bragg angle [arcsec]'])
            pylab.ylabel(['reflectivity [arb. units]'])
            titlestring = 'Takagi-Taupin curve for the ' + str(hkl) + ' reflection at %.2f' %E0 + ' keV.'
            pylab.title(titlestring)
            pylab.show()
        else:
            print( 'mode unknown, please select either \'energy\' or \'angle\'.')
            return

    def get_efficiency(self,energy=None):
        """
        Calculates the efficiency of the analyzer crystal based on the calculated reflectivity curve.
        The efficiency is calculated by averaging over the energy resolution set upon class initialization.
        energy = energy (in [keV]) for wich the efficiency is to be calculated
        """
        if not energy:
            energy = self.energy_of_refl_calculation
        # print type(energy)
        energy_resolution = self.energy_resolution * 1.0e3 # resolution in meV
        if not np.any(self.reflectivity):
            self.get_reflectivity(energy)
        dev_energy      = self.deviation_meV
        reflectivity    = self.reflectivity
        # average over the FWHM of the reflectivity curve for an estimate of the efficiency
        fwhm, x0 = xrs_utilities.fwhm(dev_energy,reflectivity)
        inds            = np.where(np.logical_and(dev_energy>=x0-fwhm/2.0,dev_energy<=x0+fwhm/2.0))[0]
        self.efficiency = np.mean(reflectivity[inds])
        self.energy_resolution_meV = energy_resolution
        return self.efficiency

class sample:
    """
    Class to describe a sample.
    """
    def __init__(self,chem_formulas,concentrations,densities,angle_tth,sample_thickness,angle_in=None,angle_out=None,shape='sphere',molar_masses=None):
        self.chem_formulas  = chem_formulas     # list of strings of chemical sum formulas
        if np.isscalar(concentrations):
            concentrations = [concentrations]
        if np.isscalar(densities):
            densities = [densities]
        self.concentrations = concentrations    # list of concentrations, should contain values between 0.0 and 1.0
        self.densities      = densities         # list of densities of the constituents [g/cm^3]
        self.molar_masses   = molar_masses      # list of molar masses of all constituents
        self.shape          = shape             # keyword, can be 'slab' or 'sphere'
        self.tth            = angle_tth         # scattering angle [deg]
        self.alpha          = angle_in          # incident beam angle in [deg] relative to sample surface normal
        self.beta           = angle_out         # beam exit angle in [deg] relatice to sample surface normal (negative for transmission geometry)
        self.thickness      = sample_thickness  # sample thickness/diameter in [cm]
        self.energy1        = []
        self.energy2        = []

    def get_formulas(self):
        return self.chem_formulas
    def get_concentrations(self):
        return self.concentrations
    def get_densities(self):
        return self.densities
    def get_average_densities(self):
        return np.sum(self.densities)/len(self.densities)
    def get_shape(self):
        return self.shape
    def get_tth(self):
        return self.tth
    def get_thickness(self):
        return self.thickness
    def get_molar_masses(self):
        return self.molar_masses
    def get_energy1(self):
        return self.energy1
    def get_energy2(self):
        return self.energy2
    def get_alpha(self):
        if self.alpha:
            return self.alpha
        else:
            print( 'alpha has not been set!')
    def get_beta(self):
        if self.beta:
            return self.beta
        else:
            print( 'beta has not been set!')

    def get_murho(self,energy1,energy2=None):
        """
        Calculates the total photoelectric absorption coefficient of the sample
        for the two energies given. Returns only one array, if only one energy axis 
        is defined.
        energy1 = numpy array of energies in [keV] 
        energy2 = numpy array of energies in [keV] (defalt is None, i.e. only one mu is returned) 
        """
        self.energy1 = energy1
        self.energy2 = energy2

        energy         = np.asarray(energy1)
        formulas       = self.get_formulas()
        concentrations = self.get_concentrations()
        E0             = energy2
        rho_formu      = self.get_densities()

        if energy2 is not None:
            return xrs_utilities.mpr_compds(energy,formulas,concentrations,E0,rho_formu) # returns mu_in and mu_out
        else:
            E0 = np.atleast_1d(energy)[-1]
            mu_tot_in, mu_tot_out = xrs_utilities.mpr_compds(energy,formulas,concentrations,E0,rho_formu)
            return mu_tot_in # returns only mu_in

    def get_absorption_correction(self,energy1,energy2,thickness=None):
        """
        Calculates the absorption correction factor for the sample
        to be multiplied with experimental data to correct for absorption effects.
        energy1 = numpy array of energies in [keV] for which the factor is to be calculated
        energy2 = numpy array of energies in [keV] for which the factor is to be calculated
        """
        alpha = self.alpha
        beta  = self.beta
        tth   = self.tth
        if not thickness:
            thickness = self.thickness

        mu_tot_in, mu_tot_out = self.get_murho(energy1,energy2)

        if isinstance(tth,list):
            ac = (mu_tot_in + mu_tot_out)/(1.0 - np.exp(-mu_tot_in*thickness -mu_tot_out*thickness))

        else:
            if self.shape == 'slab' and alpha and beta:
                if tth:
                    if self.beta<0: # transmission geometry
                        test_tth = alpha-beta
                    else: # reflection geometry
                        test_tth = 180.0 - (alpha+beta)
                    if tth != test_tth:
                        print( 'the alpha and beta values set are not congruent to the tth value set!')
                ac = xrs_utilities.abscorr2(mu_tot_in,mu_tot_out,alpha,beta,thickness)
            elif self.shape == 'sphere':
                ac = (mu_tot_in + mu_tot_out)/(1.0 - np.exp(-mu_tot_in*thickness -mu_tot_out*thickness)) 
                #1.0/np.exp(-thickness*mu_tot_in -thickness*mu_tot_out) # spherical sample just add up in and outgoing absorption
            else:
                raise ValueError('please provide either shape=\'sphere\' (default) or \'slab\' and alpha and beta')
        return ac

    def plot_inv_absorption(self,energy1,energy2,range_of_thickness = np.arange(0.0,0.5,0.01)):
        """
        Generates a figure which plots 1/Abscorr for the sample as a function of different thicknesses.
        This is usefull for finding optimum sample thicknesses for an experiment.
        energy1 = energy in [keV] at the desired edge
        energy2 = energy in [keV] at the elastic
        range_of_thickness = numpy array of sample thicknesses in [cm]

        !!! right now all samples are treates as if spherical !!!
        """
        ac_range = np.zeros_like(range_of_thickness)
        for ii in range(len(ac_range)):
            ac_range[ii] = self.get_absorption_correction(energy1,energy2,range_of_thickness[ii])
        cla()
        pylab.plot(range_of_thickness,1/ac_range)
        pylab.xlabel('sample thickness [cm]')
        pylab.ylabel('1/(absorption correction factor) [arb. units]')
        pylab.show()

class thomson:
    """
    Class to take care of the Thomson scattering cross section.
    """
    def __init__(self,omega_1,omega_2,tth,scattering_plane='vertical',polarization=0.99):
        self.omega_1 = omega_1                    # numpy array of primary energy in [keV]
        self.omega_2 = omega_2                    # analyzer energy in [keV]
        self.tth     = tth                        # scattering angle in [deg]
        self.scattering_plane  = scattering_plane # keyword to indicate scattering plane relative to lab frame ('vertical' or 'horizontal')
        self.polarization      = polarization     # degree of polarization (close to 1.0 for undulator radiation)
        self.r0 = constants.physical_constants['classical electron radius'][0]*1e2 # classical electron radius in [cm]

    def get_thomson_factor(self):
        """
        Calculates the Thomson scattering factor.
        """
        # mutiple tth values in a list
        if isinstance(self.tth,list):
            thomson = np.zeros((self.omega_1.shape[0], len(self.tth)))
            if self.scattering_plane == 'vertical':
                for ii in range(len(self.tth)):
                    thomson[:,ii] = self.polarization * 1.0 * self.omega_2/self.omega_1 * self.r0**2.0
            elif self.scattering_plane == 'horizontal':
                for ii in range(len(self.tth)):
                    thomson[:,ii] = self.polarization * np.cos(np.radians(self.tth[ii]))**2.0  * self.omega_2/self.omega_1 * self.r0**2.0
            else:
                print( 'the scattering plane can only be \'vertical\' or \'horizontal\'.')
                return
        # just one tth value
        else:
            tth = [self.tth]
            thomson = np.zeros((self.omega_1.shape[0], len(tth)))
            if self.scattering_plane == 'vertical':
                thomson[:,0] = self.omega_2/self.omega_1 * self.r0**2.0 * self.polarization
            elif self.scattering_plane == 'horizontal':
                thomson[:,0] = self.omega_2/self.omega_1 * self.r0**2.0 * self.polarization * np.cos(np.radians(tth[0]))**2.0
            else:
                print( 'the scattering plane can only be \'vertical\' or \'horizontal\'.')
                return
        return thomson

class beam:
    """
    Class to describe incident beam related things.
    """
    def __init__(self,i0_intensity,beam_height,beam_width,divergence=None):
        self.i0_intensity = i0_intensity # number of incident photons [1/sec]
        self.beam_height  = beam_height  # in micron
        self.beam_width   = beam_width   # in micron
        self.divergence   = divergence   # in milli-rad
    def get_i0_intensity(self):
        return self.i0_intensity
    def get_beam_height(self):
        return self.beam_height
    def get_beam_height_cm(self):
        return self.beam_height * 1.0e-4
    def get_beam_width(self):
        return self.beam_height
    def get_beam_width_cm(self):
        return self.beam_width * 1.0e-4
    def get_divergence(self):
        return self.divergence
    def get_beam_cross_section_area(self):
        """
        Calculates the beam cross section area.
        """
        return self.beam_height * self.beam_width # in [microns^2]

class compton_profiles:
    """
    Class to hold construct HF Compton profiles for an object of the sample class.
    """
    def __init__(self,sample_obj,eloss_range=np.arange(0.0,1000.0,0.1),E0=9.7):
        self.chem_formulas  = sample_obj.get_formulas()
        self.concentrations = sample_obj.get_concentrations()
        self.densities      = sample_obj.get_densities()
        self.mean_density   = 0.0
        if len(self.densities)>1:
            for ii in range(len(self.densities)):
                self.mean_density += self.densities[ii]*self.concentrations[ii]
        else:
            self.mean_density = self.densities
        self.E0             = E0          # elastic line energy in [keV]
        self.eloss_range    = eloss_range # desired energy loss range in [eV]
        self.sample_shape   = sample_obj.get_shape()
        self.tth            = sample_obj.get_tth()
        self.thickness      = sample_obj.get_thickness()
        self.ac_factor      = sample_obj.get_absorption_correction(E0+eloss_range,E0)

        # output
        if isinstance(self.tth,list):
            self.J = np.zeros((len(self.eloss_range),len(self.tth)))
            self.C = np.zeros((len(self.eloss_range),len(self.tth)))
            self.V = np.zeros((len(self.eloss_range),len(self.tth)))
            self.q = np.zeros((len(self.eloss_range),len(self.tth)))
        else:
            self.J = np.array([])
            self.C = np.array([])
            self.V = np.array([])
            self.q = np.array([])

    def get_E0(self):
        return self.E0
    def get_energy_in_keV(self):
        return self.eloss_range/1e3 + self.E0
    def get_tth(self):
        return self.tth

    def calc_pure_HF_profiles(self):
        if isinstance(self.tth,list):
            for tth,ii in zip(self.tth,list(range(len(self.tth)))):
                eloss,J,C,V,q = xrs_utilities.makeprofile_compds(self.chem_formulas,concentrations=self.concentrations,E0=self.E0,tth=tth)
                self.J[:,ii] = np.interp(self.eloss_range,eloss,J)
                self.C[:,ii] = np.interp(self.eloss_range,eloss,C)
                self.V[:,ii] = np.interp(self.eloss_range,eloss,V)
                self.q[:,ii] = np.interp(self.eloss_range,eloss,q)
        else:
            eloss,J,C,V,q = xrs_utilities.makeprofile_compds(self.chem_formulas,concentrations=self.concentrations,E0=self.E0,tth=self.tth)
            self.J = np.interp(self.eloss_range,eloss,J)
            self.C = np.interp(self.eloss_range,eloss,C)
            self.V = np.interp(self.eloss_range,eloss,V)
            self.q = np.interp(self.eloss_range,eloss,q)

    def calc_HF_profiles(self):
        if isinstance(self.tth,list):
            if not np.any(self.J) and not np.any(self.C) and not np.any(self.V) and not np.any(self.q):
                self.calc_pure_HF_profiles()
            for tth,ii in zip(self.tth,list(range(len(self.tth)))):
                self.J[:,ii] = self.J[:,ii]/self.ac_factor*self.mean_density
                self.C[:,ii] = self.C[:,ii]/self.ac_factor*self.mean_density
                self.V[:,ii] = self.V[:,ii]/self.ac_factor*self.mean_density
        else:
            if not np.any(self.J) and not np.any(self.C) and not np.any(self.V) and not np.any(self.q):
                self.calc_pure_HF_profiles() # calculate uncorrected profiles
            self.J = self.J/self.ac_factor*self.mean_density
            self.C = self.C/self.ac_factor*self.mean_density
            self.V = self.V/self.ac_factor*self.mean_density

    def get_HF_profiles(self):
        if not np.any(self.J) and not np.any(self.C) and not np.any(self.V) and not np.any(self.q):
            self.calc_HF_profiles() # calculate uncorrected profiles
        return self.eloss_range, self.J, self.C, self.V, self.q

    def plot_HF_profile(self):
        if not np.any(self.J) and not np.any(self.C) and not np.any(self.V) and not np.any(self.q):
            self.calc_HF_profiles()
        cla()
        pylab.plot(self.eloss_range,self.J)
        pylab.plot(self.eloss_range,self.C)
        pylab.plot(self.eloss_range,self.V)
        pylab.xlabel('energy loss [eV]')
        pylab.ylabel('intensity [1/eV]')
        pylab.title('sample absorption corrected HF Compton profile')
        pylab.show()

class absolute_cross_section:
    """
    Class to calculate an expected cross section in absolute counts using objects of the 'beam', 'sample',
    'analyzer', 'detector', 'thomson', and 'compton_profile' classes.
    """
    def __init__(self,beam_obj, sample_obj, analyzer_obj, detector_obj, thomson_obj, compton_profile_obj):
        self.eloss,self.J,self.C,self.V,self.q = compton_profile_obj.get_HF_profiles()
        self.thomson = thomson_obj.get_thomson_factor()
        self.I0      = beam_obj.get_i0_intensity()
        self.beam_h  = beam_obj.get_beam_width_cm()
        self.beam_v  = beam_obj.get_beam_height_cm()
        self.sample_tth            = sample_obj.get_tth()
        self.sample_densities      = sample_obj.get_densities()
        self.sample_molar_masses   = sample_obj.get_molar_masses()
        self.sample_formulas       = sample_obj.get_formulas()
        self.sample_thickness      = sample_obj.get_thickness()
        self.sample_concentrations = sample_obj.get_concentrations()
        self.sample_abs_in         = sample_obj.get_murho(compton_profile_obj.get_energy_in_keV())
        self.det_efficiency = detector_obj.get_efficiency(compton_profile_obj.get_E0())
        self.ana_efficiency = analyzer_obj.get_efficiency(compton_profile_obj.get_E0())
        self.ana_solid_angle = analyzer_obj.get_solid_angle()
        self.ana_energy_resolution = analyzer_obj.get_energy_resolution()
        self.energy_in_keV = compton_profile_obj.get_energy_in_keV()
        # output
        self.absolute_counts = []

    def calc_num_scatterers(self):
        """
        Calculates number of scatterers/atoms using beam size, sample thickness, sample densites,
        sample molar masses (so far does not differentiate between target atoms and random sample atoms)
        """
        sample_volume = self.beam_h * self.beam_v * self.sample_thickness
        sample_weight = np.zeros(len(self.sample_densities))
        sample_amount = np.zeros(len(self.sample_densities))
        for ii in range(len(self.sample_densities)):
            sample_weight[ii] += self.sample_densities[ii] * sample_volume * self.sample_concentrations[ii]
            sample_amount[ii] += sample_weight[ii] * self.sample_molar_masses[ii]
        return np.sum(sample_amount)*constants.physical_constants['Avogadro constant'][0]

    def calc_abs_cross_section(self):
        num_of_scatterers = self.calc_num_scatterers()
        if isinstance(self.sample_tth, list):
            self.absolute_counts = np.zeros_like(self.J)
            for ii in range(len(self.sample_tth)):
                self.absolute_counts[:,ii] = self.I0 * self.thomson[:,ii] * self.J[:,ii] * self.ana_solid_angle * self.ana_energy_resolution *  num_of_scatterers * self.sample_thickness * self.sample_abs_in * self.ana_efficiency * self.det_efficiency
        else:
            self.absolute_counts = self.I0 * self.thomson[:,0] * self.J * self.ana_solid_angle * self.ana_energy_resolution *  num_of_scatterers * self.sample_thickness * self.sample_abs_in * self.ana_efficiency * self.det_efficiency

    def plot_abs_cross_section(self):
        if not np.any(self.absolute_counts):
            self.calc_abs_cross_section()
        cla()
        pylab.plot(self.eloss,self.absolute_counts)
        pylab.xlabel('energy loss')
        pylab.ylabel('absolute counts [1/sec]')
        pylab.show()

    def save_txt(self, file_name, header=''):
        if not np.any(self.absolute_counts):
            self.calc_abs_cross_section()
        absolute_counts = np.asarray(self.absolute_counts)
        if absolute_counts.ndim == 1:
            absolute_counts = absolute_counts[:,np.newaxis]
        data = np.zeros((len(self.eloss),absolute_counts.shape[1]+1))
        data[:,0] = self.eloss
        data[:,1::] = absolute_counts
        output_dir = os.path.dirname(os.path.abspath(file_name))
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        np.savetxt(file_name, data, header=header)

    def save_hdf5( self, fname, group_name="sample1" ):
        """ **save_hdf5**
        Save the results in an HDF5 file.
        Note:
            HDF5 files are strange for overwriting files.

        Args:
            fname (str): Path and filename for the HDF5 file.
        """
        if isinstance(fname, h5py.Group):
            f=fname
        else:
            f = h5py.File(fname, "w")

        h5_group = f.require_group(group_name)
        # absolute counts for plotting
        if( len(self.absolute_counts.shape)==1):
             self.absolute_counts = np.array([self.absolute_counts] ).T
        data = np.zeros((len(self.eloss),self.absolute_counts.shape[1]+1))
        data[:,0] = self.eloss
        data[:,1::] = self.absolute_counts
        h5_group["abs_counts"] = data
        f.flush()
        f.close()      

def input_file_parser(filename):
    """
    Parses an input file, which has a structure like the example input file ('prediction.inp') provided in the
    examples/ folder. (Python lists and numpy arrays have to be profived without white spaces in their definitions,
    e.g. 'hkl = [6,6,0]' instead of 'hkl = [6, 6, 0]')
    """
    try:
        with open(filename, 'r') as handle:
            lines = handle.readlines()
    except OSError as exc:
        raise FileNotFoundError('No input file ' + filename + ' found.') from exc
    return _parse_input_lines(lines)


def _parse_input_lines(lines):
    input_parameters = {} # dictionary of input parameters
    section_names = ['detector','analyzer','sample','thomson','beam','compton_profiles']
    for name in section_names:
        input_parameters[name] = {} # cempty list for all sections, fill them up from file, add defaults for missing values
    # parse all given parameters
    lineindex = 0
    while lineindex < len(lines):
        line = lines[lineindex].strip()
        if line.startswith('####'):
            parts = line.split()
            if len(parts) < 2:
                lineindex += 1
                continue
            thekey = parts[1]
            if thekey not in input_parameters:
                input_parameters[thekey] = {}
            lineindex += 1
            while lineindex < len(lines) and not lines[lineindex].lstrip().startswith('####'):
                item = lines[lineindex].strip()
                if item and not item.startswith('#'):
                    key, value = _parse_input_assignment(item)
                    input_parameters[thekey][key] = value
                lineindex += 1
            continue
        lineindex += 1
    return input_parameters, section_names


def _parse_input_assignment(line):
    if '=' not in line:
        raise ValueError("Expected input assignment in the form 'name = value': %s" % line)
    key, value = [part.strip() for part in line.split('=', 1)]
    try:
        return key, ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return key, _parse_numpy_expression(value)


def _parse_numpy_expression(value):
    tree = ast.parse(value, mode='eval')
    allowed_calls = {
        ('np', 'arange'): np.arange,
        ('numpy', 'arange'): np.arange,
        ('np', 'linspace'): np.linspace,
        ('numpy', 'linspace'): np.linspace,
        ('np', 'logspace'): np.logspace,
        ('numpy', 'logspace'): np.logspace,
        ('np', 'array'): np.array,
        ('numpy', 'array'): np.array,
    }

    def convert(node):
        if isinstance(node, ast.Expression):
            return convert(node.body)
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = convert(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.List):
            return [convert(item) for item in node.elts]
        if isinstance(node, ast.Tuple):
            return tuple(convert(item) for item in node.elts)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name):
            call_key = (node.func.value.id, node.func.attr)
            if call_key not in allowed_calls:
                raise ValueError("Unsupported expression in prediction input: %s" % value)
            args = [convert(arg) for arg in node.args]
            kwargs = {keyword.arg: convert(keyword.value) for keyword in node.keywords}
            return allowed_calls[call_key](*args, **kwargs)
        raise ValueError("Unsupported expression in prediction input: %s" % value)

    return convert(tree)

def get_all_input(filename = 'prediction.inp'):
    """
    Adds default values if input is missing in the input-file and a default value exists for the missing one.
    """
    # parse the input file
    input_parameters, section_names = input_file_parser(filename)
    # create something for all possible inputs
    all_input = {}    
    for name in section_names:
        all_input[name] = {}
    # detector
    all_input['detector']['energy']    = 9.7  # analyzer energy in keV
    all_input['detector']['thickness'] = 500.0 # detector thickness
    all_input['detector']['material']  = 'Si'  # detector material
    all_input['detector']['pixel_size']= [256,768] #detector pixel size
    # analyzer
    all_input['analyzer']['material'] = 'Si' # analyzer crystal material 
    all_input['analyzer']['hkl']      = [6,6,0] # analyzer crystal reflection
    all_input['analyzer']['mask_d']   = 60.0    # analyzer mask diameter in mm
    all_input['analyzer']['bend_r']   = 1.0     # analyzer bending radius in m
    all_input['analyzer']['energy_resolution']   = 0.5 # resolution in eV
    all_input['analyzer']['diced']        = False    # keyword, if bent or diced analyzer is used
    all_input['analyzer']['thickness']    = 500.0    # analyzer bending radius in m
    all_input['analyzer']['database_dir'] = installation_dir  # directory to tabulated chi tables (for calculation of reflectivities)
    # sample
    all_input['sample']['chem_formulas']    = []
    all_input['sample']['concentrations']   = []
    all_input['sample']['densities']        = []
    all_input['sample']['angle_tth']        = []
    all_input['sample']['sample_thickness'] = []
    all_input['sample']['angle_in'] = None    # up to now, only spherical samples possible !!!
    all_input['sample']['angle_out'] = None   # same here
    all_input['sample']['shape'] = 'sphere'    # sample shape, right now, only this works, should be 'slab' or 'sphere' in the future
    all_input['sample']['molar_masses'] = None # this is needed for the estimation of number of scatterers. should be mandatory, acutally
    # thomson
    all_input['thomson']['omega_1'] = []
    all_input['thomson']['omega_2'] = [] 
    all_input['thomson']['tth']     = []
    all_input['thomson']['scattering_plane'] = 'vertical' # 'vertical' or 'horizontal', for polarization purposes
    all_input['thomson']['polarization'] = 0.99           # polarization factor
    # beam
    all_input['beam']['i0_intensity'] = []
    all_input['beam']['beam_height']  = []
    all_input['beam']['beam_width']   = []
    all_input['beam']['divergence'] = None # this is just a dummy parameter
    # compton_profiles
    all_input['compton_profiles']['eloss_range'] = np.arange(0.0,1000.0,0.1) # energy range in eV
    all_input['compton_profiles']['E0'] = 9.7 # analyzer energy in keV
    # if present, replace all_input variable by values provided in the input file:
    for key in input_parameters:
        for key2 in input_parameters[key]:
            all_input[key][key2] =  input_parameters[key][key2]
    return all_input


def _default_all_input():
    all_input = {}
    section_names = ['detector','analyzer','sample','thomson','beam','compton_profiles']
    for name in section_names:
        all_input[name] = {}
    # detector
    all_input['detector']['energy']    = 9.7  # analyzer energy in keV
    all_input['detector']['thickness'] = 500.0 # detector thickness
    all_input['detector']['material']  = 'Si'  # detector material
    all_input['detector']['pixel_size']= [256,768] #detector pixel size
    # analyzer
    all_input['analyzer']['material'] = 'Si' # analyzer crystal material 
    all_input['analyzer']['hkl']      = [6,6,0] # analyzer crystal reflection
    all_input['analyzer']['mask_d']   = 60.0    # analyzer mask diameter in mm
    all_input['analyzer']['bend_r']   = 1.0     # analyzer bending radius in m
    all_input['analyzer']['energy_resolution']   = 0.5 # resolution in eV
    all_input['analyzer']['diced']        = False    # keyword, if bent or diced analyzer is used
    all_input['analyzer']['thickness']    = 500.0    # analyzer bending radius in m
    all_input['analyzer']['database_dir'] = installation_dir  # directory to tabulated chi tables (for calculation of reflectivities)
    # sample
    all_input['sample']['chem_formulas']    = []
    all_input['sample']['concentrations']   = []
    all_input['sample']['densities']        = []
    all_input['sample']['angle_tth']        = []
    all_input['sample']['sample_thickness'] = []
    all_input['sample']['angle_in'] = None    # up to now, only spherical samples possible !!!
    all_input['sample']['angle_out'] = None   # same here
    all_input['sample']['shape'] = 'sphere'    # sample shape, right now, only this works, should be 'slab' or 'sphere' in the future
    all_input['sample']['molar_masses'] = None # this is needed for the estimation of number of scatterers. should be mandatory, acutally
    # thomson
    all_input['thomson']['omega_1'] = []
    all_input['thomson']['omega_2'] = [] 
    all_input['thomson']['tth']     = []
    all_input['thomson']['scattering_plane'] = 'vertical' # 'vertical' or 'horizontal', for polarization purposes
    all_input['thomson']['polarization'] = 0.99           # polarization factor
    # beam
    all_input['beam']['i0_intensity'] = []
    all_input['beam']['beam_height']  = []
    all_input['beam']['beam_width']   = []
    all_input['beam']['divergence'] = None # this is just a dummy parameter
    # compton_profiles
    all_input['compton_profiles']['eloss_range'] = np.arange(0.0,1000.0,0.1) # energy range in eV
    all_input['compton_profiles']['E0'] = 9.7 # analyzer energy in keV
    return all_input


def get_all_input_from_text(text):
    input_parameters, section_names = _parse_input_lines(text.splitlines())
    all_input = _default_all_input()
    for key in input_parameters:
        for key2 in input_parameters[key]:
            all_input[key][key2] = input_parameters[key][key2]
    return all_input


WEB_DEFAULT_PARAMETERS = {
    'chem_formulas': 'Si',
    'concentrations': '1.0',
    'densities': '2.33',
    'molar_masses': '28.0855',
    'angle_tth': 60.0,
    'sample_thickness': 0.01,
    'E0': 9.7,
    'eloss_min': 0.0,
    'eloss_max': 1000.0,
    'points': 220,
    'i0_intensity': 1.0e12,
    'beam_height': 50.0,
    'beam_width': 50.0,
    'detector_energy': 9.7,
    'detector_thickness': 500.0,
    'detector_material': 'Si',
    'analyzer_material': 'Si',
    'hkl': '6,6,0',
    'mask_d': 60.0,
    'bend_r': 1.0,
    'energy_resolution': 0.5,
    'analyzer_efficiency': 0.5,
    'use_reflectivity': False,
    'scattering_plane': 'vertical',
    'polarization': 0.99,
    'output_txt': 'prediction_output.txt',
}


def _parse_literal_value(value):
    if not isinstance(value, str):
        return value
    stripped = value.strip()
    if not stripped:
        return value
    if stripped[0] not in '[({"\'':
        return value
    try:
        return ast.literal_eval(stripped)
    except (SyntaxError, ValueError):
        return value


def _format_web_value(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return repr(list(value))
    return value


def _number_list(value, default=None, item_type=float):
    if value is None or value == '':
        value = default
    value = _parse_literal_value(value)
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return [item_type(item) for item in value]
    if isinstance(value, str):
        items = [item.strip() for item in value.replace(';', ',').split(',')]
        return [item_type(item) for item in items if item]
    return [item_type(value)]


def _string_list(value, default=None):
    if value is None or value == '':
        value = default
    value = _parse_literal_value(value)
    if value is None:
        return []
    if isinstance(value, (list, tuple, np.ndarray)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [item.strip() for item in str(value).replace(';', ',').split(',') if item.strip()]


def _scalar(value, default, value_type=float):
    if value is None or value == '':
        return value_type(default)
    return value_type(value)


def _normalize_component_lengths(formulas, concentrations, densities, molar_masses):
    component_count = len(formulas)
    if component_count == 0:
        raise ValueError('At least one chemical formula is required.')

    def expand(values, name):
        if not values:
            raise ValueError(name + ' is required.')
        if len(values) == 1 and component_count > 1:
            return values * component_count
        if len(values) != component_count:
            raise ValueError(name + ' must have one value or match the number of formulas.')
        return values

    return (
        expand(concentrations, 'concentrations'),
        expand(densities, 'densities'),
        expand(molar_masses, 'molar_masses'),
    )


def _first_value(value, default):
    if value is None:
        return default
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return default
        return value.reshape(-1)[0].item()
    if isinstance(value, (list, tuple)):
        if not value:
            return default
        return value[0]
    return value


def _has_value(value):
    if value is None:
        return False
    if isinstance(value, np.ndarray):
        return value.size > 0
    if isinstance(value, (list, tuple, str)):
        return len(value) > 0
    return True


def web_parameters_from_all_input(all_input):
    params = dict(WEB_DEFAULT_PARAMETERS)

    detector_params = all_input.get('detector', {})
    analyzer_params = all_input.get('analyzer', {})
    sample_params = all_input.get('sample', {})
    thomson_params = all_input.get('thomson', {})
    beam_params = all_input.get('beam', {})
    compton_params = all_input.get('compton_profiles', {})

    formulas = sample_params.get('chem_formulas')
    if _has_value(formulas):
        params['chem_formulas'] = _format_web_value(formulas)
    concentrations = sample_params.get('concentrations')
    if _has_value(concentrations):
        params['concentrations'] = _format_web_value(concentrations)
    densities = sample_params.get('densities')
    if _has_value(densities):
        params['densities'] = _format_web_value(densities)
    molar_masses = sample_params.get('molar_masses')
    if _has_value(molar_masses):
        params['molar_masses'] = _format_web_value(molar_masses)

    params['angle_tth'] = _format_web_value(
        sample_params.get('angle_tth', WEB_DEFAULT_PARAMETERS['angle_tth'])
    )
    params['sample_thickness'] = _first_value(
        sample_params.get('sample_thickness'),
        WEB_DEFAULT_PARAMETERS['sample_thickness'],
    )

    eloss_range = np.asarray(compton_params.get('eloss_range', []), dtype=float)
    if eloss_range.size:
        params['eloss_min'] = float(np.min(eloss_range))
        params['eloss_max'] = float(np.max(eloss_range))
        params['points'] = int(min(max(eloss_range.size, 20), 2000))
    params['E0'] = _first_value(compton_params.get('E0'), WEB_DEFAULT_PARAMETERS['E0'])

    params['i0_intensity'] = _first_value(
        beam_params.get('i0_intensity'),
        WEB_DEFAULT_PARAMETERS['i0_intensity'],
    )
    params['beam_height'] = _first_value(
        beam_params.get('beam_height'),
        WEB_DEFAULT_PARAMETERS['beam_height'],
    )
    params['beam_width'] = _first_value(
        beam_params.get('beam_width'),
        WEB_DEFAULT_PARAMETERS['beam_width'],
    )

    params['detector_energy'] = _first_value(
        detector_params.get('energy'),
        WEB_DEFAULT_PARAMETERS['detector_energy'],
    )
    params['detector_thickness'] = _first_value(
        detector_params.get('thickness'),
        WEB_DEFAULT_PARAMETERS['detector_thickness'],
    )
    params['detector_material'] = detector_params.get(
        'material',
        WEB_DEFAULT_PARAMETERS['detector_material'],
    )

    params['analyzer_material'] = analyzer_params.get(
        'material',
        WEB_DEFAULT_PARAMETERS['analyzer_material'],
    )
    params['hkl'] = _format_web_value(
        analyzer_params.get('hkl', WEB_DEFAULT_PARAMETERS['hkl'])
    )
    params['mask_d'] = _first_value(analyzer_params.get('mask_d'), WEB_DEFAULT_PARAMETERS['mask_d'])
    params['bend_r'] = _first_value(analyzer_params.get('bend_r'), WEB_DEFAULT_PARAMETERS['bend_r'])
    params['energy_resolution'] = _first_value(
        analyzer_params.get('energy_resolution'),
        WEB_DEFAULT_PARAMETERS['energy_resolution'],
    )

    params['scattering_plane'] = thomson_params.get(
        'scattering_plane',
        WEB_DEFAULT_PARAMETERS['scattering_plane'],
    )
    params['polarization'] = _first_value(
        thomson_params.get('polarization'),
        WEB_DEFAULT_PARAMETERS['polarization'],
    )
    return params


def web_parameters_from_inp_text(text):
    return web_parameters_from_all_input(get_all_input_from_text(text))


def predict_from_parameters(parameters=None):
    """
    Calculate an XRS prediction from JSON-friendly parameters.

    This helper is used by the web UI and can also be called from notebooks.
    By default it uses a fast analyzer-efficiency value so interactive plots
    update quickly. Set use_reflectivity=True to calculate the analyzer
    efficiency with the existing Takagi-Taupin path.
    """
    params = dict(WEB_DEFAULT_PARAMETERS)
    if parameters:
        params.update(parameters)

    formulas = _string_list(params.get('chem_formulas'), WEB_DEFAULT_PARAMETERS['chem_formulas'])
    concentrations = _number_list(params.get('concentrations'), WEB_DEFAULT_PARAMETERS['concentrations'])
    densities = _number_list(params.get('densities'), WEB_DEFAULT_PARAMETERS['densities'])
    molar_masses = _number_list(params.get('molar_masses'), WEB_DEFAULT_PARAMETERS['molar_masses'])
    concentrations, densities, molar_masses = _normalize_component_lengths(
        formulas, concentrations, densities, molar_masses
    )

    eloss_min = _scalar(params.get('eloss_min'), WEB_DEFAULT_PARAMETERS['eloss_min'])
    eloss_max = _scalar(params.get('eloss_max'), WEB_DEFAULT_PARAMETERS['eloss_max'])
    points = max(20, min(2000, _scalar(params.get('points'), WEB_DEFAULT_PARAMETERS['points'], int)))
    if eloss_max <= eloss_min:
        raise ValueError('eloss_max must be greater than eloss_min.')

    E0 = _scalar(params.get('E0'), WEB_DEFAULT_PARAMETERS['E0'])
    tth_values = _number_list(params.get('angle_tth'), WEB_DEFAULT_PARAMETERS['angle_tth'])
    angle_tth = tth_values[0] if len(tth_values) == 1 else tth_values
    eloss_range = np.linspace(eloss_min, eloss_max, points)

    sample_obj = sample(
        formulas,
        concentrations,
        densities,
        angle_tth,
        _scalar(params.get('sample_thickness'), WEB_DEFAULT_PARAMETERS['sample_thickness']),
        shape='sphere',
        molar_masses=molar_masses,
    )
    beam_obj = beam(
        _scalar(params.get('i0_intensity'), WEB_DEFAULT_PARAMETERS['i0_intensity']),
        _scalar(params.get('beam_height'), WEB_DEFAULT_PARAMETERS['beam_height']),
        _scalar(params.get('beam_width'), WEB_DEFAULT_PARAMETERS['beam_width']),
    )
    analyzer_obj = analyzer(
        material=str(params.get('analyzer_material', WEB_DEFAULT_PARAMETERS['analyzer_material'])),
        hkl=_number_list(params.get('hkl'), WEB_DEFAULT_PARAMETERS['hkl'], int),
        mask_d=_scalar(params.get('mask_d'), WEB_DEFAULT_PARAMETERS['mask_d']),
        bend_r=_scalar(params.get('bend_r'), WEB_DEFAULT_PARAMETERS['bend_r']),
        energy_resolution=_scalar(params.get('energy_resolution'), WEB_DEFAULT_PARAMETERS['energy_resolution']),
        diced=False,
        thickness=500.0,
        database_dir=installation_dir,
    )
    detector_obj = detector(
        energy=_scalar(params.get('detector_energy'), WEB_DEFAULT_PARAMETERS['detector_energy']),
        thickness=_scalar(params.get('detector_thickness'), WEB_DEFAULT_PARAMETERS['detector_thickness']),
        material=str(params.get('detector_material', WEB_DEFAULT_PARAMETERS['detector_material'])),
        pixel_size=[256, 768],
    )
    compton_profile_obj = compton_profiles(sample_obj, eloss_range, E0)
    eloss, J, C, V, q = compton_profile_obj.get_HF_profiles()
    thomson_obj = thomson(
        compton_profile_obj.get_energy_in_keV(),
        compton_profile_obj.get_E0(),
        compton_profile_obj.get_tth(),
        scattering_plane=str(params.get('scattering_plane', WEB_DEFAULT_PARAMETERS['scattering_plane'])),
        polarization=_scalar(params.get('polarization'), WEB_DEFAULT_PARAMETERS['polarization']),
    )

    use_reflectivity = bool(params.get('use_reflectivity'))
    if use_reflectivity:
        analyzer_efficiency = float(np.asarray(analyzer_obj.get_efficiency(E0)).reshape(-1)[0])
        efficiency_mode = 'calculated'
    else:
        analyzer_efficiency = _scalar(
            params.get('analyzer_efficiency'),
            WEB_DEFAULT_PARAMETERS['analyzer_efficiency'],
        )
        efficiency_mode = 'manual'

    detector_efficiency = float(np.asarray(detector_obj.get_efficiency(E0)).reshape(-1)[0])
    sample_abs_in = sample_obj.get_murho(compton_profile_obj.get_energy_in_keV())
    thomson_factor = thomson_obj.get_thomson_factor()

    sample_volume = beam_obj.get_beam_width_cm() * beam_obj.get_beam_height_cm() * sample_obj.get_thickness()
    sample_amount = np.zeros(len(densities))
    for index in range(len(densities)):
        sample_weight = densities[index] * sample_volume * concentrations[index]
        sample_amount[index] = sample_weight * molar_masses[index]
    num_scatterers = float(np.sum(sample_amount) * constants.physical_constants['Avogadro constant'][0])

    J_array = np.asarray(J)
    if J_array.ndim == 1:
        counts = (
            beam_obj.get_i0_intensity()
            * thomson_factor[:, 0]
            * J_array
            * analyzer_obj.get_solid_angle()
            * analyzer_obj.get_energy_resolution()
            * num_scatterers
            * sample_obj.get_thickness()
            * sample_abs_in
            * analyzer_efficiency
            * detector_efficiency
        )
        series = [{'name': 'absolute counts', 'values': counts.tolist()}]
    else:
        series = []
        for index in range(J_array.shape[1]):
            counts = (
                beam_obj.get_i0_intensity()
                * thomson_factor[:, index]
                * J_array[:, index]
                * analyzer_obj.get_solid_angle()
                * analyzer_obj.get_energy_resolution()
                * num_scatterers
                * sample_obj.get_thickness()
                * sample_abs_in
                * analyzer_efficiency
                * detector_efficiency
            )
            series.append({'name': 'tth ' + str(tth_values[index]), 'values': counts.tolist()})

    return {
        'parameters': params,
        'x': np.asarray(eloss).tolist(),
        'series': series,
        'profiles': {
            'J': np.asarray(J).tolist(),
            'C': np.asarray(C).tolist(),
            'V': np.asarray(V).tolist(),
            'q': np.asarray(q).tolist(),
        },
        'meta': {
            'detector_efficiency': detector_efficiency,
            'analyzer_efficiency': analyzer_efficiency,
            'analyzer_efficiency_mode': efficiency_mode,
            'solid_angle': float(analyzer_obj.get_solid_angle()),
            'num_scatterers': num_scatterers,
        },
    }


def prediction_txt_from_parameters(parameters=None):
    result = predict_from_parameters(parameters)
    columns = [np.asarray(result['x'], dtype=float)]
    header_names = ['energy_loss_eV']
    for series in result['series']:
        header_names.append(series['name'].replace(' ', '_') + '_per_sec')
        columns.append(np.asarray(series['values'], dtype=float))
    data = np.column_stack(columns)
    handle = io.StringIO()
    np.savetxt(handle, data, header=' '.join(header_names))
    return handle.getvalue()


def save_prediction_txt(parameters=None, output_txt='prediction_output.txt'):
    if not output_txt:
        raise ValueError('Output path is required.')
    output_txt = os.path.expanduser(str(output_txt))
    output_dir = os.path.dirname(os.path.abspath(output_txt))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    text = prediction_txt_from_parameters(parameters)
    with open(output_txt, 'w') as handle:
        handle.write(text)
    return os.path.abspath(output_txt)


PREDICTION_WEB_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>XRS Prediction Platform</title>
  <script src="https://unpkg.com/vue@3/dist/vue.global.prod.js"></script>
  <style>
    :root {
      color-scheme: light;
      --ink: #1f2528;
      --muted: #65737a;
      --line: #d7dee2;
      --panel: #ffffff;
      --field: #f8fafb;
      --accent: #0f766e;
      --accent-strong: #0b5d57;
      --warn: #a15c00;
      --bg: #eef3f4;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      min-height: 100vh;
      background: var(--bg);
      color: var(--ink);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      letter-spacing: 0;
    }
    button, input, select { font: inherit; }
    #app { min-height: 100vh; }
    .shell {
      display: grid;
      grid-template-columns: minmax(280px, 360px) minmax(0, 1fr);
      min-height: 100vh;
    }
    aside {
      background: var(--panel);
      border-right: 1px solid var(--line);
      padding: 18px;
      overflow: auto;
    }
    main {
      padding: 20px;
      overflow: auto;
    }
    h1 {
      margin: 0 0 4px;
      font-size: 21px;
      font-weight: 720;
    }
    .subhead {
      color: var(--muted);
      font-size: 13px;
      margin-bottom: 18px;
    }
    .section {
      border-top: 1px solid var(--line);
      padding: 14px 0 4px;
    }
    .section h2 {
      margin: 0 0 10px;
      font-size: 13px;
      text-transform: uppercase;
      color: #405057;
    }
    .grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
    }
    label {
      display: grid;
      gap: 5px;
      color: var(--muted);
      font-size: 12px;
      min-width: 0;
    }
    input, select {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: var(--field);
      color: var(--ink);
      min-height: 36px;
      padding: 7px 9px;
    }
    input:focus, select:focus {
      border-color: var(--accent);
      outline: 2px solid rgba(15, 118, 110, 0.16);
    }
    .wide { grid-column: 1 / -1; }
    .toggle {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 10px;
      color: var(--ink);
      min-height: 38px;
    }
    .toggle input {
      width: 18px;
      min-height: 18px;
      accent-color: var(--accent);
    }
    .actions {
      display: flex;
      gap: 8px;
      margin-top: 14px;
    }
    button {
      border: 1px solid var(--accent);
      background: var(--accent);
      color: white;
      border-radius: 6px;
      min-height: 36px;
      padding: 7px 12px;
      cursor: pointer;
    }
    button.secondary {
      background: white;
      color: var(--accent-strong);
    }
    .status {
      min-height: 22px;
      margin-top: 12px;
      font-size: 13px;
      color: var(--muted);
    }
    .status.error { color: #b42318; }
    .chart-area {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      min-height: calc(100vh - 40px);
      padding: 16px;
      display: grid;
      grid-template-rows: auto minmax(360px, 1fr) auto;
      gap: 14px;
    }
    .chart-top {
      display: flex;
      justify-content: space-between;
      gap: 16px;
      align-items: flex-start;
      flex-wrap: wrap;
    }
    .chart-title {
      display: grid;
      gap: 2px;
    }
    .chart-title strong { font-size: 18px; }
    .chart-title span { color: var(--muted); font-size: 13px; }
    .metrics {
      display: grid;
      grid-template-columns: repeat(3, minmax(115px, 1fr));
      gap: 8px;
      min-width: min(100%, 470px);
    }
    .metric {
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 9px;
      background: #fbfcfc;
    }
    .metric span {
      display: block;
      color: var(--muted);
      font-size: 11px;
      margin-bottom: 3px;
    }
    .metric strong {
      display: block;
      font-size: 14px;
      overflow-wrap: anywhere;
    }
    .svg-wrap {
      min-height: 360px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #ffffff;
      overflow: hidden;
      position: relative;
    }
    svg {
      width: 100%;
      height: 100%;
      min-height: 360px;
      display: block;
    }
    .axis { stroke: #8b9aa1; stroke-width: 1; }
    .gridline { stroke: #edf1f3; stroke-width: 1; }
    .curve {
      fill: none;
      stroke-width: 2.4;
      vector-effect: non-scaling-stroke;
    }
    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px 14px;
      color: var(--muted);
      font-size: 12px;
      min-height: 18px;
    }
    .legend-item {
      display: inline-flex;
      align-items: center;
      gap: 6px;
    }
    .legend-swatch {
      width: 18px;
      height: 3px;
      border-radius: 999px;
      flex: 0 0 auto;
    }
    .axis-label {
      fill: var(--muted);
      font-size: 12px;
    }
    .empty {
      position: absolute;
      inset: 0;
      display: grid;
      place-items: center;
      color: var(--muted);
      pointer-events: none;
    }
    .table {
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      max-height: 220px;
      overflow-y: auto;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 12px;
    }
    th, td {
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      text-align: right;
      font-variant-numeric: tabular-nums;
    }
    th:first-child, td:first-child { text-align: left; }
    th {
      position: sticky;
      top: 0;
      background: #f7f9fa;
      color: #44535a;
    }
    @media (max-width: 900px) {
      .shell { grid-template-columns: 1fr; }
      aside { border-right: 0; border-bottom: 1px solid var(--line); }
      .chart-area { min-height: 520px; }
      .metrics { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <div id="app">
    <div class="shell">
      <aside>
        <h1>XRS Prediction</h1>
        <div class="subhead">Vue parameter platform</div>

        <div class="section">
          <h2>Files</h2>
          <div class="grid">
            <label class="wide">Import .inp<input type="file" accept=".inp,.txt" @change="importInp"></label>
            <label class="wide">Output TXT path<input v-model="savePath"></label>
          </div>
          <div class="actions">
            <button class="secondary" @click="downloadTxt">Download TXT</button>
            <button class="secondary" @click="saveTxt">Save on server</button>
          </div>
        </div>

        <div class="section">
          <h2>Sample</h2>
          <div class="grid">
            <label class="wide">Formula(s)<input v-model="form.chem_formulas"></label>
            <label>Concentration<input v-model="form.concentrations"></label>
            <label>Density g/cm3<input v-model="form.densities"></label>
            <label>Molar mass<input v-model="form.molar_masses"></label>
            <label>Thickness cm<input type="number" step="0.001" v-model.number="form.sample_thickness"></label>
            <label>2 theta deg<input type="number" step="0.1" v-model.number="form.angle_tth"></label>
          </div>
        </div>

        <div class="section">
          <h2>Energy</h2>
          <div class="grid">
            <label>E0 keV<input type="number" step="0.01" v-model.number="form.E0"></label>
            <label>Points<input type="number" min="20" max="2000" step="10" v-model.number="form.points"></label>
            <label>Loss min eV<input type="number" step="1" v-model.number="form.eloss_min"></label>
            <label>Loss max eV<input type="number" step="1" v-model.number="form.eloss_max"></label>
          </div>
        </div>

        <div class="section">
          <h2>Beam</h2>
          <div class="grid">
            <label class="wide">I0 photons/sec<input type="number" step="10000000000" v-model.number="form.i0_intensity"></label>
            <label>Height micron<input type="number" step="1" v-model.number="form.beam_height"></label>
            <label>Width micron<input type="number" step="1" v-model.number="form.beam_width"></label>
          </div>
        </div>

        <div class="section">
          <h2>Analyzer</h2>
          <div class="grid">
            <label>Material<input v-model="form.analyzer_material"></label>
            <label>HKL<input v-model="form.hkl"></label>
            <label>Mask mm<input type="number" step="1" v-model.number="form.mask_d"></label>
            <label>Bend radius m<input type="number" step="0.1" v-model.number="form.bend_r"></label>
            <label>Resolution eV<input type="number" step="0.01" v-model.number="form.energy_resolution"></label>
            <label>Efficiency<input type="number" step="0.01" min="0" max="1" v-model.number="form.analyzer_efficiency"></label>
            <label class="toggle wide">Calculate reflectivity<input type="checkbox" v-model="form.use_reflectivity"></label>
          </div>
        </div>

        <div class="section">
          <h2>Detector</h2>
          <div class="grid">
            <label>Material<input v-model="form.detector_material"></label>
            <label>Energy keV<input type="number" step="0.01" v-model.number="form.detector_energy"></label>
            <label class="wide">Thickness micron<input type="number" step="10" v-model.number="form.detector_thickness"></label>
          </div>
        </div>

        <div class="actions">
          <button @click="predict">Run</button>
          <button class="secondary" @click="reset">Reset</button>
        </div>
        <div class="status" :class="{ error: error }">{{ error || status }}</div>
      </aside>

      <main>
        <section class="chart-area">
          <div class="chart-top">
            <div class="chart-title">
              <strong>Absolute Counts</strong>
              <span>{{ subtitle }}</span>
            </div>
            <div class="metrics">
              <div class="metric"><span>Detector efficiency</span><strong>{{ metric('detector_efficiency') }}</strong></div>
              <div class="metric"><span>Analyzer efficiency</span><strong>{{ metric('analyzer_efficiency') }}</strong></div>
              <div class="metric"><span>Solid angle</span><strong>{{ metric('solid_angle') }}</strong></div>
            </div>
          </div>

          <div class="svg-wrap">
            <svg viewBox="0 0 1000 520" preserveAspectRatio="none" role="img">
              <line v-for="tick in yTicks" :key="'gy' + tick.y" class="gridline" :x1="pad.l" :x2="1000 - pad.r" :y1="tick.y" :y2="tick.y"></line>
              <line v-for="tick in xTicks" :key="'gx' + tick.x" class="gridline" :x1="tick.x" :x2="tick.x" :y1="pad.t" :y2="520 - pad.b"></line>
              <polyline v-for="curve in curves" :key="curve.name" class="curve" :points="curve.points" :stroke="curve.color"></polyline>
              <line class="axis" :x1="pad.l" :x2="1000 - pad.r" :y1="520 - pad.b" :y2="520 - pad.b"></line>
              <line class="axis" :x1="pad.l" :x2="pad.l" :y1="pad.t" :y2="520 - pad.b"></line>
              <text v-for="tick in yTicks" :key="'ly' + tick.y" class="axis-label" :x="pad.l - 10" :y="tick.y + 4" text-anchor="end">{{ tick.label }}</text>
              <text v-for="tick in xTicks" :key="'lx' + tick.x" class="axis-label" :x="tick.x" y="502" text-anchor="middle">{{ tick.label }}</text>
              <text class="axis-label" x="500" y="516" text-anchor="middle">energy loss (eV)</text>
              <text class="axis-label" x="16" y="260" transform="rotate(-90 16 260)" text-anchor="middle">counts / sec</text>
            </svg>
            <div v-if="!curves.length" class="empty">No data</div>
          </div>
          <div class="legend">
            <span v-for="curve in curves" :key="'legend-' + curve.name" class="legend-item">
              <span class="legend-swatch" :style="{ background: curve.color }"></span>
              <span>{{ curve.name }}</span>
            </span>
          </div>

          <div class="table">
            <table>
              <thead>
                <tr>
                  <th>Energy loss eV</th>
                  <th v-for="series in seriesList" :key="'head-' + series.name">{{ series.name }}</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="row in previewRows" :key="row.x">
                  <td>{{ fmt(row.x) }}</td>
                  <td v-for="series in seriesList" :key="row.x + '-' + series.name">{{ fmt(series.values[row.index]) }}</td>
                </tr>
              </tbody>
            </table>
          </div>
        </section>
      </main>
    </div>
  </div>

  <script>
    const defaults = __DEFAULTS__;
    const { createApp } = Vue;

    createApp({
      data() {
        return {
          form: { ...defaults },
          savePath: defaults.output_txt || 'prediction_output.txt',
          result: null,
          status: 'Ready',
          error: '',
          timer: null,
          requestId: 0,
          colors: ['#0f766e', '#b45309', '#2563eb', '#be123c', '#6d28d9', '#15803d'],
          pad: { l: 78, r: 24, t: 24, b: 48 }
        };
      },
      computed: {
        x() { return this.result?.x || []; },
        seriesList() { return this.result?.series || []; },
        pointsBySeries() {
          return this.seriesList.map((series, seriesIndex) => ({
            name: series.name,
            color: this.colors[seriesIndex % this.colors.length],
            points: this.x.map((xValue, index) => ({
              x: Number(xValue),
              y: Number(series.values[index])
            })).filter(point => Number.isFinite(point.x) && Number.isFinite(point.y))
          }));
        },
        bounds() {
          const points = this.pointsBySeries.flatMap(series => series.points);
          if (!points.length) return null;
          const xs = points.map(point => point.x);
          const ys = points.map(point => point.y);
          const xmin = Math.min(...xs), xmax = Math.max(...xs);
          let ymin = Math.min(...ys), ymax = Math.max(...ys);
          if (ymin === ymax) { ymin -= 1; ymax += 1; }
          const margin = (ymax - ymin) * 0.08;
          return { xmin, xmax, ymin: ymin - margin, ymax: ymax + margin };
        },
        curves() {
          if (!this.bounds) return [];
          const w = 1000 - this.pad.l - this.pad.r;
          const h = 520 - this.pad.t - this.pad.b;
          return this.pointsBySeries.map(series => ({
            name: series.name,
            color: series.color,
            points: series.points.map(point => {
              const px = this.pad.l + ((point.x - this.bounds.xmin) / (this.bounds.xmax - this.bounds.xmin || 1)) * w;
              const py = this.pad.t + (1 - ((point.y - this.bounds.ymin) / (this.bounds.ymax - this.bounds.ymin || 1))) * h;
              return `${px.toFixed(2)},${py.toFixed(2)}`;
            }).join(' ')
          })).filter(curve => curve.points);
        },
        xTicks() {
          if (!this.bounds) return [];
          return Array.from({ length: 6 }, (_, index) => {
            const ratio = index / 5;
            const value = this.bounds.xmin + (this.bounds.xmax - this.bounds.xmin) * ratio;
            return {
              x: this.pad.l + (1000 - this.pad.l - this.pad.r) * ratio,
              label: this.fmt(value)
            };
          });
        },
        yTicks() {
          if (!this.bounds) return [];
          return Array.from({ length: 5 }, (_, index) => {
            const ratio = index / 4;
            const value = this.bounds.ymax - (this.bounds.ymax - this.bounds.ymin) * ratio;
            return {
              y: this.pad.t + (520 - this.pad.t - this.pad.b) * ratio,
              label: this.fmt(value)
            };
          });
        },
        previewRows() {
          const rows = [];
          const step = Math.max(1, Math.floor(this.x.length / 12));
          for (let index = 0; index < this.x.length; index += step) {
            rows.push({ index, x: this.x[index] });
          }
          return rows.slice(0, 12);
        },
        subtitle() {
          const mode = this.result?.meta?.analyzer_efficiency_mode || 'manual';
          return `${this.form.chem_formulas} at ${this.form.angle_tth} deg, analyzer efficiency ${mode}`;
        }
      },
      watch: {
        form: {
          deep: true,
          handler() {
            clearTimeout(this.timer);
            this.timer = setTimeout(() => this.predict(), 350);
          }
        }
      },
      mounted() {
        this.predict();
      },
      methods: {
        async importInp(event) {
          const file = event.target.files?.[0];
          if (!file) return;
          this.status = `Importing ${file.name}...`;
          this.error = '';
          try {
            const text = await file.text();
            const response = await fetch('/api/import-inp', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ text })
            });
            const payload = await response.json();
            if (!response.ok) throw new Error(payload.error || 'Import failed');
            this.form = { ...defaults, ...payload.parameters };
            this.savePath = file.name.replace(/\.(inp|txt)$/i, '') + '_prediction.txt';
            this.status = `Imported ${file.name}`;
          } catch (error) {
            this.error = error.message;
            this.status = '';
          } finally {
            event.target.value = '';
          }
        },
        async predict() {
          const id = ++this.requestId;
          this.status = this.form.use_reflectivity ? 'Calculating reflectivity...' : 'Calculating...';
          this.error = '';
          try {
            const response = await fetch('/api/predict', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify(this.form)
            });
            const payload = await response.json();
            if (!response.ok) throw new Error(payload.error || 'Prediction failed');
            if (id !== this.requestId) return;
            this.result = payload;
            this.status = `Updated ${new Date().toLocaleTimeString()}`;
          } catch (error) {
            if (id !== this.requestId) return;
            this.error = error.message;
            this.status = '';
          }
        },
        reset() {
          this.form = { ...defaults };
          this.savePath = defaults.output_txt || 'prediction_output.txt';
        },
        async saveTxt() {
          this.status = 'Saving...';
          this.error = '';
          try {
            const response = await fetch('/api/save-prediction', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify({ parameters: this.form, output_txt: this.savePath })
            });
            const payload = await response.json();
            if (!response.ok) throw new Error(payload.error || 'Save failed');
            this.status = `Saved ${payload.path}`;
          } catch (error) {
            this.error = error.message;
            this.status = '';
          }
        },
        downloadTxt() {
          if (!this.result) {
            this.error = 'Run a prediction before downloading.';
            return;
          }
          const text = this.resultText();
          const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
          const link = document.createElement('a');
          link.href = URL.createObjectURL(blob);
          link.download = (this.savePath || 'prediction_output.txt').split(/[\\/]/).pop() || 'prediction_output.txt';
          document.body.appendChild(link);
          link.click();
          URL.revokeObjectURL(link.href);
          link.remove();
          this.status = `Prepared ${link.download}`;
        },
        resultText() {
          const names = ['energy_loss_eV', ...this.result.series.map(item => `${item.name.replace(/\s+/g, '_')}_per_sec`)];
          const lines = [`# ${names.join(' ')}`];
          for (let index = 0; index < this.result.x.length; index += 1) {
            const row = [this.result.x[index], ...this.result.series.map(item => item.values[index])];
            lines.push(row.map(value => Number(value).toExponential(18)).join(' '));
          }
          return lines.join('\n') + '\n';
        },
        metric(key) {
          const value = this.result?.meta?.[key];
          return value === undefined ? '-' : this.fmt(value);
        },
        fmt(value) {
          if (!Number.isFinite(Number(value))) return '-';
          const number = Number(value);
          if (Math.abs(number) >= 1e4 || (Math.abs(number) > 0 && Math.abs(number) < 1e-3)) {
            return number.toExponential(3);
          }
          return number.toLocaleString(undefined, { maximumSignificantDigits: 5 });
        }
      }
    }).mount('#app');
  </script>
</body>
</html>
"""


class PredictionRequestHandler(BaseHTTPRequestHandler):
    def _send_json(self, payload, status=200):
        body = json.dumps(payload).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self):
        html = PREDICTION_WEB_HTML.replace('__DEFAULTS__', json.dumps(WEB_DEFAULT_PARAMETERS))
        body = html.encode('utf-8')
        self.send_response(200)
        self.send_header('Content-Type', 'text/html; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path in ('/', '/index.html'):
            self._send_html()
        elif path == '/api/defaults':
            self._send_json(WEB_DEFAULT_PARAMETERS)
        else:
            self._send_json({'error': 'Not found'}, status=404)

    def do_POST(self):
        path = urlparse(self.path).path
        try:
            length = int(self.headers.get('Content-Length', '0'))
            raw_body = self.rfile.read(length) if length else b'{}'
            payload = json.loads(raw_body.decode('utf-8'))
            if path == '/api/predict':
                self._send_json(predict_from_parameters(payload))
            elif path == '/api/import-inp':
                self._send_json({'parameters': web_parameters_from_inp_text(payload.get('text', ''))})
            elif path == '/api/save-prediction':
                output_path = payload.get('output_txt') or WEB_DEFAULT_PARAMETERS['output_txt']
                saved_path = save_prediction_txt(payload.get('parameters', {}), output_path)
                self._send_json({'path': saved_path})
            else:
                self._send_json({'error': 'Not found'}, status=404)
        except Exception as exc:
            self._send_json({'error': str(exc)}, status=400)

    def log_message(self, format, *args):
        return


def serve_web(host='127.0.0.1', port=8765):
    server = ThreadingHTTPServer((host, port), PredictionRequestHandler)
    print('XRS prediction platform running at http://%s:%s' % (host, port))
    print('Press Ctrl+C to stop the server.')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\nStopping XRS prediction platform.')
    finally:
        server.server_close()


def run(filename='prediction.inp', output_txt=None):
    """
    Function to create a spectrum prediction from input parameters provided in the input file filename.
    Generates a figure with the result.
    If output_txt is provided, save the calculated absolute counts to that text file.
    """
    # parse all input parameters
    inp = get_all_input(filename)
    # create all the instances
    beam_obj            = beam(inp['beam']['i0_intensity'],inp['beam']['beam_height'],inp['beam']['beam_width'],
                   inp['beam']['divergence'])
    sample_obj          = sample(inp['sample']['chem_formulas'],inp['sample']['concentrations'],inp['sample']['densities'],
                     inp['sample']['angle_tth'],inp['sample']['sample_thickness'],inp['sample']['angle_in'],
                     inp['sample']['angle_out'],inp['sample']['shape'],inp['sample']['molar_masses'])
    analyzer_obj        = analyzer(inp['analyzer']['material'],inp['analyzer']['hkl'],inp['analyzer']['mask_d'],
                       inp['analyzer']['bend_r'],inp['analyzer']['energy_resolution'], inp['analyzer']['diced'],
                       inp['analyzer']['thickness'],inp['analyzer']['database_dir'])
    detector_obj        = detector(inp['detector']['energy'],inp['detector']['thickness'],inp['detector']['material'],inp['detector']['pixel_size'])
    compton_profile_obj = compton_profiles(sample_obj,inp['compton_profiles']['eloss_range'],inp['compton_profiles']['E0'])
    thomson_obj         = thomson(compton_profile_obj.get_energy_in_keV(),compton_profile_obj.get_E0(),compton_profile_obj.get_tth())

    abs_cross_section_obj = absolute_cross_section(beam_obj, sample_obj, analyzer_obj, detector_obj, thomson_obj, compton_profile_obj)
    if output_txt:
        abs_cross_section_obj.save_txt(output_txt, header='energy_loss_eV absolute_counts_per_sec')
        print('Saved prediction results to ' + output_txt)
    abs_cross_section_obj.plot_abs_cross_section()
    return abs_cross_section_obj


def build_arg_parser():
    parser = argparse.ArgumentParser(description='Create an XRS spectrum prediction.')
    parser.add_argument('-f', '--file', dest='filename', default='prediction.inp', help='read input from FILE')
    parser.add_argument('-o', '--output-txt', dest='output_txt', default=None, help='save calculated prediction to FILE')
    parser.add_argument('--web', action='store_true', help='start the Vue parameter input platform')
    parser.add_argument('--host', default='127.0.0.1', help='web server host for --web')
    parser.add_argument('--port', type=int, default=8765, help='web server port for --web')
    return parser


def main(argv=None):
    options = build_arg_parser().parse_args(argv)
    if options.web:
        serve_web(options.host, options.port)
        return None
    return run(options.filename, options.output_txt)


if __name__ == '__main__':
    main()
