#!/usr/bin/python
# Filename: theory2.py
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

from . import xrs_utilities
import numpy as np
import math
import shelve
import os
import pylab
import h5py


from scipy import interpolate, signal, integrate, constants, optimize
from re import findall
# from pylab import *
import pylab
from optparse import OptionParser

__metaclass__ = type # new style classes

installation_dir = os.path.dirname(os.path.abspath(__file__))

def cla():
        pass

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
            self.shape == 'sphere' # list of tth only for sphere geometry
            ac = (mu_tot_in + mu_tot_out)/(1.0 - np.exp(-mu_tot_in*thickness -mu_tot_out*thickness))

        else:
            if self.shape == 'slab' and alpha and beta:
                if tth:
                    if self.beta<0: # transmission geometry
                        test_tth = alpha-beta
                    else: # reflection geometry
                        test_tth = 180.0 - (alpha+beta)
                    if tth == test_tth:
                        pass
                    else:
                        print( 'the alpha and beta values set are not congruent to the tth value set!')
                absorption_correction_factor = xrs_utilities.abscorr2(mu_tot_in,mu_tot_out,alpha,beta,thickness)
            elif self.shape == 'sphere':
                ac = (mu_tot_in + mu_tot_out)/(1.0 - np.exp(-mu_tot_in*thickness -mu_tot_out*thickness)) 
                #1.0/np.exp(-thickness*mu_tot_in -thickness*mu_tot_out) # spherical sample just add up in and outgoing absorption
            else:
                print( 'please provide either shape=\'sphere\' (default) or \'slab\' and alpha and beta!')
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
        #
        # for attr in ['eloss', 'J', 'C', 'V', 'q', 'thomson', 'I0', 'beam_h', 'beam_v', 'sample_tth',
        #                 'sample_densities', 'sample_molar_masses', 'sample_formulas', 'sample_thickness',
        #                 'sample_concentrations', 'sample_abs_in', 'det_efficiency', 'ana_efficiency',
        #                 'ana_solid_angle', 'ana_energy_resolution', 'energy_in_keV' ]:
        #     h5_group[attr] = eval( 'self.' + attr )
        f.flush()
        f.close()      

def input_file_parser(filename):
    """
    Parses an input file, which has a structure like the example input file ('prediction.inp') provided in the
    examples/ folder. (Python lists and numpy arrays have to be profived without white spaces in their definitions,
    e.g. 'hkl = [6,6,0]' instead of 'hkl = [6, 6, 0]')
    """
    try:
        lines = open(filename,'r').readlines()
        #f = open(filename,'r')
    except IOError:
        print( 'No input file ' + filename + ' found.')
        return
    input_parameters = {} # dictionary of input parameters
    section_names = ['detector','analyzer','sample','thomson','beam','compton_profiles']
    for name in section_names:
        input_parameters[name] = {} # cempty list for all sections, fill them up from file, add defaults for missing values
    # parse all given parameters
    lineindex = 0
    while lineindex < len(lines):
        line = lines[lineindex]
        if line[0:4] == '####':
            thekey = line.split()[1]
            lineindex += 1
            while lines[lineindex][0:4] != '####' and lineindex < len(lines):
                if not lines[lineindex] == '\n':
                    input_parameters[thekey][lines[lineindex].split()[0]] = eval(lines[lineindex].split()[2])
                    lineindex += 1
                else:
                    lineindex += 1
                if lineindex == len(lines):
                    break
    return input_parameters, section_names

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


import sys
if(sys.argv[0][-12:]!="sphinx-build"):

    # parse input arguments, i.e. input file
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", default="prediction.inp", help="read input from FILE", metavar="FILE")
    parser.add_option("-o", "--output-txt", dest="output_txt", default=None, help="save calculated prediction to FILE", metavar="FILE")
    (options, args) = parser.parse_args()



if __name__ == '__main__':
    run(options.filename, options.output_txt)


def main():
    run(options.filename, options.output_txt)

#######################################################
# plot q vs. matrix elements
#######################################################

class radial_wave_function:
        def __init__(self):
            self.element = None
            self.Z       = None
            self.n       = None
            self.l       = None
            self.s       = None
            self.hydrogen_like  = False
            self.spin_polarized = False
            self.R_nl_numeric   = np.array([])
            self.r              = np.array([])

        def load_from_sympy( self, Z, n, l ):
            try:
                import importlib
                hydrogen = importlib.import_module('sympy.physics.hydrogen')
                R_nl = hydrogen.R_nl
                sympy = importlib.import_module('sympy')
                var = sympy.var
            except ImportError:
                print('Did not find sympy module, will end here.')
                return
            rr, ZZ = var("rr ZZ")
            R_nl_function = R_nl(n, l, rr, ZZ)
            r = np.logspace(1e-6, 2., 1000)-1.0
            R_nl_numeric = np.zeros_like(r)
            for ii in range(r.shape[0]):
                R_nl_numeric[ii] = R_nl(n, l, rr, ZZ).evalf(subs={rr: r[ii] , ZZ:Z})
            self.r = r
            self.R_nl_numeric = R_nl_numeric
            self.n = n
            self.l = l
            self.Z = Z
            self.element = xrs_utilities.element(Z)
            self.hydrogen_like  = True
            self.spin_polarized = False

        def load_from_PP( self, Z, n, l, path='/home/christoph/programs/atomic_wavefunctions/',
                            spin_polarized=False ):
            l_str = ['s', 'p', 'd', 'f', 'g', 'h'][l]
            fname_up = path + xrs_utilities.element(Z).lower() + '/' + 'ae/wf-'+str(n)+l_str+'_up'
            fname_dn = path + xrs_utilities.element(Z).lower() + '/' + 'ae/wf-'+str(n)+l_str+'_dn'
            # wfcn in form u(r) = R(r)/r
            # for checking norm: calculate dr r^2 u(r) u*(r)
            raw_up = np.loadtxt(fname_up)
            raw_dn = np.loadtxt(fname_dn)

class matrix_element:
    def __init__( self, R1, R2 ):
        self.wfn1 = R1.R_nl_numeric
        self.wfn2 = R2.R_nl_numeric
        self.k    = np.array([])
        self.r    = np.linspace(0.0, 15.0, 1000)
        self.Mel  = np.array([])
        self.q    = np.array([])

    def compute( self, k ):
        """ **compute**
        Calculates the matrix elements for a given k or range of k.
        """
        all_k = []
        if not isinstance(k,list):
            all_k.append(k)
        else:
            all_k = k
        self.k   = np.array(all_k)
        self.q   = np.zeros_like(self.r)
        self.Mel = np.zeros((len(self.r),len(k)))
        for ii in range(len(k)):
            self.q, self.Mel[:,ii] = xrs_utilities.compute_matrix_elements( self.wfn1, self.wfn2, all_k[ii], self.r )

    def write_H5( self, filename ):
        """ **write_H5**
        Creates an HDF5 file to store the matrix elements.

        Args:
          * fname (str) : Full path and filename for the HDF5 file to be created.
        """
        if np.any(self.q):
            # check if file already exists
            if os.path.isfile( filename ):
                os.remove( filename )
            f = h5py.File( filename, "w" )
            f.require_group( "matrix_elements" )
            f["matrix_elements"]["r"] = self.r
            f["matrix_elements"]["q"] = self.q
            f["matrix_elements"]["M"] = self.Mel
            f.close()
        else:
            print('There are no matrix elements to save.')

    def write_ascii( self, filename ):
        """ **write_ascii**
        Creates an ascii-file and writes matrix elements.

        Args:
          * fname (str) : Full path and filename for the ascii file to be created.
        """
        if np.any(self.q):
            # check if file already exists
            if os.path.isfile( filename ):
                os.remove( filename )
            the_data = np.zeros((len(self.q), len(self.k)+1))
            the_data[:,0] = self.q
            for ii in range(len(self.k)):
                the_data[:,ii+1] = self.Mel[:,ii]

            np.savetxt( filename, the_data )













