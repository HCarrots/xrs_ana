#!/usr/bin/python
# Filename: xrs_utilities.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group and contains practical functions, 
# most of which are translated from Matlab functions from the University of
# Helsinki Electronic Structure Laboratory.
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

import os
import math
import copy

import numpy as np
import matplotlib.pyplot as plt
import pickle
import traceback
import sys

from matplotlib.widgets import Cursor
from itertools import groupby
from scipy.integrate import trapezoid
from scipy import interpolate, signal, integrate, constants, optimize
from re import findall
from scipy.ndimage import measurements
from scipy.optimize import leastsq, fmin, fsolve, minimize
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.integrate import odeint



# data_installation_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)),"..","..","..","..","share","xrstools","data")
# data_installation_dir = os.path.abspath('.')

# whne you test the file from its source directory you, /data sits on level above. In this case you can
# work around it by creating a link in ./ to ../data
data_installation_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)),"resources",  'data')

# os.path.join(getattr(install_cmd, 'install_lib'),"xrstools"+version,"..","..","..","..","share","xrstools","data")


def diode(current, energy, thickness=0.03):
    """ **diode**
    Calculates the number of photons incident for a Si PIPS diode.

    Args:
      * current (float): Diode current in [pA].
      * energy (float): Photon energy in [keV].
      * thickness (float): Thickness of Si active layer in [cm].

    Returns:
      * flux (float): Number of photons per second.

    Function adapted from Matlab function by S. Huotari.
    """
    t = thickness # thickness of diode in cm
    
    # Total cross-section of absorbed energy for Si (Storm & Israel)
    sx = np.array([2,3,4,5,6,8,10,15,20,30,40,50,60,80,100,150,200])
    s  = np.array([125000,44600,20600,11000,6600,2890,1510,448,186, \
                       53.3,21.9,11.2, 6.61,3.18,2.06,1.41,1.34])

    my = np.exp(np.interp(energy, sx , np.log(s)))
    my = (0.02144*2.32)*my

    n_ph= energy*(1.0-np.exp(-tauphoto(14,energy)*2.32*t))
    n_ph=n_ph + (energy)*(1.0-np.exp(-sigmainc(14,energy)*2.32*t))
    n_ph=n_ph/0.0036

    n_ph = energy*(1.0-np.exp(-my*t))/0.0036 # number of electrons/photon
    n_pA = current/1.6022e-7                     # number of electrons per pA
    print('The photon flux is: %E'%(n_pA/n_ph))
    return  n_pA/n_ph

def cshift(w1, th):
    """ **cshift**
    Calculates Compton peak position.

    Args:
      *  w1 (float, array): Incident energy in [keV].
      *  th (float): Scattering angle in [deg].

    Returns:
      *  w2 (foat, array): Energy of Compton peak in [keV].

    Funktion adapted from Keijo Hamalainen.

    """
    return w1/(1+w1/510.967*(1-np.cos(th/180*np.pi)))


def tauphoto(Z, energy, logtablefile=os.path.join(data_installation_dir,'logtable.dat')):
    """ **tauphoto**
    Calculates Photoelectric Cross Section in cm^2/g using Log-Log Fit.

    Args:
      * z (int or string): Element number or elements symbol.
      * energy (float or array): Energy (can be number or vector)

    Returns:
      * tau (float or array): Photoelectric cross section in [cm**2/g]

    Adapted from original Matlab function of Keijo Hamalainen.
    """

    en = np.array([])
    en = np.append(en,energy)
    logtable = np.loadtxt(logtablefile)
    # find the right places in logtable
    if not isinstance(Z,int):
        Z = element(Z)
    try:
        ind = list(logtable[:,0]).index(Z)
    except:
        print( 'no such element in logtable.dat')

    c = np.array(logtable[ind:ind+5,:]) # 5 lines that corresponds to the element
    tau_i = np.zeros((4, len(en)))
    for ii in range(4):
        for jj in range(4):
            tau_i[ii,:] = tau_i[ii,:] + c[jj+1,ii]*np.log(en)**(jj)
    tau2 = np.zeros(len(en))
    tau2 = (en<c[0,1])*1*tau_i[0,:] + \
            np.logical_and(en>=c[0,1] , en<c[0,2])*1*tau_i[1,:] + \
            np.logical_and(en>=c[0,2] , en<c[0,3])*1*tau_i[2,:] + \
            (en>=c[0,3])*1*tau_i[3,:]
    return np.exp(tau2)*0.6022/c[0,4]

def sigmainc(Z, energy, logtablefile=os.path.join(data_installation_dir,'logtable.dat')):
    """ **sigmainc**
    Calculates the Incoherent Scattering Cross Section 
    in cm^2/g using Log-Log Fit.

    Args:
      * z (int or string): Element number or elements symbol.
      * energy (float or array): Energy (can be number or vector)

    Returns:
      * tau (float or array): Photoelectric cross section in [cm**2/g]

    Adapted from original Matlab function of Keijo Hamalainen.

    """

    en = np.array([])
    en = np.append(en,energy)
    logtable = np.loadtxt(logtablefile)
    # find the right places in logtable
    if not isinstance(Z,int):
        Z = element(Z)
    try:
        ind = list(logtable[:,0]).index(Z)
    except:
        print( 'no such element in logtable.dat')

    c = np.array(logtable[ind:ind+5,:]) # 5 lines that corresponds to the element
    sigmai=0
    for jj in range(4):
        sigmai = sigmai + c[jj+1,5]*np.log(energy)**(jj)

    return np.exp(sigmai)*0.6022/c[0,4]

def Rx(chi, degrees=True):
    """ **Rx**
    Rotation matrix for vector rotations around the [1,0,0]-direction.

    Args:
      * chi   (float) : Angle of rotation.
      * degrees(bool) : Angle given in radians or degrees.

    Returns:
       * 3x3 rotation matrix.
    """
    if degrees:
        chi = np.radians(chi)
    return np.array([[1,0,0],[0, np.cos(chi), -np.sin(chi)], [0, np.sin(chi), np.cos(chi)]])

def Ry(phi, degrees=True):
    """ **Ry**
    Rotation matrix for vector rotations around the [0,1,0]-direction.

    Args:
      * phi   (float) : Angle of rotation.
      * degrees(bool) : Angle given in radians or degrees.

    Returns:
       * 3x3 rotation matrix.
    """
    if degrees:
        phi = np.radians(phi)
    return np.array([[np.cos(phi), 0, np.sin(phi)],[0, 1, 0],[-np.sin(phi), 0, np.cos(phi)]])

def Rz(omega, degrees=True):
    """ **Rz**
    Rotation matrix for vector rotations around the [0,0,1]-direction.

    Args:
      * omega (float) : Angle of rotation.
      * degrees(bool) : Angle given in radians or degrees.

    Returns:
      * 3x3 rotation matrix.
    """
    if degrees:
        omega = np.radians(omega)
    return np.array([[np.cos(omega), -np.sin(omega), 0],[np.sin(omega), np.cos(omega), 0],[0,0,1]])


def Phi(phi, degrees=True):
    """ rotation around (0,1,0), neg sense
    """
    if degrees:
        phi = np.radians(phi)
    return np.array([[np.cos(phi), 0, np.sin(phi)],[0, 1, 0],[-np.sin(phi), 0, np.cos(phi)]])
    
def Chi(chi, degrees=True):
    """ rotation around (1,0,0), pos sense
    """
    if degrees:
        chi = np.radians(chi)
    return np.array([[1,0,0],[0, np.cos(chi), -np.sin(chi)], [0, np.sin(chi), np.cos(chi)]])
    
def Omega(omega, degrees=True):
    """ rotation around (0,0,1), pos sense
    """
    if degrees:
        omega = np.radians(omega)
    return np.array([[np.cos(omega), -np.sin(omega), 0],[np.sin(omega), np.cos(omega), 0],[0,0,1]])

def get_UB_Q(tthv, tthh, phi, chi, omega, **kwargs):
    """ **get_UB_Q**
    Returns the momentum transfer and scattering vectors for given
    FOURC spectrometer and sample angles. U-, B-matrices and 
    incident/scattered wavelength are passed as keyword-arguments.

    Args:
      * tthv (float): Spectrometer vertical 2Theta angle.
      * tthh (float): Spectrometer horizontal 2Theta angle.
      * chi (float): Sample rotation around x-direction.
      * phi (float): Sample rotation around y-direction.
      * omega (float): Sample rotation around z-direction.

      * kwargs (dict): Dictionary with key-word arguments:
          * kwargs['U'] (array): 3x3 U-matrix Lab-to-sample transformation.
          * kwargs['B'] (array): 3x3 B-matrix reciprocal lattice 
            to absolute units transformation.
          * kwargs['lambdai'] (float): Incident x-ray wavelength in Angstrom.
          * kwargs['lambdao'] (float): Scattered x-ray wavelength in Angstrom.

    Returns:
      * Q_sample  (array): Momentum transfer in sample coordinates.
      * Ki_sample (array): Incident beam direction in sample coordinates.
      * Ko_sample (array): Scattered beam direction in sample coordinates.

    """
    U       = kwargs['U']
    B       = kwargs['B']
    Lab     = kwargs['Lab']
    beam_in = kwargs['beam_in']
    lambdai = kwargs['lambdai']
    lambdao = kwargs['lambdao']
    # scattering vectors in laboratory frame
    Ki_test = 2.0*np.pi/lambdai * beam_in/np.linalg.norm(beam_in)
    Ko_test = 2.0*np.pi/lambdao * Rz(tthh).dot(Ry(tthv)).dot(beam_in/np.linalg.norm(beam_in))
    # h_lab = Omega*Chi*Phi*U*B*h_cryst (Busing equ. 19)
    # invert
    # h_cryst = B_inv*U_inv*Phi_inv*Chi_inv*Omega_inv*h_lab
    Q_test    = Ko_test - Ki_test
    Phi_inv   = Phi(phi).T #np.linalg.inv(Phi(phi))
    print(Phi_inv)
    Chi_inv   = Chi(chi).T #np.linalg.inv(Chi(chi))
    Omega_inv = Omega(omega).T #np.linalg.inv(Omega(omega))
    U_inv     = np.linalg.inv(U)
    B_inv     = np.linalg.inv(B)
    Q_sample  = np.matmul(B_inv ,np.matmul(U_inv , np.matmul( Phi_inv , np.matmul(Chi_inv, np.matmul( Omega_inv, Q_test)))))
    Ki_sample = np.matmul(B_inv ,np.matmul(U_inv , np.matmul( Phi_inv , np.matmul(Chi_inv, np.matmul( Omega_inv, Ki_test)))))
    Ko_sample = np.matmul(B_inv ,np.matmul(U_inv , np.matmul( Phi_inv , np.matmul(Chi_inv, np.matmul( Omega_inv, Ko_test)))))
    return Q_sample, Ki_sample, Ko_sample
    
def find_diag_angles(q, x0, U, B, Lab, beam_in, lambdai, lambdao, tol=1e-8, method='BFGS'):
    """ **find_diag_angles**
    Finds the FOURC spectrometer and sample angles for a desired q.

    Args:
      * q (array): Desired momentum transfer in Lab coordinates.
      * x0 (list): Guesses for the angles (tthv, tthh, chi, phi, omega).
      * U (array): 3x3 U-matrix Lab-to-sample transformation.
      * B (array): 3x3 B-matrix reciprocal lattice to absolute units transformation.
      * lambdai (float): Incident x-ray wavelength in Angstrom.
      * lambdao (float): Scattered x-ray wavelength in Angstrom.
      * tol (float): Toleranz for minimization (see scipy.optimize.minimize)
      * method (str): Method for minimization (see scipy.optimize.minimize)

    Returns:
       * ans (array): tthv, tthh, phi, chi, omega
    """
    # put UB matrix and energies into keyword argument for minimization
    kwargs = {'U': U, 'B': B, 'Lab': Lab, 'beam_in':beam_in, 'lambdai': lambdai, 'lambdao': lambdao}
    # least square minimization between wanted and guessed q
    fitfctn = lambda x: np.sum(( q - get_UB_Q(x[0], x[1], x[2], x[3], x[4], \
                                                  **kwargs)[0]  )**2)
    ans=minimize(fitfctn, x0, bounds=((0.,0.),(-10.,110.),(-7.,7.),(-7.,7.),(None,None)), tol=tol, method=method)
    print( ans )
    return ans.x

def get_gnuplot_rgb( start=None, end=None, length=None ):
    """ **get_gnuplot_rgb**
    Prints out a progression of RGB hex-keys to use in Gnuplot.

    Args:
      * start (array): RGB code to start from (must be numbers out of [0,1]).
      * end   (array): RGB code to end at (must be numbers out of [0,1]).
      * length  (int): How many colors to print out.
    """
    if start==None and end==None and length==None:
        rgb =  [[0,      0,           1], [0,      0.2353, 0.7647], \
                [0,      0.4706, 0.5294], [0,      0.7059, 0.2941], \
                [0,      0.9412, 0.0588], [0.1765,  0.8235, 0    ], \
                [0.4118,  0.5882, 0    ], [0.6471,  0.3529, 0    ], \
                [1,      0,      0     ] ]
    else:
        rgb = np.zeros((length,3))
        rgb[:,0] = np.linspace(start[0], end[0], length)
        rgb[:,1] = np.linspace(start[1], end[1], length)
        rgb[:,2] = np.linspace(start[2], end[2], length)
    #
    for ii in range(len(rgb)):
        print( 'set style line %d lt -1 lc rgb \'#%02x%02x%02x\' lw 1.5'%(ii+1, rgb[ii][0]*255., rgb[ii][1]*255., rgb[ii][2]*255.) )

class maxipix_det:
    """
    Class to store some useful values from the detectors used. To be used for arranging the ROIs.
    """
    def __init__(self,name,spot_arrangement):
        self.name = name
        assert spot_arrangement in ['3x4','vertical'], 'unknown ROI arrangement, select \'3x4\' or \'vertical\'.'
        self.spot_arrangement = spot_arrangement

        self.roi_indices   = []
        self.roi_x_indices = []
        self.roi_y_indices = []
        self.roi_x_means   = []
        self.roi_y_means   = []
        self.pixel_size    = [256,256]
        self.PIXEL_RANGE   = {'VD': [0,256,0,256],  'VU': [0,256,256,512],  'VB': [0,256,512,768],
                            'HR': [256,512,0,256],'HL': [256,512,256,512],'HB': [256,512,512,768]}
    def get_pixel_range(self):
        return self.PIXEL_RANGE[self.name]

    def get_det_name(self):
        return self.name


class bragg_refl:
    """
    Dynamical theory of diffraction.
    """
    def __init__(self, crystal, hkl, alpha=0.0 ):

        # constants
        self.hc = constants.h*constants.c
        self.C  = 1.0
        self.r0 = constants.e**2/ \
          (4.0*np.pi*constants.epsilon_0* \
            constants.m_e*constants.c**2)* \
            1.0e10 # classical electron radius expressed in Angstrom
        self.P = 1 # polarization factor
        self.alpha = alpha
            
        # params
        self.hkl     = hkl
        self.crystal = crystal
        self.dspace  = dspace( self.hkl, self.crystal )
        self.ff_energy, self.f1_energy, self.f2_energy = \
          self.get_nff()

    def get_reflectivity_bent(self, energy, delta_theta, R):
        _, theta_B_deg, _ = self.get_reflectivity(energy, delta_theta)
        theta_B = np.radians(theta_B_deg)
        refl,e,dev,e0 = taupgen( energy, self.hkl, crystals=self.crystal, R=R,
                    dev=delta_theta, alpha=self.alpha )
        return refl, theta_B*180/np.pi
        
    def get_reflectivity(self, energy, delta_theta, case='sigma'):
        energy = energy*1e3
        wavelength = self.hc/(energy*constants.e)*1e10
        print(wavelength)
        theta_B = np.arcsin(wavelength/(2.0*self.dspace)) # Bragg angle
        theta_inc = theta_B + np.radians(self.alpha) # incidence angle
        theta_ref = theta_B - np.radians(self.alpha) # reflection angle
        b = -np.sin(theta_ref)/np.sin(theta_inc) # asymmetry factor
        # polarization factor
        P = self.get_polarization_factor(2.*theta_B, case=case)
        # chi
        chi_0, chi_h, chi_hbar = self.get_chi(energy, self.crystal, self.hkl)
        #    index of refraction
        n = 1 + chi_0/2
        n_delta = 1 - np.real(n)
        n_beta = -np.imag(n)
        # linear absorption coefficient
        mu = -2*np.pi*np.imag(chi_0)/wavelength 
        # Bragg angle correction
        theta_B_correction = -chi_0*(1 - b)/(2*np.sin(2*theta_B)) 
        # width of the total reflection domain
        delta = np.abs(P)*np.sqrt(np.abs(b)*chi_h*chi_hbar)/np.sin(2*theta_B) 
        Darwin_width = 2*np.real(delta)
        lambda_B = wavelength*np.abs(np.sin(theta_inc))/ \
          (2*np.pi*np.real(delta)*np.sin(2*theta_B))
        delta_theta = delta_theta/1e6 # delta_theta is input in microrad
        eta = (delta_theta - theta_B_correction)/delta # deviation parameter
        # reflectivity curve
        reflectivity_curve = np.abs(eta - np.sign(np.real(eta))* \
                                        np.sqrt(eta**2 - 1))**2
        return reflectivity_curve, theta_B*180/np.pi, theta_B_correction
        
    def get_nff(self,
        nff_path = os.path.join(data_installation_dir,'atomic_form_factors')):
        fname    = os.path.join(nff_path, self.crystal.lower()+'.nff')
        table = np.loadtxt(fname, unpack = True, skiprows = 1)
        table = np.transpose(table)
        ff_energy = table[:,0]
        f1_energy = table[:,1]
        f2_energy = table[:,2]
        return ff_energy, f1_energy, f2_energy

    def get_chi(self, energy, crystal=None, hkl=None):
        path       =  os.path.join(data_installation_dir,'chitable_')
        hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
        filestring = path + crystal.lower() + hkl_string + '.dat'
        self.chi        = np.loadtxt(filestring)
        self.chi_0 = complex(np.interp(energy, self.chi[:,0], self.chi[:,1]), \
                    np.interp(energy, self.chi[:,0],self.chi[:,2]))
        self.chi_h = complex(np.interp(energy, self.chi[:,0], self.chi[:,3]), \
                    np.interp(energy, self.chi[:,0], self.chi[:,4]))
        self.chi_hbar = np.conj(self.chi_h)
        return self.chi_0, self.chi_h, self.chi_hbar

    def get_polarization_factor(self, tth, case='sigma'):
        """
        Calculate polarization factor.
        """
        if case == 'sigma':
            P = 1.0
            self.P = P
            return P
        elif case == 'pi':
            P = np.cos(tth)**2
            self.P = P
            return P
        elif case == None:
            P = (1 + np.cos(tth)**2)/2.0
            self.P = P
            return P

        
        
class dtxrd:
    """ class to hold all things dynamic theory of diffraction.
    """
    def __init__(self, hkl, energy, crystal='Si', asym_angle=0.0, angular_range=[-0.5e-3, 0.5e-3] , angular_step=1e-8 ):

        # constants
        self.hc = 12.3984191
        self.C  = 1.0

        # params
        self.hkl     = hkl
        self.energy  = energy
        self.lam     = self.hc/self.energy
        self.crystal = crystal
        self.dspace  = dspace( self.hkl, self.crystal )

        # load Chi from tables:
        path       =  os.path.join(data_installation_dir,'chitable_')
        hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
        filestring = path + crystal.lower() + hkl_string + '.dat'
        self.chi        = np.loadtxt(filestring)
        self.chi0 = complex(np.interp(self.energy, self.chi[:,0], self.chi[:,1]), \
                    np.interp(self.energy, self.chi[:,0],self.chi[:,2]))
        self.chih = complex(np.interp(self.energy, self.chi[:,0], self.chi[:,3]), \
                    np.interp(self.energy, self.chi[:,0], self.chi[:,4]))
        self.chihbar = np.conj(self.chih)

        # set Bragg angle:
        self.thetab  = float( bragg( hkl, energy , xtal=crystal ) )
        self.thetabd = float( braggd( hkl, energy , xtal=crystal ) )

        # set asymmetry:
        self.set_asymmetry(asym_angle)

        # set reduced deviation parameter
        self.angular_range = angular_range
        self.angular_step  = angular_step
        self.get_eta(angular_range, angular_step)

        self.lam_ext = None
        self.mu0     = None
        self.mus     = None
        self.R       = None
        self.omega_h = None
        self.omega_0 = None

    def set_asymmetry(self, alpha):
        """
        negative alpha -> more grazing incidence
        """
        self.alpha  = alpha
        self.gammah = -np.sin(self.thetab + alpha) # Krisch et al. convention 
        self.gamma0 =  np.sin(self.thetab - alpha) # Krisch et al. convention 
        self.gamma  =  self.gamma0/np.abs(self.gammah)                             
        self.beta   =  self.gammah/self.gamma0

    def set_hkl(self, hkl):
        self.hkl = hkl

    def set_energy(self, energy):
        self.energy = energy

    def get_reflectivity(self, angular_range=None, angular_step=None):
        if not angular_range:
            angular_range = self.angular_range
        if not angular_step:
            angular_step = self.angular_step

        self.get_eta(angular_range, angular_step)
        
        pre_factor = (np.sqrt(self.chih*self.chihbar))/(self.chihbar) * \
          1.0/np.sqrt(np.abs(self.gamma)) * self.gammah/np.abs(self.gammah) * self.C/np.abs(self.C)
          
        self.R = np.abs(pre_factor * np.piecewise(self.eta, self.eta.real>0, \
                                          [lambda x: x - np.sqrt( x**2 + self.gammah/np.abs(self.gammah)  )  , \
                                            lambda x: x + np.sqrt( x**2 + self.gammah/np.abs(self.gammah)  ) ]))
            
    def get_eta(self, angular_range, angular_step=1e-8):
        self.theta = np.arange( self.thetab+angular_range[0], self.thetab+angular_range[1], angular_step )
        self.eta = ( (self.theta-self.thetab)*np.sin(2.*self.thetab) + self.chi0/2.0*(1-self.gamma)  ) / \
          (np.abs(self.C)*np.sqrt(np.abs(self.gamma))* np.sqrt(self.chih*self.chihbar))

    def get_anomalous_absorption(self, energy=None):
        if not energy:
            energy = self.energy

        # photoelectric absorption
        mu0 = myprho( energy, self.crystal )[0][0][0]*0.602252/ \
          myprho( energy, self.crystal )[2]

        omega = np.arctan( (1.0 -  np.abs(self.R)**2) / (1.0 + np.abs(self.R)**2) * np.tan(self.thetab) )

        mus = mu0 * np.cos(omega)/np.cos(self.thetab) * ( 1.0 + 2.0*np.abs(self.C) * \
                    self.chih.imag/self.chi0.imag * self.R.real/(1.0 + np.abs(self.R)**2)  )
        self.mus = mus

    def get_extinction_length(self, energy=None):

        if energy:
            lam = energy/self.hc
        else:
            lam = self.lam

        kappa = -np.abs( (self.chih.imag)/(self.chih.real) )
        
        self.lam_ext = (lam*np.sqrt(self.gamma0 * np.abs(self.gammah))* np.cos(kappa) ) / \
          (2.0*np.pi*np.abs(self.C) * np.sqrt( np.abs(self.chih.real*self.chihbar.real)  ))

    def get_reflection_width(self):
        if not self.lam_ext:
            self.get_extinction_length()
        
        self.omega_h = self.lam*np.abs(self.gammah)/( np.pi*self.lam_ext*np.sin(2.0*self.thetab) )
        self.omega_0 = self.lam*self.gamma0/( np.pi*self.lam_ext*np.sin(2.0*self.thetab) )

"""
def dtxrd_reflectivity( energy, hkl, alpha=0.0, crystal='Si', angular_range=np.arange(-0.5e-3, 0.5e-3) ):

    # constants
    hc, lam, dsp, thetab, alpha, C = \
      get_dtxrd_constants( energy, hkl, crystal=crystal, alpha=alpha )   

    # theta scale
    theta  = np.arange( thetab+angular_range[0], \
                        thetab+angular_range[1], 1e-8 )

    # load chi
    chi0, chih, chihbar = get_dtxrd_chi( energy, hkl, crystal=crystal )

    # asymmetry parameter
    gamma, gamma0, gammah, beta = get_dtxrd_assymmetry_params(thetab, alpha)

    # calculate eta
    eta = get_dtxrd_eta(theta, thetab, chi0, chih, chihbar)

    # calculate reflectivity
    pre_factor = (np.sqrt(chih*chihbar))/(chihbar) * \
        1.0/np.sqrt(np.abs(gamma)) * gammah/np.abs(gammah) * C/np.abs(C)
    r = pre_factor * np.piecewise(eta, eta.real>0, \
        [lambda x: x - np.sqrt( x**2 + gammah/np.abs(gammah)  )  , \
         lambda x: x + np.sqrt( x**2 + gammah/np.abs(gammah)  ) ])

    return theta, r
"""   

"""
def dtxrd_anomalous_absorption( energy, hkl, alpha=0.0, crystal='Si', angular_range=np.arange(-0.5e-3, 0.5e-3) ):

    # constants
    hc     = 12.3984191
    lam    = hc/energy * 1e-10
    dsp    = dspace(hkl, crystal)
    thetab = float( bragg( hkl, energy , xtal=crystal ) )
    alpha  = np.radians(alpha)
    C      = 1.0

    # theta scale
    theta  = np.arange( thetab+angular_range[0], \
                        thetab+angular_range[1], 1e-8 )

    # load chi
    path       = os.path.join(data_installation_dir,'chitable_')
    hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
    filestring = path + crystal.lower() + hkl_string + '.dat'
    chi        = np.loadtxt(filestring)
    chi0 = complex(np.interp(energy,chi[:,0],chi[:,1]),
                   np.interp(energy,chi[:,0],chi[:,2]))
    chih = complex(np.interp(energy,chi[:,0],chi[:,3]),
                   np.interp(energy,chi[:,0],chi[:,4]))
    chihbar = np.conj(chih)
    
    # asymmetry parameter
    gammah = -np.sin(thetab + alpha) # Krisch et al. convention
    gamma0 = np.sin(thetab - alpha) # Krisch et al. convention
    beta   = gamma0/np.abs(gammah)
    gamma  = gammah/gamma0

    # calculate eta
    eta = ( (theta-thetab)*np.sin(2.*thetab) + chi0/2.0*(1-gamma)  ) / \
    (np.abs(C)*np.sqrt(np.abs(gamma))* np.sqrt(chih*chihbar))

    # calculate reflectivity
    pre_factor = (np.sqrt(chih*chihbar))/(chihbar) * \
        1.0/np.sqrt(np.abs(gamma)) * gammah/np.abs(gammah) * C/np.abs(C)
    r = pre_factor * np.piecewise(eta, eta.real>0, \
        [lambda x: x - np.sqrt( x**2 + gammah/np.abs(gammah)  )  , \
         lambda x: x + np.sqrt( x**2 + gammah/np.abs(gammah)  ) ])

    # photoelectric absorption
    mu0 = myprho( energy, crystal )[0][0][0]*0.602252/ \
        myprho( energy, crystal )[2]
    mus = mu0 * np.cos(omega)/np.cos(thetab) * ( 1.0 + 2.0*np.abs(C) * chih.imag/chi0.imag * r.real/(1.0 + np.abs(r)**2)  )
    return theta, mus
"""
def dtxrd_extinction_length( energy, hkl, alpha=0.0, crystal='Si' ):
    pass

def delE_dicedAnalyzerIntrinsic(E, Dw, Theta):
    """Calculates the intrinsic energy resolution of a diced crystal
    analyzer.

    Args:
        E     (float): Working energy in [eV].
        Dw    (float): Darwin width of the used reflection [microRad].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Intrinsic energy resolution of a perfect analyzer crystal.
    """
    Dw = Dw/1000000.0 # conversion to radians
    return E * (Dw)/(np.tan(np.radians(Theta)))

def delE_JohannAberration(E, A, R, Theta):
    """Calculates the Johann aberration of a spherical analyzer crystal.

    Args:
        E     (float): Working energy in [eV].
        A     (float): Analyzer aperture [mm].
        R     (float): Radius of the Rowland circle [mm].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Johann abberation in [eV].
    """
    return E/2.0 * ((A)/(2.0*R*np.tan(np.radians(Theta))))**2

def delE_pixelSize(E, p, R, Theta):
    """Calculates the pixel size contribution to the resolution function
    of a diced analyzer crystal.

    Args:
        E     (float): Working energy in [eV].
        p     (float): Pixel size in [mm].
        R     (float): Radius of the Rowland circle [mm].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Pixel size contribution in [eV] to the energy resolution for a diced analyzer
        crystal.
    """
    return E * (p/(2.0*R*np.sin(np.radians(Theta))))/(np.tan(np.radians(Theta)))

def delE_sourceSize(E, s, R, Theta):
    """Calculates the source size contribution to the resolution function.

    Args:
        E     (float): Working energy in [eV].
        s     (float): Source size in [mm].
        R     (float): Radius of the Rowland circle [mm].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Source size contribution in [eV] to the energy resolution.
    """
    return E * (s/(R*np.sin(np.radians(Theta))))/(np.tan(np.radians(Theta)))

def delE_offRowland(E, z, A, R, Theta):
    """Calculates the off-Rowland contribution of a spherical analyzer crystal.

    Args:
        E     (float): Working energy in [eV].
        z     (float): Off-Rowland distance [mm].
        A     (float): Analyzer aperture [mm].
        R     (float): Radius of the Rowland circle [mm].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Off-Rowland contribution in [eV] to the energy resolution.

    """
    return E * (z*A)/( (R*np.sin(np.radians(Theta)) + z)**2 ) * (1.0)/(np.tan(np.radians(Theta)))

def delE_stressedCrystal(E, t, v, R, Theta):
    """Calculates the stress induced contribution to the resulution function
    of a spherically bent crystal analyzer.

    Args:
        E     (float): Working energy in [eV].
        t     (float): Absorption length in the analyzer material [mm].
        v     (float): Poisson ratio of the analyzer material.
        R     (float): Radius of the Rowland circle [mm].
        Theta (float): Analyzer Bragg angle [degree].

    Returns:
        Stress-induced contribution in [eV] to the energy resolution.

    """
    return E * t/R * np.absolute((1.0)/(np.tan(np.radians(Theta))**2) - 2.0*v)

def get_num_of_MD_steps(time_ps,time_step):
    """Calculates the number of steps in an MD simulation for a desired 
    time (in ps) and given step size (in a.u.)

    Args:
        time_ps   (float): Desired time span (ps).
        time_step (float): Chosen time step (a.u.).

    Returns:
        The number of steps required to span the desired time span.
    """
    return time_ps / time_step /1.0e12 / 2.418884326505e-17

def nonzeroavg(y=None):
    dim1 = y.shape[0]
    yavg = np.zeros(dim1,float)
    for ii in range(dim1): 
        length = 0
        rowsum = 0.
        for ij in range(y.shape[1]):        
            if (not np.isnan(y[ii, ij])):                
                length += 1
                rowsum += y[ii, ij]            
        yavg[ii] = rowsum / float(length)
    yavg = yavg * y.shape[1]
    return(yavg)

def fermi(rs):
    """ **fermi**
    Calculates the plasmon energy (in eV), Fermi energy (in eV), Fermi 
    momentum (in a.u.), and critical plasmon cut-off vector (in a.u.).

    Args:
     * rs (float): electron separation parameter

    Returns:
       * wp (float): plasmon energy (in eV)
       * ef (float): Fermi energy (in eV)
       * kf (float): Fermi momentum (in a.u.)
       * kc (float): critical plasmon cut-off vector (in a.u.)

    Based on Matlab function from A. Soininen.
    """
    au   = 27.212
    alfa = (9.0*np.pi/4.0)**(1.0/3.0)
    kf = alfa/rs
    ef = kf*kf/2.0
    wp = np.sqrt(3.0/rs/rs/rs)
    kc = kf * (np.sqrt(1.0+wp/ef)-1.0)
    wp = wp*au
    ef = ef*au
    return wp, ef, kf, kc

def lindhard_pol(q,w,rs=3.93,use_corr=False, lifetime=0.28):
    """ **lindhard_pol**
    Calculates the Lindhard polarizability function (RPA) for 
    certain q (a.u.), w (a.u.) and rs (a.u.).

    Args:
      * q (float): momentum transfer (in a.u.)
      * w (float): energy (in a.u.)
      * rs (float): electron parameter
      * use_corr (boolean): if True, uses Bernardo's calculation for n(k) instead of the Fermi function.
      * lifetime (float): life time (default is 0.28 eV for Na).

    Based on Matlab function by S. Huotari.
    """
    if type(w) in [float, int]:
        w = np.array([w])
    wp, ef, kf, kc = fermi(rs)
    ef     = ef/27.212
    gammal = lifetime/27.212  # lifetime  (0.28 eV for Na)
    th     = np.arange( 0.0, np.pi, np.pi/700.0 )
    k      = np.arange( 0.0, 2.0*kf, kf/1000.0 )
    [K,TH] = np.meshgrid(k,th)

    ek     = K**2/2.0
    ekq    = ( K**2+q**2+2*q*K*np.cos(TH) )/2.0
    if not use_corr:
        fek = np.zeros(np.shape(ek)) 
        fek[ek<=ef]=1.0
        fekq=np.zeros(np.shape(ekq)) 
        fekq[ekq<=ef] = 1.0
    if use_corr:
        print('Not implemented yet!')
    x = np.zeros_like(w, dtype='complex')
    for ii in range(len(w)):
        y=np.sin(TH)*(fek-fekq)/(w[ii]+ek-ekq+np.complex(0,1)*gammal)
        y=np.trapezoid( y, th, axis=0 )
        y=np.trapezoid( k**2.0*y, k, axis=0 )
        x[ii]=y
    x = 4.0*np.pi*x
    x = x/(2.0*np.pi)**3
    return x 

def energy(d,ba):
    """
    % ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
    %         function e=energy(dspace,bragg_angle)
    %
    %      dspace for reflection
    %      bragg_angle in DEG
    %
    %         KH 28.09.93
    """
    hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
    return (2.0*d*np.sin(ba/180.0*np.pi)/hc)**(-1)

def dspace(hkl=[6,6,0],xtal='Si'):
    """
    % DSPACE Gives d-spacing for given xtal
    %     d=dspace(hkl,xtal)
    %     hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
    %     xtal='Si','Ge','LiF','InSb','C','Dia','Li' (case insensitive)
    %     if xtal is number this is user as a d0
    %
    %     KH 28.09.93 
    %        SH 2005
    %
    """
    # create a database of lattice constants (could be a shelf)
    xtable = {}
    xtable['SI'] = 5.43102088
    xtable['GE'] = 5.657
    xtable['SIXOP'] = 5.430919
    xtable['SIKOH'] = 5.430707
    xtable['LIF'] = 4.027
    xtable['INSB'] = 6.4784
    xtable['C'] = 6.708
    xtable['DIA'] = 3.57
    xtable['LI'] = 3.41
    xtable['TCNE'] = 9.736
    xtable['CU'] = 3.61
    xtable['PB'] = 4.95
    xtable['NA'] = 4.2906
    xtable['AL'] = 4.0495

    if isinstance(xtal,str):
        try:
            a0 = xtable[xtal.upper()]
        except KeyError:
            print( 'Lattice constant is not in database')
            return
    else: 
        a0 = xtal # if number is provided, it's taken as lattice constant

    return a0/np.sqrt(np.sum(np.array(hkl)**2.0))

def bragg(hkl,e,xtal='Si'):
    """
    % BRAGG  Calculates Bragg angle for given reflection in RAD
    %      output=bangle(hkl,e,xtal)
    %        hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
    %      e=energy in keV
    %      xtal='Si', 'Ge', etc. (check dspace.m) or d0 (Si default)
    %
    %      KH 28.09.93
    %
    """
    hc = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
    return np.real(np.arcsin((2.0*dspace(hkl,xtal)*e/hc)**(-1.0)))

def braggd(hkl,e,xtal='Si'):
    """
    # BRAGGD  Calculates Bragg angle for given reflection in deg
    #      Call BRAGG.M
    #      output=bangle(hkl,e,xtal)
    #        hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
    #      e=energy in keV
    #      xtal='Si', 'Ge', etc. (check dspace.m) or d0 (Si default)
    #
    #      KH 28.09.93
    """
    return bragg(hkl,e,xtal)/np.pi*180.0

def addch(xold,yold,n,n0=0,errors=None):
    """
    # ADDCH     Adds contents of given adjacent channels together
    #
    #           [x2,y2] = addch(x,y,n,n0)
    #           x  = original x-scale  (row or column vector)
    #           y  = original y-values (row or column vector)
    #           n  = number of channels to be summed up
    #            n0 = offset for adding, default is 0
    #           x2 = new x-scale 
    #           y2 = new y-values
    #
    #           KH 17.09.1990
    #        Modified 29.05.1995 to include offset
    """
    n0=int(n0-np.fix(n0/n)*n)
    if n0<0:
         n0 = (n + n0)
    datalen = np.floor( (len(xold) - n0) / n)

    xnew = np.zeros(int(np.min([datalen,len(xold)])))
    ynew = np.zeros(int(np.min([datalen,len(xold)])))
    errnew = np.zeros(int(np.min([datalen,len(xold)])))

    for i in range(int(datalen)):
        xnew[i] = np.sum(xold[i*n+n0:i*n+n+n0])/n
        ynew[i] = np.sum(yold[i*n+n0:i*n+n+n0])/n
        if np.any(errors):
            errnew[i] = np.sqrt(np.sum(errors[i*n+n0:i*n+n+n0]**2.0))
            return xnew, ynew, errnew
    return xnew, ynew

def fwhm(x,y):
    """
    finds full width at half maximum of the curve y vs. x
    returns 
    f  = FWHM
    x0 = position of the maximum
    """
    if x[-1] < x[0]:
        x = np.flipud(x)
        y = np.flipud(y)

    y0 = np.amax(y)
    i0 = np.where(y == y0)
    x0 = x[i0]

    i1 = np.where(np.logical_and(y>y/3.0, x<x0))[0]
    i2 = np.where(np.logical_and(y>y/3.0, x>x0))[0]

    if len(y[i1])==0 or len(y[i2])==0:
        return 0,0
    #f  = interpolate.interp1d(y[i1],x[i1], bounds_error=False, fill_value=0.0)
    #x1 = f(y0/2.0)
    #f  = interpolate.interp1d(y[i2],x[i2], bounds_error=False, fill_value=0.0)
    #x2 = f(y0/2.0)
    x1 = np.interp(y0/2.0,y[i1],x[i1])
    x2 = np.interp(y0/2.0,np.flipud(y[i2]),np.flipud(x[i2]))
    fwhm = x2 - x1
    x0 = np.mean([x2, x1])
    return fwhm, x0

def gauss(x,x0,fwhm):
    # area-normalized gaussian
    sigma = fwhm/(2*np.sqrt(2*np.log(2)));
    y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
    return y

def convg(x,y,fwhm):
    """
    Convolution with Gaussian
    x  = x-vector
    y  = y-vector
    fwhm = fulll width at half maximum of the gaussian with which y is convoluted
    """
    dx = np.min(np.absolute(np.diff(x)))
    x2 = np.arange(np.min(x)-1.5*fwhm, np.max(x)+1.5*fwhm, dx)
    xg = np.arange(-np.floor(2.0*fwhm/dx)*dx, np.floor(2.0*fwhm/dx)*dx, dx)
    yg = gauss(xg,0,fwhm)
    yg = yg/np.sum(yg)
    y2 = spline2(x,y,x2)
    c  = np.convolve(y2,yg, mode='full')
    n  = int( np.floor(np.max(np.shape(xg))/2))
    c  = c[n:len(c)-n+1] # not sure about the +- 1 here
    f  = interpolate.interp1d(x2,c)
    return f(x)

def interpolate_M(xc, xi, yi, i0):
    """
    Linear interpolation scheme after Martin Sundermann that conserves
    the absolute number of counts.

    ONLY WORKS FOR EQUALLY/EVENLY SPACED XC, XI! 

    Args:
        xc (np.array): The x-coordinates of the interpolated values.
        xi (np.array): The x-coordinates of the data points, must be increasing.
        yi (np.array): The y-coordinates of the data points, same length as `xp`.
        i0 (np.array): Normalization values for the data points, same length as `xp`.

    Returns:
        ic (np.array): The interpolated and normalized data points.

from scipy.interpolate import Rbf
x = arange(20)
d = zeros(len(x))
d[10] = 1
xc = arange(0.5,19.5)
rbfi = Rbf(x, d)
di = rbfi(xc)


    """
    assert len(xi)==len(yi) and len(xi)==len(i0), "xi, yi, and i0 must have the same length."

    xc        = np.array(xc)
    xi        = np.array(xi)
    yi        = np.array(yi)
    i0        = np.array(i0)
    dx        = (xi-xc[:len(xi)])*(len(xc)-1)/(xc[-1]-xc[0])
    yc        = 0.0*xc
    ic        = 0.0*xc
    for i in np.unique(np.floor(np.sort(dx))):
        dxi    = [x-i if(((x-i)>0) and ((x-i)<1)) else -1 for x in dx]
        if((i>=0) and (i+i+len(xi)+1)<=len(xc)): # if yi and i0 lay inside the grid range add them
            yc[i:i+len(xi)]        = yc[i:i+len(xi)]        + yi*(1-np.absolute(dxi))
            ic[i:i+len(xi)]        = ic[i:i+len(xi)]        + i0*(1-np.absolute(dxi))
            yc[i+1:i+len(xi)+1]    = yc[i+1:i+len(xi)+1]    + yi*np.maximum(dxi,0)
            ic[i+1:i+len(xi)+1]    = ic[i+1:i+len(xi)+1]    + i0*np.maximum(dxi,0)
        else:
            if (min(i+len(xi),len(xc))-max(i,0)) >= 0: # if yi and i0 lay only partial inside the grid range add only the overlapping region
                yc[max(i,0):min(i+len(xi),len(xc))]        = yc[max(i,0):min(i+len(xi),len(xc))]        + (yi*(1-np.absolute(dxi)))[max(-i,0):len(xi)+min(len(xc)-i-len(xi),0)]
                ic[max(i,0):min(i+len(xi),len(xc))]        = ic[max(i,0):min(i+len(xi),len(xc))]        + (i0*(1-np.absolute(dxi)))[max(-i,0):len(xi)+min(len(xc)-i-len(xi),0)]
            if (min(i+len(xi)+1,len(xc))-max(i+1,0)) >= 0:
                yc[max(i+1,0):min(i+len(xi)+1,len(xc))]    = yc[max(i+1,0):min(i+len(xi)+1,len(xc))]    + (yi*np.maximum(dxi,0))[max(-i-1,0):len(xi)+min(len(xc)-i-len(xi)-1,0)]
                ic[max(i+1,0):min(i+len(xi)+1,len(xc))]    = ic[max(i+1,0):min(i+len(xi)+1,len(xc))]    + (i0*np.maximum(dxi,0))[max(-i-1,0):len(xi)+min(len(xc)-i-len(xi)-1,0)]
    return yc, ic

def spline2(x,y,x2):
    """
    Extrapolates the smaller and larger valuea as a constant
    """
    xmin = np.min(x)
    xmax = np.max(x)
    imin = x == xmin
    imax = x == xmax
    f  = interpolate.interp1d(x,y, bounds_error=False, fill_value=0.0)
    y2 = f(x2)
    i     = np.where(x2<xmin)
    y2[i] = y[imin]
    i     = np.where(x2>xmax)
    y2[i] = y[imax]
    return y2

def pz2e1(w2,pz,th):
    """Calculates the incident energy for a specific scattered photon and momentum value.

    Returns the incident energy for a given photon energy and scattering angle.
    This function is translated from Keijo Hamalainen's Matlab implementation (KH 29.05.96).

    Args:
      * w2 (float): scattered photon energy in [keV]
      * pz (np.array): pz scale in [a.u.]
      * th (float): scattering angle two theta in [deg]

    Returns:
      * w1 (np.array): incident energy in [keV]
    """
    pz  = np.array(pz)
    w   = np.array(np.arange(np.array(w2)/4.0,4.0*np.array(w2),np.array(w2)/5000.0))
    p   = e2pz(w,w2,th)[0]
    if ( p[1]-p[0] <0) :
        tck = interpolate.UnivariateSpline(p[::-1],w[::-1])
    else:
        tck = interpolate.UnivariateSpline(p,w)
    w1  = tck(pz)
    return w1

def e2pz(w1,w2,th):
    """Calculates the momentum scale and the relativistic Compton cross section 
    correction according to P. Holm, PRA 37, 3706 (1988).

    This function is translated from Keijo Hamalainen's Matlab implementation (KH 29.05.96).

    Args:
      * w1 (float or np.array): incident energy in [keV]
      * w2 (float or np.array): scattered energy in [keV]
      * th (float): scattering angle two theta in [deg]
    returns:
      * pz (float or np.array): momentum scale in [a.u.]
      * cf (float or np.array): cross section correction factor such that: J(pz) = cf * d^2(sigma)/d(w2)*d(Omega) [barn/atom/keV/srad]
    """
    w1  = np.array(w1)    # make sure arrays are used
    w2  = np.array(w2)           
    m   = constants.value('electron mass energy equivalent in MeV')*1e3 #511.003      # Mass in natural units
    th  = math.radians(th) # th/180.0*np.pi  # Angle to radians
    alp = constants.value('fine-structure constant') #1.0/137.036  # Fine structure constant
    r0  = constants.value('classical electron radius') #2.8179e-15   # Electron radius
    q   = np.sqrt(w1**2.0 + w2**2.0-2.0*w1*w2*np.cos(th))                        # Momentum transfer    
    pz  = q/2.0 - (w1-w2) * np.sqrt(1.0/4.0 + m**2.0/(2.0*w1*w2*(1.0-np.cos(th)))) # In natural units
    E   = np.sqrt(m**2.0+pz**2.0)
    A   = ((w1-w2)*E-w1*w2*(1.0-np.cos(th)))/q
    D   = (w1-w2*np.cos(th))*A/q
    R   = w1*(E-D)
    R2  = R-w1*w2*(1-np.cos(th))
    chi = R/R2 + R2/R + 2.0*m**2.0 * (1.0/R-1.0/R2) + m**4.0 * (1.0/R-1.0/R2)**2.0
    cf  = 2.0*w1*q*E/(m**2.0*r0**2.0*w2*chi)
    cf  = cf*(1.0e-28*(m*alp)) # Cross section now in barns/atom/keV/srad
    pz  = pz/(m*alp)           # pz to atomic units (a.u.)
    return pz, cf

def momtrans_au(e1,e2,tth):
    """ Returns the momentum transfer (in a.u.).

    Calculates the momentum transfer in atomic units for two given
    energies e1 and e1 (in keV) and the scattering angle tth (two theta).

    Args:
      *e1 (float or np.array): incident energy in [keV], can be a single value or a vector
      *e2 (float or np.array): scattered energy in [keV], can be a single value or a vector
      *tth (float): scattering angle two theta in [deg]

    Returns:
      * q (float or np.array): momentum transfer [a.u.], single value or vector depending on input
    """
    e1    = np.array(e1*1.0e3/13.60569172/2.0)
    e2    = np.array(e2*1.0e3/13.60569172/2.0)
    th    = math.radians(tth)#tth/180.0*np.pi
    hbarc = 137.03599976
    q     = 1/hbarc*np.sqrt(e1**2.0+e2**2.0-2.0*e1*e2*np.cos(th));
    return q

def vrot(v,vaxis,phi):
    """ **vrot**
    Rotates a vector around a given axis.

    Args:
      * v (np.array): vector to be rotated
      * vaxis (np.array): rotation axis
      * phi (float): angle [deg] respecting the right-hand rule 

    Returns:
     * v2 (np.array): new rotated vector

    Function by S. Huotari (2007) adopted to Python.
    """
    h = vaxis[0]
    k = vaxis[1]
    l = vaxis[2]
    alpha = np.arctan2(k,h)
    if np.absolute(alpha)>np.finfo(float).eps:
        h2 = np.cos(alpha)*(h+k*np.tan(alpha))
    else:
        h2 = h
    v2 = np.array([h2, 0.0, l])
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    R1 = np.array([[ca, sa, 0.0], [-sa, ca, 0.0], [0.0, 0.0, 1.0]])
    beta = np.radians(vangle(v2,np.array([0.0, 0.0, 1.0])))
    cb = np.cos(beta)
    sb = np.sin(beta)
    R2 = np.array([[cb, 0.0, -sb], [0.0, 1.0, 0.0], [sb, 0.0, cb]])
    phi = np.radians(phi)
    cp = np.cos(phi)
    sp = np.sin(phi)
    R3 = np.array([[cp, -sp, 0.0], [sp, cp, 0.0], [0.0, 0.0, 1.0]])
    v2 = np.dot(R3,np.dot(R2,np.dot(R1,v))) 
    v2 = np.dot(np.linalg.inv(R1),np.dot(np.linalg.inv(R2),v2))
    return v2    

def vrot2(vector1,vector2,angle):
    """ **rotMatrix**
    Rotate vector1 around vector2 by an angle.
    """
    theta = np.radians(angle)
    R=np.array([[vector2[0]**2+(1.0-vector2[0]**2)*np.cos(theta), (1.0-np.cos(theta))*vector2[0]*vector2[1]-np.sin(theta)*vector2[2], (1.0-np.cos(theta))*vector2[0]*vector2[2]+np.sin(theta)*vector2[1]], [(1.0-np.cos(theta))*vector2[0]*vector2[1]+np.sin(theta)*vector2[2], vector2[1]**2+(1.0-vector2[1]**2)*np.cos(theta), (1.0-np.cos(theta))*vector2[1]*vector2[2]-np.sin(theta)*vector2[0]],[(1.0-np.cos(theta))*vector2[0]*vector2[2]-np.sin(theta)*vector2[1], (1.0-np.cos(theta))*vector2[1]*vector2[2]+np.sin(theta)*vector2[0], vector2[2]**2+(1.0-vector2[2]**2)*np.cos(theta)]])
    return np.dot(R,vector1)

def vangle(v1, v2):
    """ **vangle**
    Calculates the angle between two cartesian vectors v1 and v2 in degrees.

    Args:
      * v1 (np.array): first vector.
      * v2 (np.array): second vector.

    Returns:
      * th (float): angle between first and second vector.

    Function by S. Huotari, adopted for Python.
    """
    return np.arccos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))/np.pi*180.0;

def convtoprim(hklconv):
    """ **convtoprim**
    converts diamond structure reciprocal lattice expressed in conventional
    lattice vectors to primitive one (Helsinki -> Palaiseau conversion)
    from S. Huotari
    """
    return hklconv[2]*np.array([0.5,0.5,0.0]) + hklconv[1]*np.array([0.5,0.0,0.5]) + hklconv[0]*np.array([0.0,0.5,0.5])

def primtoconv(hklprim):
    """ **primtoconv**
    converts diamond structure reciprocal lattice expressed in primitive basis
    to the conventional basis (Palaiseau -> Helsinki conversion)
    from S. Huotari
    """
    a = np.array([0.0, 0.5, 0.5])
    b = np.array([0.5, 0.0, 0.5])
    c = np.array([0.5, 0.5, 0.0])
    Gp = np.linalg.inv([a,b,c]).T
    ap = Gp[0,:]
    bp = Gp[1,:]
    cp = Gp[2,:]
    return hklprim[0]*ap + hklprim[1]*bp + hklprim[2]*cp

def householder(b,k):
    """
    function H = householder(b, k)
    % H = householder(b, k)
    % Atkinson, Section 9.3, p. 611
    % b is a column vector, k an index < length(b)
    % Constructs a matrix H that annihilates entries
    % in the product H*b below index k

    % $Id: householder.m,v 1.1 2008-01-16 15:33:30 mike Exp $
    % M. M. Sussman
    """
    n = len(b)
    d = b[k:n]

    if d[0] >= 0.0:
        alpha = -np.linalg.norm(d)
    else:
        alpha = np.linalg.norm(d)

    if alpha == 0.0:
        H = np.eye(n)
        return

    lenD = len(d)
    v = np.zeros(lenD)

    v[0] = np.sqrt(0.5*(1.0-d[0]/alpha))
    p = -alpha*v[0]
    v[1:lenD] = d[1:lenD]/(2.0*p)
    w = np.append(  np.zeros((k,1)) ,v).reshape(n,1)
    H = np.eye(n)-2.0 * np.dot(w,w.T)
    return H

def svd_my(M,maxiter=100,eta=0.1):
    sind = 0
    import copy
    import scipy as sp

    # initialize U,S,V
    X = copy.deepcopy(M)
    m,n = np.shape(X)
    k   = np.amin([m,n])
    U   = np.random.rand(m,k)
    V   = np.random.rand(n,k)
    S   = np.random.rand(k,k)

    # orthogonalize U,V
    #U = sp.linalg.orth(U)
    #V = sp.linalg.orth(V)

    # compute S
    #S = np.dot(np.dot(U.T,X),V)

    # compute cost J0
    J0 = 0.5*np.linalg.norm(X - np.dot(np.dot(U,S),V.T) )**2
    J  = J0
    dJ = J
    while sind <= maxiter:
        sind += 1
        # update U and V
        U = U + eta*(np.dot(X,V) + U.dot(V.T).dot(X.T).dot(U)  ).dot(S)
        V = V + eta*(np.dot(X.T,U) + V.dot(U.T).dot(X).dot(V) ).dot(S)
        # compute S
        S = U.T.dot(X).dot(V)
        # make S_ii positive
        V = np.dot(V,np.sign(S))
        S = np.abs(S)
        Jnew = 0.5*np.linalg.norm(X - np.dot(np.dot(U,S),V.T) )**2
        dJ   = Jnew - J
        J    = Jnew
        print( Jnew)
    return U,S,V

def bidiag_reduction(A):
    """
    function [U,B,V]=bidiag_reduction(A)
    % [U B V]=bidiag_reduction(A)
    % Algorithm 6.5-1 in Golub & Van Loan, Matrix Computations
    % Johns Hopkins University Press
    % Finds an upper bidiagonal matrix B so that A=U*B*V'
    % with U,V orthogonal.  A is an m x n matrix
    """
    import copy
    m,n = np.shape(A)
    B = copy.deepcopy(A)
    U = np.eye(m)
    V = np.eye(n)
    for k in range(n):
        # eliminate non-zeros below the diagonal
        H = householder(B[:,k],k)
        B = np.dot(H,B)
        U = np.dot(U,H)
        # eliminate non-zeros to the right of the 
        # superdiagonal by working with the transpose
        if k<n-1:
            H = householder(B[k,:].T,k+1)
            B = np.dot(B,H.T)
            V = np.dot(H,V)
    return U, B, V

def cixsUBgetQ_primo(tthv, tthh, psi):
    G = np.array([-1.,-1.,-1.])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    lambdao = hc/wo
    wi = wo
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U
    th = braggd(G,wo)
    xxx = vrot(np.array([0.0,-1.0, 1.0]),np.array([-2.0,1.0,1.0]),th)
    yyy = vrot(np.array([-2.0,1.0, 1.0]),np.array([-2.0,1.0,1.0]),th)
    zzz = vrot(G,np.array([-2.0,1.0,1.0]),th)
    #xxx = vrot(np.array([2.0,-1.0,-1.0]),np.array([0.0,-1.0,1.0]),th)
    #yyy = vrot(np.array([0.0,-1.0, 1.0]),np.array([0.0,-1.0,1.0]),th)
    #zzz = vrot(G,np.array([0.0,-1.0,1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # axis of rotation of psi
    v = np.array([-np.sin(np.radians(th)), 0.0, np.cos(np.radians(th))])
    Ki_test = 2.0*np.pi/lambdai*X
    Ko_test = 2.0*np.pi/lambdao*vrot(vrot(X,Y,-tthv) ,Z, tthh)
    Q_test = np.dot(np.linalg.lstsq(B,U)[0],vrot(Ki_test-Ko_test,v,-psi))
    return Q_test, Ki_test, Ko_test

def cixsUBgetAngles_primo(Q):
    G = np.array([-1.,-1.,-1.])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    print( wo)
    lambdao = hc/wo
    wi = wo + 0.10
    print( wi)
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U for Si111-crystal
    th = braggd(G,wo)
    xxx = vrot(np.array([0.0,-1.0, 1.0]),np.array([-2.0,1.0,1.0]),th)
    yyy = vrot(np.array([-2.0,1.0, 1.0]),np.array([-2.0,1.0,1.0]),th)
    zzz = vrot(G,np.array([-2.0,1.0,1.0]),th)
    #xxx = vrot(np.array([2.0,-1.0,-1.0]),np.array([0.0,-1.0,1.0]),th)
    #yyy = vrot(np.array([0.0,-1.0, 1.0]),np.array([0.0,-1.0,1.0]),th)
    #zzz = vrot(G,np.array([0.0,-1.0,1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # lab-to-sample reference system transformation matrix U for Si220-crystal
    #th = braggd(G,wo)
    #xxx = vrot(np.array([1.0,-1.0,-0.0]),np.array([0.0,0.0,1.0]),th)
    #yyy = vrot(np.array([0.0,0.0, 1.0]),np.array([0.0,0.0,1.0]),th)
    #zzz = vrot(G,np.array([0.0,0.0,1.0]),th)
    #U = np.zeros((3,3))
    #U[:,0] = xxx/np.linalg.norm(xxx)
    #U[:,1] = yyy/np.linalg.norm(yyy)
    #U[:,2] = zzz/np.linalg.norm(zzz)

    # lab-to-sample reference system transformation matrix U for Si1-11-crystal
    #th = braggd(G,wo)
    #xxx = vrot(np.array([2.0,1.0,-1.0]),np.array([0.0,1.0,1.0]),th)
    #yyy = vrot(np.array([0.0,1.0, 1.0]),np.array([0.0,1.0,1.0]),th)
    #zzz = vrot(G,np.array([0.0,1.0,1.0]),th)
    #U = np.zeros((3,3))
    #U[:,0] = xxx/np.linalg.norm(xxx)
    #U[:,1] = yyy/np.linalg.norm(yyy)
    #U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # desired momentum in the laboratory reference system before any rotation is applied
    v_c = np.dot(B,Q)
    Q_lab = np.linalg.lstsq(U,v_c)[0]
    
    #$[angles,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) UBfind(x, G, Q_lab), [0 45 0]);
    lab_angles = optimize.fsolve(cixsUBfind, [20.5, 15.0, 5.0], args=(G,Q_lab,wi,wo,lambdai,lambdao))

    tthv = lab_angles[1]
    tthh = lab_angles[0]
    psi  = lab_angles[2]
    return tthv, tthh, psi

def cixsUBgetAngles_secondo(Q):
    G = np.array([-2.,-2.,0.0])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    lambdao = hc/wo
    wi = wo
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U for Si220-crystal
    th = braggd(G,wo)
    xxx = vrot(np.array([1.0,-1.0,0.0]),np.array([0.0,0.0,1.0]),th)
    yyy = vrot(np.array([0.0,0.0, 1.0]),np.array([0.0,0.0,1.0]),th)
    zzz = vrot(G,np.array([0.0,0.0,1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # desired momentum in the laboratory reference system before any rotation is applied
    v_c = np.dot(B,Q)
    Q_lab = np.linalg.lstsq(U,v_c)[0]
    
    #$[angles,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) UBfind(x, G, Q_lab), [0 45 0]);
    lab_angles = optimize.fsolve(cixsUBfind, [25.5, 0.0, 0.0], args=(G,Q_lab,wi,wo,lambdai,lambdao), xtol=1.49012e-12,maxfev=1000000)

    tthv = lab_angles[1]
    tthh = lab_angles[0]
    psi  = lab_angles[2]
    #if psi <= -360.0:
    #    psi += 360.0
    #if psi >= 360.0:
    #    psi -= 360.0

    return tthv, tthh, psi

def cixsUBgetQ_secondo(tthv, tthh, psi):
    G = np.array([-2.,-2.,0.0])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    lambdao = hc/wo
    wi = wo
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U
    th = braggd(G,wo)
    xxx = vrot(np.array([1.0,-1.0, 0.0]),np.array([0.0,0.0,1.0]),th)
    yyy = vrot(np.array([0.0, 0.0, 1.0]),np.array([0.0,0.0,1.0]),th)
    zzz = vrot(G,np.array([0.0,0.0,1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # axis of rotation of psi
    v = np.array([-np.sin(np.radians(th)), 0.0, np.cos(np.radians(th))])
    Ki_test = 2.0*np.pi/lambdai*X
    Ko_test = 2.0*np.pi/lambdao*vrot(vrot(X,Y,-tthv) ,Z, tthh)
    Q_test = np.dot(np.linalg.lstsq(B,U)[0],vrot(Ki_test-Ko_test,v,-psi))
    return Q_test


def cixsUBgetAngles_terzo(Q):
    G = np.array([-1.0,-1.0,-1.0])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    lambdao = hc/wo
    wi = wo
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U for Si220-crystal
    th = braggd(G,wo)
    #xxx = vrot(np.array([0.0,-1.0,1.0]),np.array([-2.0,1.0,1.0]),th)
    #yyy = vrot(np.array([-2.0,1.0,1.0]),np.array([-2.0,1.0,1.0]),th)
    #zzz = vrot(G,np.array([-2.0,1.0,1.0]),th)
    xxx = vrot(np.array([0.0,1.0,-1.0]),np.array([2.0,-1.0,-1.0]),th)
    yyy = vrot(np.array([2.0,-1.0,-1.0]),np.array([2.0,-1.0,-1.0]),th)
    zzz = vrot(G,np.array([2.0,-1.0,-1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # desired momentum in the laboratory reference system before any rotation is applied
    v_c = np.dot(B,Q)
    Q_lab = np.linalg.lstsq(U,v_c)[0]
    
    #$[angles,FVAL,EXITFLAG,OUTPUT] = fsolve(@(x) UBfind(x, G, Q_lab), [0 45 0]);
    lab_angles = optimize.fsolve(cixsUBfind, [55., 20.0, 0.0], args=(G,Q_lab,wi,wo,lambdai,lambdao), xtol=1.49012e-12,maxfev=1000000)

    tthv = lab_angles[1]
    tthh = lab_angles[0]
    psi  = lab_angles[2]
    #if psi <= -360.0:
    #    psi += 360.0
    #if psi >= 360.0:
    #    psi -= 360.0

    return tthv, tthh, psi

def cixsUBgetQ_terzo(tthv, tthh, psi):
    G = np.array([-1.0,-1.0,-1.0])
    # incoming/outgoing energy/wavelength
    hc = 12.3984191
    bragg_ang = 86.5
    wo = energy(dspace([4., 4., 4.]),bragg_ang)
    lambdao = hc/wo
    wi = wo
    lambdai = hc/wi

    # lattice parameters
    lattice = np.array([5.43095, 5.43095, 5.43095])
    angles  = np.radians(np.array([90.0, 90.0, 90.0])) # in radians !!!
    a = np.array([lattice[0], 0, 0])
    b = np.array([lattice[0]*np.cos(angles[2]), lattice[1]*np.sin(angles[2]), 0])
    c = np.array([lattice[2]*np.cos(angles[1]), lattice[2]*(-np.cos(angles[1])*np.arctan(angles[2])+np.cos(angles[0])*(1.0/np.sin(angles[2]))), lattice[2]/np.sqrt(2.0)*np.sqrt((1.0/np.sin(angles[2]))*((4.0*np.cos(angles[0])*np.cos(angles[1])*np.arctan(angles[2])-(1.0 + np.cos(2.0*angles[0])+np.cos(2.0*angles[1])+np.cos(2.0*angles[2]))*(1.0/np.sin(angles[2])))))])

    # lab-to-sample reference system transformation matrix U
    th = braggd(G,wo)
    #xxx = vrot(np.array([0.0,-1.0,1.0]),np.array([-2.0,1.0,1.0]),th)
    #yyy = vrot(np.array([-2.0,1.0,1.0]),np.array([-2.0,1.0,1.0]),th)
    #zzz = vrot(G,np.array([-2.0,1.0,1.0]),th)
    xxx = vrot(np.array([0.0,1.0,-1.0]),np.array([2.0,-1.0,-1.0]),th)
    yyy = vrot(np.array([2.0,-1.0,-1.0]),np.array([2.0,-1.0,-1.0]),th)
    zzz = vrot(G,np.array([2.0,-1.0,-1.0]),th)
    U = np.zeros((3,3))
    U[:,0] = xxx/np.linalg.norm(xxx)
    U[:,1] = yyy/np.linalg.norm(yyy)
    U[:,2] = zzz/np.linalg.norm(zzz)

    # reciprocal lattice to absolute units transformation matrix
    a_star = 2.0*np.pi*np.cross(b,c)/np.dot(a,np.cross(b,c))
    b_star = 2.0*np.pi*np.cross(c,a)/np.dot(a,np.cross(b,c))
    c_star = 2.0*np.pi*np.cross(a,b)/np.dot(a,np.cross(b,c))
    angles_star = np.array([np.arccos(np.dot(b_star,c_star)/np.linalg.norm(b_star)/np.linalg.norm(c_star)), np.arccos(np.dot(c_star,a_star)/np.linalg.norm(c_star)/np.linalg.norm(a_star)), np.arccos(np.dot(a_star,b_star)/np.linalg.norm(a_star)/np.linalg.norm(b_star))])
    B = np.zeros((3,3))
    B[:,0] = np.array([np.linalg.norm(a_star), np.linalg.norm(b_star)*np.cos(angles_star[2]), np.linalg.norm(c_star)*np.cos(angles_star[1])])
    B[:,1] = np.array([0.0, np.linalg.norm(b_star)*np.sin(angles_star[2]), -np.linalg.norm(c_star)*np.sin(angles_star[1])*np.cos(angles[0])])
    B[:,2] = np.array([0.0, 0.0, 2.0*np.pi/np.linalg.norm(c)])

    # laboratory reference frame
    X = np.array([1.0, 0.0, 0.0])
    Y = np.array([0.0, 1.0, 0.0])
    Z = np.array([0.0, 0.0, 1.0])

    # axis of rotation of psi
    v = np.array([-np.sin(np.radians(th)), 0.0, np.cos(np.radians(th))])
    Ki_test = 2.0*np.pi/lambdai*X
    Ko_test = 2.0*np.pi/lambdao*vrot(vrot(X,Y,-tthv) ,Z, tthh)
    Q_test = np.dot(np.linalg.lstsq(B,U)[0],vrot(Ki_test-Ko_test,v,-psi))
    return Q_test

def cixsUBfind(x,G,Q_sample,wi,wo,lambdai,lambdao):
    """ **cixsUBfind**
    """    
    tthh = x[0]
    tthv = x[1]
    psi  = x[2]
    X = np.array([1, 0, 0])
    Y = np.array([0, 1, 0])
    Z = np.array([0, 0, 1])
    Ki = 2.0*np.pi/lambdai*X
    Ko = 2.0*np.pi/lambdao* vrot(vrot(X,Y,-tthv ),Z,tthh)
    Q = Ki-Ko
    th = braggd(G,wo)
    v  = np.array([-np.sin(np.radians(th)), 0.0, np.cos(np.radians(th))])
    y = Q - vrot(Q_sample, v, psi)
    tthh = y[0]
    tthv = y[1]
    psi  = y[2]
    return tthh, tthv, psi

def cixs_primo(tthv,tthh,psi,anal_braggd=86.5):
    """ **cixs_primo**
    """
    import copy
    lattice_a = dspace([1., 0., 0.]) # Si lattice constant
    # crystal vectors
    crystVec1 = np.array([-1.,-1.,-1.])/np.linalg.norm(np.array([-1.,-1.,-1.])) # "z-axis"
    crystVec2 = np.array([ 0.,-1., 1.])/np.linalg.norm(np.array([ 0.,-1., 1.])) # "x-axis"
    crystVec3 = np.array([-2., 1., 1.])/np.linalg.norm(np.array([-2., 1., 1.])) # "y-axis"
    # rotate x- and y-vectors about G by the miscut of PRIMO
    crystVec2 = vrot(crystVec2,crystVec1,-39.8)
    crystVec3 = vrot(crystVec3,crystVec1,-39.8)
    # calculate energies and wavelengths
    hc      = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
    E_out   = energy(dspace(np.array([4., 4., 4.])),anal_braggd)
    lam_out = hc/E_out
    E_in    = E_out #+0.02; % if want to be precise, E=Eout-20 eV @ plasmon peak
    lam_in  = hc/E_in
    # initially k0 is along crystVec2,
    # then rotate k0 about crystVec3 by the Bragg angle
    k0 = vrot(crystVec2,crystVec3,braggd(np.array([1., 1., 1.]),E_in))
    k0 = k0/np.linalg.norm(k0)*2.0*np.pi/lam_in
    # define lab coordinates
    hutch_x = copy.deepcopy(k0) # k0 is along the beam
    hutch_y = copy.deepcopy(crystVec2) # perpendicular to beam/untouched so far
    hutch_z = vrot(crystVec1,crystVec3,braggd(np.array([1., 1., 1.]),E_in)) # toward hutch ceiling (if k0 rotates, z has to rotate with it)
    # rotate the crystal abouts its G vector
    k0      = vrot(k0,crystVec1,psi)
    hutch_x = copy.deepcopy(k0) # hutch_x is always along k0
    hutch_y = vrot(hutch_y,crystVec1,psi) # perpendicular to beam
    hutch_z = vrot(hutch_z,crystVec1,psi) # toward hutch ceiling
    # calculate kh using G-vector
    kh = k0 + np.array([-1.,-1.,-1.])/np.linalg.norm(np.array([-1.,-1.,-1.]))
    # rotate vertical
    kprime = vrot(k0,hutch_y,-tthv) # we can rotate vertical tth from 0 to 90 (eta from 0 to 90)
    kprime = vrot(kprime,hutch_z,tthh) # we can rotate horizontal tth from 0 to 90
    kprime = kprime/np.linalg.norm(kprime)*2.0*np.pi/lam_out
    # calculate momentum transfer
    qh = kh-kprime
    q0 = k0-kprime
    return q0, qh, kprime #hutch_x, hutch_y, hutch_z

def cixs_secondo(tthv,tthh,psi,anal_braggd=86.5):
    """ **cixs_secondo**
    """
    import copy
    lattice_a = dspace([1., 0., 0.]) # Si lattice constant
    # crystal vectors
    crystVec1 = np.array([-2.,-2., 0.])/np.linalg.norm(np.array([-2.,-2., 0.])) # "z-axis"
    crystVec2 = np.array([ 1.,-1., 0.])/np.linalg.norm(np.array([ 1.,-1., 0.])) # "x-axis"
    crystVec3 = np.array([ 0., 0., 1.])/np.linalg.norm(np.array([ 0., 0., 1.])) # "y-axis"
    # rotate x- and y-vectors about G by the miscut of PRIMO
    crystVec2 = vrot(crystVec2,crystVec1,0.0)
    crystVec3 = vrot(crystVec3,crystVec1,0.0)
    # calculate energies and wavelengths
    hc      = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
    E_out   = energy(dspace(np.array([4., 4., 4.])),anal_braggd)
    lam_out = hc/E_out
    E_in    = E_out #+0.02; % if want to be precise, E=Eout-20 eV @ plasmon peak
    lam_in  = hc/E_in
    # initially k0 is along crystVec2,
    # then rotate k0 about crystVec3 by the Bragg angle
    k0 = vrot(crystVec2,crystVec3,braggd(np.array([1., 1., 1.]),E_in))
    k0 = k0/np.linalg.norm(k0)*2.0*np.pi/lam_in
    # define lab coordinates
    hutch_x = copy.deepcopy(k0) # k0 is along the beam
    hutch_y = copy.deepcopy(crystVec2) # perpendicular to beam/untouched so far
    hutch_z = vrot(crystVec1,crystVec3,braggd(np.array([2., 2., 0.]),E_in)) # toward hutch ceiling (if k0 rotates, z has to rotate with it)
    # rotate the crystal abouts its G vector
    k0      = vrot(k0,crystVec1,psi)
    hutch_x = copy.deepcopy(k0) # hutch_x is always along k0
    hutch_y = vrot(hutch_y,crystVec1,psi) # perpendicular to beam
    hutch_z = vrot(hutch_z,crystVec1,psi) # toward hutch ceiling
    # calculate kh using G-vector
    kh = k0 + np.array([-2.,-2.,0.])/np.linalg.norm(np.array([-2.,-2.,0.]))
    # rotate vertical
    kprime = vrot(k0,hutch_y,-tthv) # we can rotate vertical tth from 0 to 90 (eta from 0 to 90)
    kprime = vrot(kprime,hutch_z,tthh) # we can rotate horizontal tth from 0 to 90
    kprime = kprime/np.linalg.norm(kprime)*2.0*np.pi/lam_out
    # calculate momentum transfer
    qh = kh-kprime
    q0 = k0-kprime
    return q0, qh, hutch_x, hutch_y, hutch_z

def cixs_terzo(tthv,tthh,psi,anal_braggd=86.5):
    """ **cixs_terzo**
    """
    hc = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
    zz = np.array([-1., -1., -1.])
    G  = 2.0*np.pi*zz/dspace(np.array([1., 0., 0.]))

    xx = vrot(np.array([0., 1., -1.,]),np.array([-1., -1., -1.]),90-81.1)
    xx = vrot(xx,G,psi)
    yy = vrot(xx,zz,90.0)

    a = dspace(np.array([1., 0., 0.]))
    Eout = energy(dspace(np.array([4., 4., 4.])),anal_braggd)
    lambdaout = hc/Eout
    E = Eout #+0.02;
    lambdain = hc/E
    k0 = vrot(xx,yy,braggd(zz,E))
    k0 = k0/np.linalg.norm(k0)*2.0*np.pi/lambdain
    nn = vrot(zz,yy,braggd(zz,E)) # nn is our spectrometer (hutch) vertical coordinate
    kh = k0 + G
    kprime = vrot(k0,yy,-tthv) # we can rotate vertical tth from 0 to 90 (eta from 0 to 90)
    kprime = kprime/np.linalg.norm(kprime)*2.0*np.pi/lambdaout

    kprime = vrot(kprime,nn,tthh) # we can rotate horizontal tth from 0
    q0 = k0-kprime
    qh = kh-kprime
    return q0, qh

# def constrained_nnmf(A,W_ini,H_ini,W_up,H_up,max_iter=10000,verbose=False):
#     """ **constrained_nnmf**
#     Approximate non-negative matrix factorization with constrains.
    
#     function [W H]=johannes_nnmf_ALS(A,W_ini,H_ini,W_up,H_up)
#     % *****************************************************************
#     % *****************************************************************
#     % ** [W H]=johannes_nnmf(A,W_ini,H_ini,W_up,H_up)   **
#     % ** performs A=WH approximate matrix factorization,             **
#     % ** where A(n*m), W(n*k), and H(k*m) are non-negative matrices, **
#     % ** and k<min(n,m). Masking arrays W_up(n*k), H_up(k*m) = 0,1   **
#     % ** control elements of W and H to be updated (1) or not (0).   **
#     % ** This fact can be used to set constraints.                   **
#     % **                                                             **
#     % **         Johannes Niskanen 13.10.2015                        **
#     % **                                                             **
#     % *****************************************************************
#     % *****************************************************************
#     by Johannes Niskanen
#     """
#     # initialize matrices
#     H = H_ini
#     W = W_ini
    
#     # initial cost
#     J = np.sum(np.sum(0.5 * (A-np.dot(W,H))*(A-np.dot(W,H))))
#     print('Initial cost J = %1.4f at step 0'%J)
#     dJ = -0.1

#     sind = 0
#     while sind <= max_iter:
#         sind += 1
#         # check singularity
#         if np.isnan(np.linalg.det(np.dot(H,H.T))) or np.abs( np.linalg.det(np.dot(H,H.T))) < 1.0e-12:
#             print('H is singular, will break here.')
#             return

#         # solve W from (H*H')*W'=H*A'
#         W = np.linalg.lstsq( np.dot(H,H.T),np.dot(H,A.T) )[0].T

#         # make W nonnegative
#         inds = W < 0.0
#         W[inds] = 0.0

#         # restore fixed components
#         inds = W_up==0.0
#         W[inds] = W_ini[inds]

#         # check singularity
#         if np.isnan( np.linalg.det(np.dot(W.T,W)) ) or np.abs( np.linalg.det(np.dot(W.T,W)) ) < 1.0e-12:
#             W = np.zeros(np.shape(W))
#             H = np.zeros(np.shape(H))
#             return

#         # solve H from: (W'*W)*H=W'*A
#         H = np.linalg.lstsq( np.dot(W.T,W),np.dot(W.T,A) )[0]

#         # make H non-negative
#         inds = H < 0.0
#         H[inds] = 0.0

#         # restore fixed components
#         inds = H_up == 0.0
#         H[inds] = H_ini[inds]

#         # formalize spectra and coefficients
#         W = W/(np.dot(np.ones((np.shape(W)[0],1)),np.sum(W,axis=0).reshape(1,len(np.sum(W,axis=0))) ))
#         H = H/(np.dot(np.ones((np.shape(H)[0],1)),np.sum(H,axis=0).reshape(1,len(np.sum(H,axis=0))) ))

#         # print some progression
#         if sind % 100 == 0 and verbose:
#             Jnew = np.sum(np.sum(0.5 * (A-np.dot(W,H))*(A-np.dot(W,H))))
#             dJ   = Jnew-J
#             J    = Jnew
#             print('Iteration %1d J = %1.4f') %(sind,J)
#             print('dJ = %5.3f') % dJ
#             print('Fnorm = %5.3f') % np.mean(np.sum(W))
#             print('Cnorm = %5.3f') % np.mean(np.sum(H))

#     return W, H




def mat2con(W,H,W_up,H_up):
    x = W[W_up == 1]
    x = np.append(x, H[H_up == 1])
    return x

def con2mat(x,W,H,W_up,H_up):
    W[W_up == 1] = x[0:len(W[W_up == 1])]
    H[H_up == 1] = x[len(W[W_up == 1]):len(W[W_up == 1])+len(H[H_up == 1])]
    return W, H

def mat2vec(F,C,F_up,C_up,n,k,m):
    idxs=np.where(F_up == 1)
    nF=len(idxs[0])
    if nF>0:
        x = F[idxs]
    else:
        x = np.array([])

    idxs=np.where(C_up == 1)
    nC=len(idxs[0])
    if nC>0:
        x = np.hstack([x, C[idxs]])

    return x


def vec2mat(x,F,C,F_up,C_up,n,k,m):
    idxs=np.where(F_up == 1)
    nF=len(idxs[0])
    if idxs:
        F[idxs] = x[:nF]

    idxs=np.where(C_up == 1)
    nC=len(idxs[0])
    if idxs:
      C[idxs] = x[nF:]

    F=F.reshape(n,k)
    C=C.reshape(k,m)
    return F, C 



def NNMFcost_old(x,A,W,H,W_up,H_up):
    """ **NNMFcost**
    Returns cost and gradient for NNMF with constraints.
    """
    # calculate W, H
    W, H = con2mat(x,W,H,W_up,H_up)
    # calculate cost and gradient
    J = np.sum(np.sum(0.5*(A-np.dot(W,H))*(A-np.dot(W,H))))
    gradW = -(np.dot((A-np.dot(W,H)),H.T))
    gradH = -(np.dot((A-np.dot(W,H)).T,W)).T
    # return constraint only for updates
    xgrad = mat2con(gradW,gradH,W_up,H_up)
    return J, xgrad



def NNMFcost(x,A,F,C,F_up,C_up,n,k,m):
    """ **NNMFcost**
    Returns cost and gradient for NNMF with constraints.
    """
    # calculate W, H
    F, C = vec2mat(x,F,C,F_up,C_up,n,k,m)
    # calculate cost and gradient
    J = np.sum(np.sum(0.5*(A-np.dot(F,C))*(A-np.dot(F,C))))
    gradF = -(np.dot((A-np.dot(F,C)),C.T))
    gradC = -(np.dot((A-np.dot(F,C)).T,F)).T
    # return gradient only for updates
    gradF=gradF.reshape(1,n*k)[0]
    gradC=gradC.reshape(1,k*m)[0]
    xgrad = mat2vec(gradF,gradC,F_up,C_up,n,k,m)
    return J


def bootstrapCNNMF(A,F_ini, C_ini, F_up, C_up, Niter):
    """ **bootstrapCNNMF**
    Constrained non-negative matrix factorization with bootstrapping
    for error estimates.
    """
    n,m = A.shape
    k   = np.shape(F_ini)[1]
    print( 'NNMF problem of dimension: ' + str(n) + 'x' + str(k) + 'x' + str(m) )
    F1s = []
    C1s = []
    for ii in range(Niter):
        A1 = copy.deepcopy(A)
        # add random noise
        #A1 += np.random.random((n,m))*Aerr
        F1 = F_ini.reshape(1,n*k)[0]
        C1 = C_ini.reshape(1,k*m)[0]
        F1up=F_up.reshape(1,n*k)[0]
        C1up=C_up.reshape(1,k*m)[0]
        #print( len(np.where(F1up==1)[0]))
        #print( len(np.where(C1up==1)[0]))
        # get starting values
        x0 = mat2vec(F1,C1,F1up,C1up,n,k,m)
        # set limits
        bnds = [(0.0,1.0) for ii in x0]
        # minimize cost
        res=optimize.minimize(NNMFcost,x0,args=(A1,F1,C1,F1up,C1up,n,k,m),
                         bounds=bnds,
                         method='SLSQP',
                         options={'maxiter': 100000, 'disp': True},
                         tol=1e-10)
        #translate vector of parameters bacj to matrices
        Fbs1, Cbs1 = vec2mat(res.x,F1,C1,F1up,C1up,n,k,m)
        # Normalize, translated to python by Risto & Johannes
        Fbs1=Fbs1/(np.dot(np.ones((n,1)),np.sum(Fbs1,axis=0).reshape(1,k)))
        Cbs1=Cbs1*(np.dot(np.sum(Fbs1,axis=0).reshape(k,1), np.ones((1,m))))
        # store meaningful data
        F1s.append(Fbs1.copy())
        C1s.append(Cbs1.copy())

    if Niter>1:
        # stantard deviation
        Cerr=np.squeeze(np.std(np.array(C1s),axis=0))
        Ferr=np.squeeze(np.std(np.array(F1s),axis=0))
        # average
        C=np.squeeze(np.mean(np.array(C1s),axis=0))
        F=np.squeeze(np.mean(np.array(F1s),axis=0))
    else:
        # only one of each matrix
        F=np.array(F1s[0])
        C=np.array(C1s[0])
        Ferr=np.zeros(np.shape(F))
        Cerr=np.zeros(np.shape(C))
    return F, C, Ferr, Cerr 


def bootstrapCNNMF_old(A,k,Aerr, F_ini, C_ini, F_up, C_up, Niter=100):
    """ **bootstrapCNNMF**
    Constrained non-negative matrix factorization with bootstrapping
    for error estimates.
    """
    n,m = A.shape
    import copy
    F1s = np.zeros((Niter,n,k))
    C1s = np.zeros((Niter,C_ini.shape[0],C_ini.shape[1]))
    for ii in range(Niter):
        A1 = copy.deepcopy(A)
        # add random noise
        A1 += np.random.random((n,m))*Aerr
        F1 = F_ini*(1.0-F_up) + F_up*np.random.random((n,k))
        C1 = C_ini*(1.0-C_up) + C_up*np.random.random((k,m))
        F1[F1<0.0]=0.0
        C1[C1<0.0]=0.0
        # minimize with trust-region-algorithm
        # get starting values
        x0 = mat2con(F1,C1,F_up,C_up)
        cons = ({'args': (A1,F1,C1,F_up,C_up)})
        bnds = [(0.0,1.0) for ii in x0]
        costfun = lambda x:NNMFcost(x,A1,F1,C1,F_up,C_up) #[0]
        #gradfun = lambda x:NNMFcost(x,A1,F1,C1,F_up,C_up)[1]
        x=minimize(NNMFcost_old,x0,args=(A1,F1,C1,F_up,C_up), method='Newton-CG', tol=1e-5, jac=True, bounds=bnds,options={'maxiter' : 1e6, 'disp': True} ).x
        Fbs1, Cbs1 = con2mat(x,F1,C1,F_up,C_up)
        # store meaningful data 
        print( Fbs1.shape)
        print( Cbs1.shape)
        F1s[ii,:,:] = Fbs1/(np.dot(np.ones((np.shape(Fbs1)[0],1)), np.sum(Fbs1,axis=0).reshape(1,len(np.sum(Fbs1,axis=0))) ))
        print( )
        C1s[ii,:,:] = Cbs1 * ( np.dot(np.sum(Fbs1,axis=0).reshape(k,1), np.ones(( 1, A1.shape[1] )) ) )
    # do RMS
    print( F1s.shape, C1s.shape)
    Cerr=np.squeeze(np.std(C1s,axis=0))
    Ferr=np.squeeze(np.std(F1s,axis=0))
    # average
    C=np.squeeze(np.mean(C1s,axis=0))
    F=np.squeeze(np.mean(F1s,axis=0))
    return F, C, Ferr, Cerr









def constrained_svd(M,U_ini,S_ini,VT_ini,U_up,max_iter=10000,verbose=False):
    """ **constrained_nnmf**
    Approximate singular value decomposition with constraints.
    
    function [U, S, V] = constrained_svd(M,U_ini,S_ini,V_ini,U_up,max_iter=10000,verbose=False)
    """
    # initialize matrices
    # M = [n x m]
    U  = U_ini # [n x n] (unitary)
    S  = S_ini # [n x m] (diagonal matrix)
    VT = VT_ini # [m x m] (unitary)

    n,m = np.shape(M)

    # initial cost
    J = np.sum(np.sum(0.5 * (M-np.dot(np.dot(U,S),VT))*(M-np.dot(np.dot(U,S),VT))))
    print('Initial cost J = %1.4f at step 0') % J
    dJ = -0.1

    sind = 0
    while sind <= max_iter:
        sind += 1

        # solve S from: U*S = M*(VT)^-1
        S = np.linalg.lstsq( U, np.dot(M, np.linalg.pinv(VT)))[0]
        # make S diagonal
        for ii in range(S.shape[0]):
            for jj in range(S.shape[1]):
                if ii != jj:
                    S[ii,jj] = 0.0

        # solve VT from: U*S*V=M
        VT = np.linalg.lstsq( np.dot(U,S),M )[0]

        # solve U from: VT.T*S.T*U.T = M.T
        U = np.linalg.lstsq( np.dot(VT.T,S.T) , M.T  )[0].T

        # restore fixed components
        inds = U_up==0.0
        U[inds] = U_ini[inds]

        # formalize spectra and coefficients
        U  = U/(np.dot(np.ones((np.shape(U)[0],1)),np.sum(U,axis=0).reshape(1,len(np.sum(U,axis=0))) ))
        VT = VT/(np.dot(np.ones((np.shape(VT)[0],1)),np.sum(VT,axis=0).reshape(1,len(np.sum(VT,axis=0))) ))


        # print some progression
        if sind % 100 == 0 and verbose:
            Jnew = np.sum(np.sum(0.5 * (M-np.dot(np.dot(U,S),VT))*(M-np.dot(np.dot(U,S),VT))))
            dJ   = Jnew-J
            J    = Jnew
            print('Iteration %1d J = %1.4f') %(sind,J)
            print('dJ = %5.3f') % dJ

    return U, S, VT

def unconstrained_mf(A,numComp=3, maxIter=1000, tol=1.0e-8):
    """ **unconstrained_mf**
    Returns main components from an off-diagonal Matrix (energy-loss x angular-departure),
    using the power method iteratively on the different main components.
    """
    # initialize random coefficient matrix
    coeff  = np.random.random((A.shape[1],numComp))
    W      = np.random.random((numComp,A.shape[0]))

    # normalize W
    for ii in range(numComp):
        W[ii,:] /= np.linalg.norm(W[ii,:])

    ind  = 0
    err  = 1.0e8

    # start looping:
    while ind <= maxIter or dJ <= tol:
        # update coefficient matrix
        abc   = np.linalg.lstsq( W.T,A)[0].T
        coeff = np.copy(abc)
        for comp in range(numComp):
            # updatea coefficients and
            # set one of the coefficient vectors to zero
            coeff[:,comp] = np.zeros_like(abc[:,comp])
            # calculate error matrix
            errM = A - np.dot(coeff,W).T
            # initialize power method
            V =  np.random.random((len(W[comp,:]),1))
            V /= np.linalg.norm(V)
            for jj in range(1000):
                vnew = np.dot(errM, errM.T).dot(V)
                vnew /= np.linalg.norm(vnew)
                V = vnew
            V /= np.linalg.norm(V)
            W[comp,:] = V.reshape(W[comp,:].shape)
            # set the zeroed coefficients back to orig
            coeff[:,comp] = abc[:,comp]
        # calculate error
        newerr = np.linalg.norm(A - np.dot(coeff,W).T)
        dJ = err - newerr
        err = newerr
        ind += 1
    return W, coeff, err

def constrained_mf(A, W_ini, W_up, coeff_ini, coeff_up,  maxIter=1000, tol=1.0e-8, maxIter_power=1000):
    """ **cfactorizeOffDiaMatrix**
    constrained version of factorizeOffDiaMatrix
    Returns main components from an off-diagonal Matrix (energy-loss x angular-departure).
    """
    numComp = coeff_ini.shape[1]
    # initialize random coefficient matrix
    coeff = np.copy(coeff_ini)
    W     = np.copy(W_ini)
    # normalize W
    for ii in range(numComp):
        W[:,ii] /= np.linalg.norm(W[:,ii])
    # looping index
    ind  = 0
    err  = 1.0e8
    # find columns to be updated
    W_up_cols     = []
    coeff_up_cols = []
    for ii in range(numComp):
        if np.all(W_up[:,ii] == 1):
            W_up_cols.append(ii) 
        if np.all(coeff_up[:,ii] == 1):
            coeff_up_cols.append(ii)
    # start looping:
    while ind <= maxIter:
        # update coefficient matrix where desired
        abc   = np.linalg.lstsq( W,A)[0].T
        coeff = np.copy(abc)
        coeff[:,coeff_up_cols] = abc[:,coeff_up_cols]
        for col in W_up_cols:
            # set one of the coefficient vectors to zero
            coeff[:,col] = np.zeros_like(coeff[:,col])
            # calculate error matrix
            errM = A - np.dot(coeff,W.T).T
            # initialize power method
            V =  np.random.random((len(W[:,col]),1))
            V /= np.linalg.norm(V)
            for jj in range(maxIter_power):
                vnew = np.dot(errM, errM.T).dot(V)
                vnew /= np.linalg.norm(vnew)
                V = vnew
            V /= np.linalg.norm(V)
            W[:,col] = V.reshape(W[:,col].shape)
            # set the zeroed coefficients back to orig
            coeff[:,col] = abc[:,col]
        # calculate error
        newerr = np.linalg.norm(A - np.dot(coeff,W.T).T)
        dJ = err - newerr
        err = newerr
        ind += 1
    return W, coeff, err













def readbiggsdata(filename,element):
    """
    Reads Hartree-Fock Profile of element 'element' from values tabulated 
    by Biggs et al. (Atomic Data and Nuclear Data Tables 16, 201-309 (1975))
    as provided by the DABAX library (http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/ComptonProfiles.dat).
    input:
    filename = path to the ComptonProfiles.dat file (the file should be distributed with this package)
    element  = string of element name
    returns:

      * data     = the data for the according element as in the file:
          * #UD  Columns: 
          * #UD  col1: pz in atomic units 
          * #UD  col2: Total compton profile (sum over the atomic electrons
          * #UD  col3,...coln: Compton profile for the individual sub-shells

      * occupation = occupation number of the according shells
      * bindingen  = binding energies of the accorting shells
      * colnames   = strings of column names as used in the file
    """
    elementid = '#S'
    sizeid    = '#N'
    occid     = '#UOCCUP'
    bindingid = '#UBIND'
    colnameid = '#L'
    data = []
    f = open(filename,'r')
    istrue = True
    while istrue:
        line = f.readline()
        if line[0:2] == elementid:
            if line.split()[-1] == element:
                line = f.readline()
                while line[0:2] != elementid:
                    if line[0:2] == sizeid:
                        arraysize = int(line.split()[-1])
                        line = f.readline()
                    if line[0:7] == occid:
                        occupation = line.split()[1:]
                        line = f.readline()
                    if line[0:6] == bindingid:
                        bindingen = line.split()[1:]    
                        line = f.readline()
                    if line[0:2] == colnameid:
                        colnames = line.split()[1:]
                        line = f.readline()
                    if line[0]== ' ':
                        data.append([float(n) for n in line.strip().split()])
                        #data = np.zeros((31,arraysize))
                        line = f.readline()
                break
    length = len(data)
    data = (np.reshape(np.array(data),(length,arraysize)))
    return data, occupation, bindingen, colnames

def makepzprofile(element,filename=os.path.join(data_installation_dir,'ComptonProfiles.dat')):
    """
    constructs compton profiles of element 'element' on pz-scale 
    (-100:100 a.u.) from the Biggs tables provided in 'filename'

    input:
      * element   = element symbol (e.g. 'Si', 'Al', etc.)
      * filename  = path and filename to tabulated profiles

    returns:
      * pzprofile = numpy array of the CP:
        *  1. column: pz-scale
        *  2. ... n. columns: compton profile of nth shell
        * binden     = binding energies of shells
        * occupation = number of electrons in the according shells
    """
    theory,occupation,binden,colnames = readbiggsdata(filename,element)
    # first spline onto a rough grid:
    roughpz = np.logspace(0.01,2,65)-1
    roughtheory      = np.zeros((len(roughpz),len(binden)+2))
    roughtheory[:,0] = roughpz
    for n in range(len(binden)+1):
        intf               = interpolate.pchip(theory[:,0], theory[:,2]) # interpolate.interp1d(theory[:,0],theory[:,n+1])
        roughtheory[:,n+1] = intf(roughpz)
    pzscale   = np.linspace(-100,100,num=4000)
    pzprofile      = np.zeros((len(pzscale),len(binden)+1))
    pzprofile[:,0] = pzscale     
    # mirror, spline onto fine grid
    for n in range(len(binden)):
        intf             = interpolate.splrep(roughtheory[:,0],roughtheory[:,n+2],s=0.000000001,k=2) # skip the column with the total J for now #try interp1d with bounds_error=False and fill_value=0.0
        pzprofile[:,n+1] = interpolate.splev(abs(pzscale),intf,der=0)
    # normalize to one electron, multiply by number of electrons
    for n in range(len(binden)):
        normval = integrate.trapezoid(pzprofile[:,n+1],pzprofile[:,0])
        pzprofile[:,n+1] = pzprofile[:,n+1]/normval*int(occupation[n])
    binden     = [float(en) for en in binden]
    occupation = [float(val) for val in occupation]
    return pzprofile, binden, occupation

def makeprofile(element,filename=os.path.join(data_installation_dir,'ComptonProfiles.dat'),E0=9.69,tth=35.0,correctasym=None):
    """
    takes the profiles from 'makepzprofile()', converts them onto eloss 
    scale and normalizes them to S(q,w) [1/eV]
    input:
    element  = element symbol (e.g. 'Si', 'Al', etc.)
    filename = path and filename to tabulated profiles
    E0       = scattering energy [keV]
    tth      = scattering angle  [deg]
    returns:
    enscale = energy loss scale
    J = total CP
    C = only core contribution to CP
    V = only valence contribution to CP
    q = momentum transfer [a.u.]
    """
    pzprofile,binden,occ = makepzprofile(element,filename)
    # convert to eloss scale
    enscale = ((np.flipud(pz2e1(E0,pzprofile[:,0],tth))-E0)*1e3)
    q = momtrans_au(enscale/1000.0+E0,E0,tth)
    # add asymmetry if needed (2p1/2 and 2p3/2 for Z > 35 (Br))
    asymmetry = np.flipud(HRcorrect(pzprofile,occ,q));  # asymmetry flipped for conversion to e-loss scale (???)
    if correctasym:
        pzprofile[:,1:4] = pzprofile[:,1:4] + asymmetry*correctasym

    # discard profiles below zero
    hfprofile = pzprofile[np.nonzero(enscale.T>=0)[0],:]
    q         = q[np.nonzero(enscale.T>=0)[0]] #q[:,np.nonzero(enscale.T>=0)[0]]
    enscale   = enscale[np.nonzero(enscale.T>=0)[0]] #enscale[:,np.nonzero(enscale.T>=0)[0]]
    hfprofile[:,0] = enscale
    # cut at edges
    for n in range(len(binden)):
        hfprofile[np.where(enscale<binden[n]),n+1] = 0 
    # convert J(pz) to S(q,w) via J(pz)=N_electrons*hartree*q*S(q,w) and
    # normalize using the f-sum rule (sum(S(q,w)*w)=f)
    # convert to a.u.
    hartree = 1.0/constants.physical_constants['electron volt-hartree relationship'][0]
    enscaleh = enscale/hartree # eloss in a.u.
    # normalize to one then multiply by N_el*q**2/2
    for n in range(len(binden)):
        hfprofile[:,n+1] = hfprofile[:,n+1]/(integrate.trapezoid(np.multiply(hfprofile[:,n+1],enscaleh),enscaleh))
        hfprofile[:,n+1] = np.multiply(hfprofile[:,n+1],(q**2.0)/2.0)*occ[n]
    # convert back to [1/eV] and sum up
    # total profile J and valence V (all edges )
    J = np.zeros((len(enscale)))
    V = np.zeros((len(enscale)))
    for n in range(len(binden)):
        if binden[n] < enscale[-1]:
            J += hfprofile[:,n+1]/hartree
            if binden[n] < 10:
                V += hfprofile[:,n+1]/hartree
    C = J - V
    return enscale,J,C,V,q

def makeprofile_comp(formula,filename=os.path.join(data_installation_dir,'ComptonProfiles.dat'),E0=9.69,tth=35,correctasym=None):
    """
    returns the compton profile of a chemical compound with formula 'formula'
    input:
    formula = string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)
    filename = path and filename to tabulated profiles
    E0       = scattering energy [keV]
    tth      = scattering angle  [deg]
    returns:
    eloss = energy loss scale
    J = total CP
    C = only core contribution to CP
    V = only valence contribution to CP
    q = momentum transfer [a.u.]
    """
    elements,stoichiometries = parseformula(formula)
    
    if not np.any(correctasym):
        correctasym = np.zeros(len(elements))
        
    eloss,J,C,V,q = makeprofile(elements[0],filename,E0,tth,correctasym[0])
    J *= stoichiometries[0]
    C *= stoichiometries[0]
    V *= stoichiometries[0]

    for n in range(len(elements[1:])):
        eloss,j,c,v,q = makeprofile(elements[n+1],filename,E0,tth,correctasym[n+1])
        J += j*stoichiometries[n+1]
        C += c*stoichiometries[n+1]
        V += v*stoichiometries[n+1]
    return eloss, J,C,V,q


#os.path.join(data_installation_dir,'data/ComptonProfiles.dat')

def makeprofile_compds(formulas,concentrations=None,filename=os.path.join( data_installation_dir,'ComptonProfiles.dat'),E0=9.69,tth=35.0,correctasym=None):
    """
    returns sum of compton profiles from a lost of chemical compounds weighted by the given concentration
    """
    # if correctasym is not given, no HR correction is applied 
    if not np.any(concentrations):
        concentrations = np.ones(len(formulas))/len(formulas)
    if not np.any(correctasym):
        correctasym = []
        for formula in formulas:
            elements,stoichiometries = parseformula(formula)
            correctasym.append(np.zeros(len(elements)))
    
    eloss,J,C,V,q = makeprofile_comp(formulas[0],filename,E0,tth,correctasym[0])
    if len(formulas)>1:
        J = J*concentrations[0]
        C = C*concentrations[0]
        V = V*concentrations[0]
        for n in range(len(formulas[1:])):
            eloss,j,c,v,q = makeprofile_comp(formulas[n+1],filename,E0,tth,correctasym[n+1])
            J += j*concentrations[n+1]
            C += c*concentrations[n+1]
            V += v*concentrations[n+1]
    return eloss,J,C,V,q

def HRcorrect(pzprofile,occupation,q):
    """ Returns the first order correction to filled 1s, 2s, and 2p Compton profiles.

    Implementation after Holm and Ribberfors (citation ...).

    Args: 
      * pzprofile (np.array): Compton profile (e.g. tabulated from Biggs) to be corrected (2D matrix). 
      * occupation (list): electron configuration.
      * q (float or np.array): momentum transfer in [a.u.].

    Returns:
       asymmetry (np.array):  asymmetries to be added to the raw profiles (normalized to the number of electrons on pz scale)
    """
    # prepare output matrix
    if len(occupation) == 1:
        asymmetry = np.zeros((len(pzprofile[:,0]),1))
    elif len(occupation) == 2:
        asymmetry = np.zeros((len(pzprofile[:,0]),2))
    elif len(occupation) >= 3:
        asymmetry = np.zeros((len(pzprofile[:,0]),3))

    # take care for the cases where 2p levels have spin-orbit split taken into account in the Biggs table
    if len(occupation)>3 and occupation[2]==2 and occupation[3]==4:
        pzprofile[:,3] = pzprofile[:,3] + pzprofile[:,4]
        occupation[2] = 6
    
    # 1s 
    if occupation[0] < 2:
        pass
    else:
        # find gamma1s lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
        fitfct  = lambda a: (np.absolute(np.max(pzprofile[:,1])-np.max(occupation[0]*8.0*a**5.0/3.0/np.pi/(a**2.0+pzprofile[:,0]**2.0)**3.0)))
        res = optimize.leastsq(fitfct,np.sum(occupation))
        gamma1s = res[0][0]
        # calculate j0 and j1
        j0 = occupation[0]*8.0*gamma1s**5.0/3.0/np.pi/((gamma1s**2.0+pzprofile[:,0]**2.0)**3.0)
        j1 = 2.0*gamma1s*np.arctan2(pzprofile[:,0],gamma1s)-3.0/2.0*pzprofile[:,0] 
        j1 = j1/q*j0
        asymmetry[:,0] = j1
    # 2s
    if len(occupation)>1:
        if occupation[1] < 2:
            pass
        else:
            # find gamma2s
            fitfct  = lambda a: (np.absolute(np.max(pzprofile[:,2])-np.max(occupation[1]*((a**4.0-10.0*a**2.0*pzprofile[:,0]**2 + 40.0*pzprofile[:,0]**4.0)*128.0*a**5.0/15.0/np.pi/(a**2.0 + 4.0*pzprofile[:,0]**2.0)**5.0))))
            res = optimize.leastsq(fitfct,np.sum(occupation)*2.0/3.0)
            gamma2s = res[0][0]
            # calculate j0 and j1
            j0 = occupation[1]*(gamma2s**4.0-10.0*gamma2s**2.0*pzprofile[:,0]**2.0+40.0*pzprofile[:,0]**4.0)*128.0*gamma2s**5.0/15.0/np.pi/(gamma2s**2.0 + 4.0*pzprofile[:,0]**2.0)**5.0
            j1 = 2.0*gamma2s*np.arctan2(2.0*pzprofile[:,0],gamma2s)-5.0/4.0*(gamma2s**4.0+48.0*pzprofile[:,0]**4.0)/(gamma2s**4.0-10.0*gamma2s**2.0*pzprofile[:,0]**2.0+40.0*pzprofile[:,0]**4.0)*pzprofile[:,0] 
            j1 = j1/q*j0
            asymmetry[:,1] = j1
    # 2p
    if len(occupation)>2:
        if occupation[2] < 6:
            pass
        else:
            forgamma = 3.0*pzprofile[:,3]/np.trapezoid(pzprofile[:,3],pzprofile[:,0]) # 2p correction is defined for 3 electrons in the 2p shell
            # find gamma2p
            fitfct = lambda a: (np.absolute(np.max(forgamma)-np.max(((a**2.0+20.0*pzprofile[:,0]**2.0)*64.0*a**7.0/5.0/np.pi/(a**2.0+4.0*pzprofile[:,0]**2.0)**5.0))))
            res = optimize.leastsq(fitfct,np.sum(occupation)*1.0/3.0)
            gamma2p = res[0][0]
            # calculate j0 and j1
            j0 = 2.0*(gamma2p**2.0+20.0*pzprofile[:,0]**2.0)*64.0*gamma2p**7.0/5.0/np.pi/(gamma2p**2.0+4.0*pzprofile[:,0]**2.0)**5.0
            j1 = 2.0*gamma2p*np.arctan2(2.0*pzprofile[:,0],gamma2p)-2.0/3.0*pzprofile[:,0]*(10.0*gamma2p**2.0+60.0*pzprofile[:,0]**2.0)/(gamma2p**2.0+20.0*pzprofile[:,0]**2.0)
            j1 = j1/q*j0
            asymmetry[:,2] = j1
    return asymmetry

def parseformula(formula):
    """Parses a chemical sum formula.

    Parses the constituing elements and stoichiometries from a given 
    chemical sum formula.

    Args:
      * formula (string): string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)

    Returns:
      * elements (list): list of strings of constituting elemental symbols.
      * stoichiometries (list): list of according stoichiometries in the same order as 'elements'.
    """
    elements = []
    stoichiometries = []
    splitted = findall(r'([A-Z][a-z]*)(\d*)',formula)
    elements.extend([element[0] for element in splitted])
    stoichiometries.extend([(int(element[1]) if element[1] else 1) for element in splitted])
    return elements,stoichiometries

def element(z):
    """Converts atomic number into string of the element symbol and vice versa.

    Returns atomic number of given element, if z is a string of the 
    element symbol or string of element symbol of given atomic number z.

    Args:
     * z (string or int): string of the element symbol or atomic number. 

    Returns:
       * Z (string or int): string of the element symbol or atomic number.
    """
    zs = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
              'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
              'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo',
              'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba',
              'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
              'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf',
              'Es','Fm','Md','No','Lr','Ku']
    if isinstance(z,str):
        try:
            Z = zs.index(z)+1
        except:
            Z = None
            print( 'Given element ' + z + ' unknown.')
    elif isinstance(z,int):
        if z > 0 and z < 105:
            Z = zs[z-1]
        else:
            print( 'Element Z = '+ str(z) +' unknown.')
    else:
        print( 'type '+ str(type(z)) + 'not supported.'    )
    return Z

#os.path.join(data_installation_dir,'data/logtable.dat')

def myprho(energy,Z,logtablefile=os.path.join(data_installation_dir,'logtable.dat') ):
    """Calculates the photoelectric, elastic, and inelastic absorption of 
    an element Z 

    Calculates the photelectric , elastic, and inelastic absorption of an element Z. 
    Z can be atomic number or element symbol.

    Args:
      * energy (np.array): energy scale in [keV].
      * Z (string or int): atomic number or string of element symbol.

    Returns:
      * murho (np.array): absorption coefficient normalized by the density.
      * rho (float): density in UNITS?
      * m (float): atomic mass in UNITS?

    """
    en = np.array([])
    en = np.append(en,energy) 
    logtable = np.loadtxt(logtablefile)
    # find the right places in logtable
    if not isinstance(Z,int):
        Z = element(Z)
    try:
        ind = list(logtable[:,0]).index(Z)
    except:
        print( 'no such element in logtable.dat')
    c     = np.array(logtable[ind:ind+5,:]) # 5 lines that corresponds to the element
    le    = np.log(en) # logarithm of the energy
    mr    = np.exp(c[1,3]+le*(c[2,3]+le*(c[3,3]+le*c[4,3])))  # extract mu from loglog table
    i     = np.where(en<=c[0,3])
    l     = le[i]
    mr[i] = np.exp(c[1,2]+l*(c[2,2]+l*(c[3,2]+l*c[4,2])))
    i     = np.where(en<c[0,2])
    l     = le[i]
    mr[i] = np.exp(c[1,1]+l*(c[2,1]+l*(c[3,1]+l*c[4,1])))
    i     = np.where(en<c[0,1])
    l     = le[i]
    # mu
    mu      = np.zeros((len(mr),3))
    mu[:,0] = mr
    mu[i,0] = np.exp(c[1,0]+l*(c[2,0]+l*(c[3,0]+l*c[4,0]))) # photoelectric absorption
    mu[:,1] = np.exp(c[1,4]+le*(c[2,4]+le*(c[3,4]+le*c[4,4]))) # elastic absorption
    mu[:,2] = np.exp(c[1,5]+le*(c[2,5]+le*(c[3,5]+le*c[4,5]))) # inelastic abssorption
    #
    m = c[0,4] # atomic mass
    murho = mu*0.602252/m # mu/rho
    rho = c[0,5]
    return murho, rho, m

def mpr(energy,compound):
    """Calculates the photoelectric, elastic, and inelastic absorption of 
    a chemical compound.

    Calculates the photoelectric, elastic, and inelastic absorption of a
    chemical compound.

    Args:
      * energy (np.array): energy scale in [keV].
      * compound (string): chemical sum formula (e.g. 'SiO2')

    Returns:
      * murho (np.array): absorption coefficient normalized by the density.
      * rho (float): density in UNITS?
      * m (float): atomic mass in UNITS?

    """
    en   = np.array([])
    en   = np.append(en,energy) # turn energy into an iterable array
    z,w  = parseformula(compound)
    mr   = np.zeros((len(en),3)) # 1. photoelectric absorption, 2. elastic absorption, 3. inelastic absorption
    rhov = np.zeros((len(z),1))
    mv   = np.zeros((len(z),1))
    for i in range(len(z)):
        tmp,rho,m = myprho(en,z[i])
        m         = m*w[i] # weigh atomic masses by stoichiometry.
        mv[i]     = m
        rhov[i]   = rho
        mr       += tmp*m # sum up individual mu/rho
    mtot = sum(mv)
    mr   = mr/mtot
    mr      = np.sum(mr,1)
    return mr, rhov, mv

def mpr_compds(energy,formulas,concentrations,E0,rho_formu):
    """Calculates the photoelectric, elastic, and inelastic absorption of 
    a mix of compounds.

    Returns the photoelectric absorption for a sum of different chemical 
    compounds.

    Args:
      * energy (np.array): energy scale in [keV].
      * formulas (list of strings): list of chemical sum formulas

    Returns:
      * murho (np.array): absorption coefficient normalized by the density.
      * rho (float): density in UNITS?
      * m (float): atomic mass in UNITS?

    """
    en  = np.array([]) # turn energy into an iterable array
    en  = np.append(en,energy)
    e0  = np.array([])
    e0  = np.append(e0,E0)
    mu_tot_in  = np.zeros((len(en)))
    mu_tot_out = np.zeros((len(e0))) # should also work for series of E0's 
    for n in range(len(formulas)):
        mu_tot_in += mpr(en,formulas[n])[0]*concentrations[n]*rho_formu[n]
        mu_tot_out += mpr(e0,formulas[n])[0]*concentrations[n]*rho_formu[n]
    return mu_tot_in, mu_tot_out

def abscorr2(mu1,mu2,alpha,beta,samthick):
    """Calculates absorption correction for given mu1 and mu2.
    Multiply the measured spectrum with this correction factor.

    This is a translation of Keijo Hamalainen's Matlab function (KH 30.05.96).

    Args:
      * mu1 (np.array): absorption coefficient for the incident energy in [1/cm].
      * mu2 (np.array): absorption coefficient for the scattered energy in [1/cm].
      * alpha (float): incident angle relative to plane normal in [deg].
      * beta (float): exit angle relative to plane normal [deg] (for transmission geometry use beta < 0).
      * samthick (float): sample thickness in [cm].

    Returns:
      * ac (np.array): absorption correction factor. Multiply this with your measured spectrum.

    """
    cosa = math.cos(math.radians(alpha))
    cosb = math.cos(math.radians(beta))
    if beta >= 0: # reflection geometry
        ac =  cosa*(mu1/cosa + mu2/cosb)/(1.0 - np.exp(-mu1*samthick/cosa - mu2*samthick/cosb))
    elif np.absolute(mu1/cosa - mu2/cosb).any() > np.spacing(1): # transmission geometry
        ac = -cosa*(mu1/cosa - mu2/cosb)/(np.exp(-mu1*samthick/cosa) - np.exp(-mu2*samthick/cosb))
    else:
        ac = cosa/(samthick*np.exp(-mu1*samthick/cosa))
    return ac

def absCorrection(mu1,mu2,alpha,beta,samthick,geometry='transmission'):
    """
          **absCorrection**

          Calculates absorption correction for given mu1 and mu2.
          Multiply the measured spectrum with this correction factor.
          This is a translation of Keijo Hamalainen's Matlab function (KH 30.05.96).

          Args
            * mu1 : np.array  Absorption coefficient for the incident energy in [1/cm].
            * mu2 : np.array Absorption coefficient for the scattered energy in [1/cm].
            * alpha : float Incident angle relative to plane normal in [deg].
            * beta : float  Exit angle relative to plane normal [deg].
            * samthick : float  Sample thickness in [cm].
            * geometry : string, optional
              Key word for different sample geometries ('transmission', 'reflection', 'sphere'). 
              If *geometry* is set to 'sphere', no angular dependence is assumed.
      
          Returns
            * ac : np.array
              Absorption correction factor. Multiply this with your measured spectrum.

    """
    cosa = np.cos(math.radians(alpha))
    cosb = np.cos(math.radians(beta))

        # reflection geometry
    if geometry == 'reflection':
        if beta >= 90.0:
            print('WARNING: are you sure about the beta angle?')
        ac =  cosa*(mu1/cosa + mu2/cosb)/(1.0 - np.exp(-mu1*samthick/cosa - mu2*samthick/cosb))

        # transmission geometry
    elif geometry == 'transmission' and np.absolute(mu1/cosa - mu2/cosb).any() > np.spacing(1):
        ac = -cosa*(mu1/cosa - mu2/cosb)/(np.exp(-mu1*samthick/cosa) - np.exp(-mu2*samthick/cosb))
    elif geometry == 'transmission' and np.absolute(mu1/cosa - mu2/cosb).any() <= np.spacing(1):
        ac = cosa/(samthick*np.exp(-mu1*samthick/cosa))

        # spherical sample
    elif geometry == 'sphere':
        ac = (mu1 + mu2)/(1.0 - np.exp(-mu1*samthick -mu2*samthick))
    return ac

def gettransmission(energy,formulas,concentrations,densities,thickness):
    """
    returns the transmission through a sample composed of chemical formulas 
    with certain densities mixed to certain concentrations, and a thickness
    """
    en  = np.array([]) # turn energy into an iterable array
    en  = np.append(en,energy)
    if not isinstance(formulas,list):
        theformulas = []
        theformulas.append(formulas)
    else:
        theformulas = formulas
    if not isinstance(concentrations,list):
        theconcentrations = []
        theconcentrations.append(concentrations)
    else:
        theconcentrations = concentrations
    if not isinstance(densities,list):
        thedensities = []
        thedensities.append(densities)
    else:
        thedensities = densities
    # get mu
    mu_tot = np.zeros((len(en)))
    for n in range(len(theformulas)):
         mu_tot += mpr(en,theformulas[n])[0]*theconcentrations[n]*thedensities[n]
    return np.exp(-mu_tot*thickness)

def plottransmission(energy,formulas,concentrations,densities,thickness):
    """
    opens a plot with the transmission plotted along the given energy vector
    """
    if not isinstance(formulas,list):
        theformulas = []
        theformulas.append(formulas)
    else:
        theformulas = formulas
    if not isinstance(concentrations,list):
        theconcentrations = []
        theconcentrations.append(concentrations)
    else:
        theconcentrations = concentrations
    if not isinstance(densities,list):
        thedensities = []
        thedensities.append(densities)
    else:
        thedensities = densities
    transmission = gettransmission(energy,formulas,concentrations,densities,thickness)
    plt.plot(energy,transmission)
    titlestring = 'transmission of: ' + ' '.join(formulas)
    plt.title(titlestring)
    plt.xlabel('energy [keV]')
    plt.ylabel('transmission [%]')
    plt.grid(False)
    plt.show()

def getpenetrationdepth(energy,formulas,concentrations,densities):
    """
    returns the penetration depth of a mixture of chemical formulas
    with certain concentrations and densities
    """
    en  = np.array([]) # turn energy into an iterable array
    en  = np.append(en,energy)
    if not isinstance(formulas,list):
        theformulas = []
        theformulas.append(formulas)
    else:
        theformulas = formulas
    if not isinstance(concentrations,list):
        theconcentrations = []
        theconcentrations.append(concentrations)
    else:
        theconcentrations = concentrations
    if not isinstance(densities,list):
        thedensities = []
        thedensities.append(densities)
    else:
        thedensities = densities
    # get mu
    mu_tot = np.zeros((len(en)))
    for n in range(len(theformulas)):
         mu_tot += mpr(en,theformulas[n])[0]*theconcentrations[n]*thedensities[n]
    return 1.0/mu_tot

def plotpenetrationdepth(energy,formulas,concentrations,densities):
    """
    opens a plot window of the penetration depth of a mixture of chemical formulas
    with certain concentrations and densities plotted along the given energy vector
    """
    if not isinstance(formulas,list):
        theformulas = []
        theformulas.append(formulas)
    else:
        theformulas = formulas
    if not isinstance(concentrations,list):
        theconcentrations = []
        theconcentrations.append(concentrations)
    else:
        theconcentrations = concentrations
    if not isinstance(densities,list):
        thedensities = []
        thedensities.append(densities)
    else:
        thedensities = densities
    pendepth = getpenetrationdepth(energy,formulas,concentrations,densities)
    plt.plot(energy,pendepth)
    titlestring = 'penetration depth of: ' + ' '.join(formulas)
    plt.title(titlestring)
    plt.xlabel('energy [keV]')
    plt.ylabel('penetration depth [cm]')
    plt.grid(False)
    plt.show()

def sumx(A):
    """
    Short-hand command to sum over 1st dimension of a N-D matrix (N>2) and to squeeze it to N-1-D matrix.
    """
    return np.squeeze(np.sum(A,axis=0))

def momtrans_au(e1,e2,tth):
    """
    Calculates the momentum transfer in atomic units
    input: 
    e1  = incident energy  [keV]    
    e2  = scattered energy [keV]
    tth = scattering angle [deg]
    returns:
    q   = momentum transfer [a.u.] (corresponding to sin(th)/lambda)
    """
    e1    = np.array(e1*1e3/13.60569172/2)
    e2    = np.array(e2*1e3/13.60569172/2)
    th    = np.radians(tth)#tth/180.0*numpy.pi
    hbarc = 137.03599976
    q     = 1/hbarc*np.sqrt(e1**2.0+e2**2.0-2.0*e1*e2*np.cos(th));
    return q

def momtrans_inva(e1,e2,tth):
    """
    Calculates the momentum transfer in inverse angstrom
    input: 
    e1  = incident energy  [keV]    
    e2  = scattered energy [keV]
    tth = scattering angle [deg]
    returns:
    q   = momentum transfer [a.u.] (corresponding to sin(th)/lambda)
    """
    e = 1.602e-19
    c = 2.9979e8 
    hbar = 6.626e-34/2/np.pi

    e1    = np.array(e1*1e3*e/c/hbar)
    e2    = np.array(e2*1e3*e/c/hbar)
    th    = np.radians(tth)
    q     = np.sqrt(e1**2+e2**2-2*e1*e2*np.cos(th))/1e10
    return q

def energy_monoangle(angle,d=5.4307/np.sqrt(11)):
    """
    % ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
    %         function e=energy(dspace,bragg_angle)
    %
    %         dspace for reflection (defaulf for Si(311) reflection)
    %         bragg_angle in DEG
    %
    %         KH 28.09.93
    %
    """
    hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
    e  = (2.0*d*np.sin(angle/180.0*np.pi)/hc)**(-1.0)
    return e

def find_center_of_mass(x,y):
    """
    Returns the center of mass (first moment) for the given curve y(x)
    """
    deno = np.trapezoid(y,x)
    if deno==0.0:
        return 0.0
        # print "*** print_tb:"
        # traceback.print_stack()
        # print " DENO==0!"
        # return 0.0
    return np.trapezoid(y*x,x)/deno

def is_allowed_refl_fcc(H):
    """ **is_allowed_refl_fcc**
    Check if given reflection is allowed for a FCC lattice.

    Args:
      * H (array, list, tuple): H=[h,k,l]
    
    Returns:
      * boolean
    
    """
    h = H[0]
    k = H[1]
    l = H[2]
    if h%2==0.0 and k%2==0.0 and l%2==0.0:
        answer = True
    elif (h+k+l)/4.0%1==0.0:
        answer = True
    elif h%2==1.0 and k%2==1.0 and l%2==1.0:
        answer = True
    else:
        answer = False
    return answer

def TTsolver1D(el_energy, hkl=[6,6,0], crystal='Si', R=1.0, dev=np.arange(-50.0,150.0,1.0), alpha=0.0, chitable_prefix='/home/christoph/sources/XRStools/data/chitables/chitable_'):
    """ **TTsolver**
    Solves the Takagi-Taupin equation for a bent crystal.
    
    This function is based on a Matlab implementation by S. Huotari of M. Krisch's
    Fortran programs.

    Args:
      * el_energy (float): Fixed nominal (working) energy in keV.
      * hkl (array): Reflection order vector, e.g. [6, 6, 0]
      * crystal (str): Crystal used (can be silicon 'Si' or 'Ge')
      * R (float): Crystal bending radius in m.
      * dev (np.array): Deviation parameter (in arc. seconds) for 
        which the reflectivity curve should be calculated.
      * alpha (float): Crystal assymetry angle.

    Returns:
      * refl (np.array): Reflectivity curve.
      * e (np.array): Deviation from Bragg angle in meV.
      * dev (np.array): Deviation from Bragg angle in microrad.

    """
    # load dielectric susceptibility data (tabulated)
    chi = np.loadtxt(chitable_prefix + crystal.lower() + str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2])) + '.dat')
    if len(chi[:,0]) == 1:
        print( 'Will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!')
    else:
        if el_energy < np.min(chi[:,0]) or el_energy > np.max(chi[:,0]):
            print( 'Energy outside of values defined in Chi-table.')
            return
    # interpolate 
    chi0 = complex(np.interp(el_energy,chi[:,0],chi[:,1]),np.interp(el_energy,chi[:,0],chi[:,2]))
    chih = complex(np.interp(el_energy,chi[:,0],chi[:,3]),np.interp(el_energy,chi[:,0],chi[:,4]))
    # set the stress tensor values depending on crystal used
    if crystal.upper() == 'SI':
        s13 = -0.278
    elif crystal.upper() == 'GE':
        s13 = -0.273
    else:
        print( 'Poisson ratio for this crystal not defined')
        return
    s15 = -0.0 # s15/s11
    # scattering angle in degree
    th  = braggd(hkl,el_energy,crystal)
    # dspace in m
    dsp = dspace(hkl,crystal)/10.0*1e-9
    # wavelength in m
    lam = 12.3984191/el_energy/10.0*1e-9
    # debye-waller factor
    dwf = 1.0 # dwf = 0.899577 
    # meridional bending radius
    radius = R 
    # sagittal bending radius
    rsag   = R*np.sin(np.radians(th))**2.0 
    # thickness in m
    thick  = 500.0*1e-6 
    # asymmetry in radians
    alpha    = np.radians(alpha)
    # deviation parameter in arcsec
    dev      = dev/3600.0/180.0*np.pi
    # gamma0,gammah = cos<n,K_0,h > , n = inward normal of crystal surface, K_0,h = wave vector
    gammah = -np.sin(np.arcsin(lam/(2.0*dsp)) + alpha) # Krisch et al. convention
    gamma0 =  np.sin(np.arcsin(lam/(2.0*dsp)) - alpha) # Krisch et al. convention
    gamma  = gammah/gamma0
    a0 = np.sqrt(1-gamma0**2.0)
    ah = np.sqrt(1-gammah**2.0)
    beta  = gamma0/np.abs(gammah)
    # polarization factor
    cpol   = 1.0 
    # penetration depth
    mu     = -2.0*np.pi/lam*chi0.imag
    tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))
    lex    = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chih.real)
    y0     = chi0.imag*(1.0+beta)/(2.0*np.sqrt(beta)*chih.real)
    c1     = cpol*dwf* complex(1.0,-chih.imag/chih.real)
    #abbreviation concerning the deviation parameter y
    abb0 = -np.sqrt(beta)/2.0/chih.real
    abb1 = chi0.real*(1.0+beta)/(2.0*np.sqrt(beta)*chih.real)
    #abbreviations concerning the deformation field
    abb2 = gamma0*gammah*(gamma0-gammah)
    abb3 = 1.0 + 1.0/(gamma0*gammah)
    abb4 = s13*(1.0 + radius/rsag)
    abb5 = (ah - a0)/(gamma0 - gammah)*s15 
    abb6 = 1.0/(np.abs(cpol)*chih.real*np.cos(np.arcsin(lam/(2.0*dsp)))*radius)
    abb7 = 2.0*np.abs(cpol)*chih.real*np.cos(np.arcsin(lam/(2.0*dsp)))/gamma0
    # spherical diced crystal, 1-m bending radius, nearly backscattering conditions, strain gradient
    sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))
    # number of steps along reflectivity curve
    nstep=len(dev)
    # reflectivity curve
    refl  = np.zeros_like(dev)

    OLDMETHOD = 0
    if OLDMETHOD:
        # loop over all steps along the reflectivity curve
        for l in range(nstep):
            # deviation parameter
            abb8   = -2.0*np.sin(2.0*np.arcsin(lam/(2.0*dsp)))*dev[l]
            T = np.arange(np.max([-10.0*tdepth, -thick]),0.0,1e-8)
            Y = odeint(odefctn,np.array([0.0, 0.0]),T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 
            # normalized reflectivity at this point
            refl[l] = np.sum(Y[-1,:]**2.0)
    else:
        # deviation parameter
        abb8   = -2.0*np.sin(2.0*np.arcsin(lam/(2.0*dsp)))*dev
        # dev axis (complex)
        YY = np.zeros([nstep],"D")
        # small step size
        ministep = tdepth/10000.0
        ssrk     = ministep/2
        xpoints  = np.arange(np.max([-10.0*tdepth, -thick]),0, ministep)
        for xpos in xpoints[:-1]:
            Yp0  = odefctn_CN( YY,              xpos+0*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Yp1  = odefctn_CN( YY+1.0*ssrk*Yp0, xpos+1*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Yp2  = odefctn_CN( YY+1.0*ssrk*Yp1, xpos+1*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Ypb  = odefctn_CN( YY+2.0*ssrk*Yp2, xpos+2*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            YY   =  YY  + ministep*(  Yp0 + 2.0*(Yp1+Yp2) + Ypb  )/6.0
            refl = YY.real*YY.real+YY.imag*YY.imag

    # deviation in degree
    dev = dev/4.848136811e-06/3600.0 
    # deviation in meV
    e_meV   = (energy(dspace(hkl,crystal),th+dev)-el_energy)*1.0e6
    # deviation in microrad
    dev = dev*1.e3
    return refl, e_meV, dev


def odefctn_CN(yCN,t,abb0,abb1,abb7,abb8N,lex,sgbeta,y0,c1):
    fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8N + abb7*sgbeta*t) + abb1) + complex(0,y0))*(yCN) + c1*(1.0 + yCN* yCN) )
    return fcomp



def odefctn(y,t,abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1):
    """
    #%    [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
    #%    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
    #%    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
    #%    no options are set.   
    """
    #print 'shape of y is ' , np.shape(y), np.shape(t)
    fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8 + abb7*sgbeta*t) + abb1) + complex(0,y0))*(y[0] + complex(0,y[1])) + c1*(1.0 + (y[0] + complex(0,y[1]))**2.0))
    return fcomp.real,fcomp.imag


def odefctn_CN(yCN,t,abb0,abb1,abb7,abb8N,lex,sgbeta,y0,c1):
    fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8N + abb7*sgbeta*t) + abb1) + complex(0,y0))*(yCN) + c1*(1.0 + yCN* yCN) )
    return fcomp

def taupgen(e, hkl = [6,6,0], crystals = 'Si', R = 1.0, dev = np.arange(-50.0,150.0,1.0), alpha = 0.0):
    """
    % TAUPGEN          Calculates the reflectivity curves of bent crystals
    %
    % function [refl,e,dev]=taupgen_new(e,hkl,crystals,R,dev,alpha);
    %
    %              e = fixed nominal energy in keV
    %            hkl = reflection order vector, e.g. [1 1 1]
    %       crystals = crystal string, e.g. 'si' or 'ge'
    %              R = bending radius in meters
    %            dev = deviation parameter for which the 
    %                  curve will be calculated (vector) (optional)
    %          alpha = asymmetry angle 
    % based on a FORTRAN program of Michael Krisch
    % Translitterated to Matlab by Simo Huotari 2006, 2007
    % Is far away from being good matlab writing - mostly copy&paste from
    % the fortran routines. Frankly, my dear, I don't give a damn. 
    % Complaints -> /dev/null
    """
    path =  os.path.join(data_installation_dir,'chitable_') # prefix + 'data/chitables/chitable_' # path to chitables
    # load the according chitable (tabulated)
    hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
    filestring = path + crystals.lower() + hkl_string + '.dat'
    chi = np.loadtxt(filestring)

    # good for 1 m bent crystals in backscattering
    ystart = -50.0 # start value of angular range in arcsecs
    yend   = 150.0 # end value of angular range in arcsecs
    ystep  = 1.0   # step width in arcsecs

    if len(chi[:,0]) == 1:
        print( ' I will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!')
    else:
        if e < np.min(chi[:,0]) or e > np.max(chi[:,0]):
            print( 'Energy outside of the range in ' + filestring)
            return

        chi0r = np.interp(e,chi[:,0],chi[:,1])
        chi0i = np.interp(e,chi[:,0],chi[:,2])
        chihr = np.interp(e,chi[:,0],chi[:,3])
        chihi = np.interp(e,chi[:,0],chi[:,4])

    th = braggd(hkl,e,crystals)
    lam = 12.3984191/e/10.0 # wavelength in nm

    reflcorr = 0.0
    chi0 = complex(chi0r,chi0i)
    chih = complex(chihr,chihi)

    if crystals.upper() == 'SI':
        s13 = -0.278
    elif crystals.upper() == 'GE':
        s13 = -0.273
    else:
        print( 'Poisson ratio for this crystal not defined')
        return

    s15 = -0.0 # s15/s11
    dsp = dspace(hkl,crystals)/10.0 # dspace

    dwf    = 1.0 # dwf = 0.899577 # debye-waller factor
    radius = R # meridional bending radius
    rsag   = R*np.sin(np.radians(th))**2.0 # sagittal bending radius
    thick  = 500.0 # thickness in micrometers #rsag = R

    lam      = lam*1e-9
    dsp      = dsp*1e-9
    alpha    = np.radians(alpha) # alpha in rad
    thick    = thick*1e-6
    ystart   = ystart/3600.0/180.0*np.pi
    yend     = yend/3600.0/180.0*np.pi
    ystep    = ystep/3600.0/180*np.pi
    dev      = dev/3600.0/180.0*np.pi
    reflcorr = reflcorr/3600.0/180.0*np.pi

    thetab = np.arcsin(lam/(2.0*dsp))
    cpol   = 1.0 # cpol=0.5*(1+cos(2*thetab).^2) # cpol=cos(2*thetab).^2

    # gamma0 = sin(thetab+alpha) # normal convention
    # gammah = -sin(thetab-alpha) # normal convention
    gammah = -np.sin(thetab + alpha) # Krisch et al. convention (really!)
    gamma0 = np.sin(thetab - alpha) # Krisch et al. convention (I'm not kidding!!)

    beta  = gamma0/np.abs(gammah)
    gamma = gammah/gamma0

    a0 = np.sqrt(1-gamma0**2.0)
    ah = np.sqrt(1-gammah**2.0)

    mu = -2.0*np.pi/lam*chi0i

    tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))

    lex = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chihr)

    y0 = chi0i*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    pfried = -chihi/chihr

    c1 = cpol*dwf* complex(1.0,pfried)

    #abbreviation concerning the deviation parameter y
    abb0 = -np.sqrt(beta)/2.0/chihr
    abb1 = chi0r*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    #abbreviations concerning the deformation field

    abb2 = gamma0*gammah*(gamma0-gammah)
    abb3 = 1.0 + 1.0/(gamma0*gammah)
    abb4 = s13*(1.0 + radius/rsag)
    abb5 = (ah - a0)/(gamma0 - gammah)*s15 
    abb6 = 1.0/(np.abs(cpol)*chihr*np.cos(thetab)*radius)
    abb7 = 2.0*np.abs(cpol)*chihr*np.cos(thetab)/gamma0

    #   a spectrometer based on a spherical diced analyzer crystal with a 1-m bending radius in nearly backscattering conditions utilizing a strain gradient beta
    sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))

    nstep=len(dev)
    eta  = np.zeros_like(dev)
    abb8z = np.zeros_like(dev)
    refl  = np.zeros_like(dev)
    refl1 = np.zeros_like(dev)
    refl2 = np.zeros_like(dev)

    OLDMETHOD = 1

    if OLDMETHOD :
            for l in range(nstep):
                    # actual value of the deviation angle
                    # dev[l] = ystart + (l - 1)*ystep

                    # deviation parameter
                    abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
                    eta[l] = (dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/(np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih))
                    eta[l] = eta[l].real

                    ndiff = 2
                    xend = 0
                    x = np.max([-10.0*tdepth, -thick])
                    y = np.array([0.0, 0.0])
                    h = xend
                    abb8z[l] = abb8

                    # in this point call the subroutine
                    #     [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
                    #    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
                    #    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
                    #    no options are set.   
                    #print 'the fucking shape of y is ', np.shape(y)
                    T = np.arange(x,xend,1e-8)
                    Y = odeint(odefctn,y,T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 

                    # normalized reflectivity at this point
                    refl[l] = np.sum(Y[-1,:]**2.0)
                    refl1[l] = Y[-1,0]
                    refl2[l] = Y[-1,1]

    else:
            YY = np.zeros([nstep],"D")
            for l in range(nstep):
                    abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
                    eta[l] = ((dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/
                              (np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih)))
                    eta[l] = eta[l].real
                    abb8z[l] = abb8


            xend = 0
            x = np.max([-10.0*tdepth, -thick])
            ministep = tdepth/1000.0 ## when delta/beta is 100 we have still 10
            xpoints = np.arange(x,0, ministep)
            substep_RungeKutta = ministep/2
            ssrk  = substep_RungeKutta

            for xpos in xpoints[:-1]:
                    Yp0 = odefctn_CN( YY, xpos+0*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Yp1 = odefctn_CN( YY+1*ssrk*Yp0, xpos+1*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Yp2 = odefctn_CN( YY+1*ssrk*Yp1, xpos+1*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Ypb = odefctn_CN( YY+2*ssrk*Yp2, xpos+2*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)

                    YY =  YY  + ministep*(  Yp0 + 2*(Yp1+Yp2) + Ypb  )/6.0

            refl1 = YY.real
            refl2 = YY.imag
            refl  = refl1*refl1+refl2*refl2


    de = dev * e * 1.0e6 /np.tan(thetab)

    lam    = lam *1.0e+09        
    dsp    = dsp*1.0e+09        
    alpha  = alpha/np.pi*180.0        
    ystart = ystart/4.848136811e-06
    yend   = yend/4.848136811e-06   
    ystep  = ystep/4.848136811e-06
    dev    = dev/4.848136811e-06 # dev in arcsecs
    
    dev = dev/3600.0 # in degrees
    thb = th
    th  = thb + dev
    e0  = e
    e   = energy(dspace(hkl,crystals),th)-e0
    e = e*1e6

    dev = dev*3600.0 # back to arcsecs

    return refl,e,dev,e0

def taupgen_amplitude(e, hkl = [6,6,0], crystals = 'Si', R = 1.0, dev = np.arange(-50.0,150.0,1.0), alpha = 0.0):
    """
    % TAUPGEN          Calculates the reflectivity curves of bent crystals
    %
    % function [refl,e,dev]=taupgen_new(e,hkl,crystals,R,dev,alpha);
    %
    %              e = fixed nominal energy in keV
    %            hkl = reflection order vector, e.g. [1 1 1]
    %       crystals = crystal string, e.g. 'si' or 'ge'
    %              R = bending radius in meters
    %            dev = deviation parameter for which the 
    %                  curve will be calculated (vector) (optional)
    %          alpha = asymmetry angle 
    % based on a FORTRAN program of Michael Krisch
    % Translitterated to Matlab by Simo Huotari 2006, 2007
    % Is far away from being good matlab writing - mostly copy&paste from
    % the fortran routines. Frankly, my dear, I don't give a damn. 
    % Complaints -> /dev/null
    """
    path =  os.path.join(data_installation_dir,'chitable_') # prefix + 'data/chitables/chitable_' # path to chitables
    # load the according chitable (tabulated)
    hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
    filestring = path + crystals.lower() + hkl_string + '.dat'
    chi = np.loadtxt(filestring)

    # good for 1 m bent crystals in backscattering
    ystart = -50.0 # start value of angular range in arcsecs
    yend   = 150.0 # end value of angular range in arcsecs
    ystep  = 1.0   # step width in arcsecs

    if len(chi[:,0]) == 1:
        print( ' I will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!')
    else:
        if e < np.min(chi[:,0]) or e > np.max(chi[:,0]):
            print( 'Energy outside of the range in ' + filestring)
            return

        chi0r = np.interp(e,chi[:,0],chi[:,1])
        chi0i = np.interp(e,chi[:,0],chi[:,2])
        chihr = np.interp(e,chi[:,0],chi[:,3])
        chihi = np.interp(e,chi[:,0],chi[:,4])

    th = braggd(hkl,e,crystals)
    lam = 12.3984191/e/10.0 # wavelength in nm

    reflcorr = 0.0
    chi0 = complex(chi0r,chi0i)
    chih = complex(chihr,chihi)

    if crystals.upper() == 'SI':
        s13 = -0.278
    elif crystals.upper() == 'GE':
        s13 = -0.273
    else:
        print( 'Poisson ratio for this crystal not defined')
        return

    s15 = -0.0 # s15/s11
    dsp = dspace(hkl,crystals)/10.0 # dspace

    dwf    = 1.0 # dwf = 0.899577 # debye-waller factor
    radius = R # meridional bending radius
    rsag   = R*np.sin(np.radians(th))**2.0 # sagittal bending radius
    thick  = 500.0 # thickness in micrometers #rsag = R

    lam      = lam*1e-9
    dsp      = dsp*1e-9
    alpha    = np.radians(alpha) # alpha in rad
    thick    = thick*1e-6
    ystart   = ystart/3600.0/180.0*np.pi
    yend     = yend/3600.0/180.0*np.pi
    ystep    = ystep/3600.0/180*np.pi
    dev      = dev/3600.0/180.0*np.pi
    reflcorr = reflcorr/3600.0/180.0*np.pi

    thetab = np.arcsin(lam/(2.0*dsp))
    cpol   = 1.0 # cpol=0.5*(1+cos(2*thetab).^2) # cpol=cos(2*thetab).^2

    # gamma0 = sin(thetab+alpha) # normal convention
    # gammah = -sin(thetab-alpha) # normal convention
    gammah = -np.sin(thetab + alpha) # Krisch et al. convention (really!)
    gamma0 = np.sin(thetab - alpha) # Krisch et al. convention (I'm not kidding!!)

    beta  = gamma0/np.abs(gammah)
    gamma = gammah/gamma0

    a0 = np.sqrt(1-gamma0**2.0)
    ah = np.sqrt(1-gammah**2.0)

    mu = -2.0*np.pi/lam*chi0i

    tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))

    lex = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chihr)

    y0 = chi0i*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    pfried = -chihi/chihr

    c1 = cpol*dwf* complex(1.0,pfried)

    #abbreviation concerning the deviation parameter y
    abb0 = -np.sqrt(beta)/2.0/chihr
    abb1 = chi0r*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    #abbreviations concerning the deformation field

    abb2 = gamma0*gammah*(gamma0-gammah)
    abb3 = 1.0 + 1.0/(gamma0*gammah)
    abb4 = s13*(1.0 + radius/rsag)
    abb5 = (ah - a0)/(gamma0 - gammah)*s15 
    abb6 = 1.0/(np.abs(cpol)*chihr*np.cos(thetab)*radius)
    abb7 = 2.0*np.abs(cpol)*chihr*np.cos(thetab)/gamma0

    #   a spectrometer based on a spherical diced analyzer crystal with a 1-m bending radius in nearly backscattering conditions utilizing a strain gradient beta
    sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))

    nstep=len(dev)
    eta  = np.zeros_like(dev)
    abb8z = np.zeros_like(dev)
    refl  = np.zeros_like(dev)
    refl1 = np.zeros_like(dev)
    refl2 = np.zeros_like(dev)

    OLDMETHOD = 0

    if OLDMETHOD :
            for l in range(nstep):
                    # actual value of the deviation angle
                    # dev[l] = ystart + (l - 1)*ystep

                    # deviation parameter
                    abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
                    eta[l] = (dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/(np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih))
                    eta[l] = eta[l].real

                    ndiff = 2
                    xend = 0
                    x = np.max([-10.0*tdepth, -thick])
                    y = np.array([0.0, 0.0])
                    h = xend
                    abb8z[l] = abb8

                    # in this point call the subroutine
                    #     [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
                    #    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
                    #    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
                    #    no options are set.   
                    #print 'the fucking shape of y is ', np.shape(y)
                    T = np.arange(x,xend,1e-8)
                    Y = odeint(odefctn,y,T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 

                    # normalized reflectivity at this point
                    refl[l] = np.sum(Y[-1,:]**2.0)
                    refl1[l] = Y[-1,0]
                    refl2[l] = Y[-1,1]

    else:
            YY = np.zeros([nstep],"D")
            for l in range(nstep):
                    abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
                    eta[l] = ((dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/
                              (np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih)))
                    eta[l] = eta[l].real
                    abb8z[l] = abb8


            xend = 0
            x = np.max([-10.0*tdepth, -thick])
            ministep = tdepth/1000.0 ## when delta/beta is 100 we have still 10
            xpoints = np.arange(x,0, ministep)
            substep_RungeKutta = ministep/2
            ssrk  = substep_RungeKutta

            for xpos in xpoints[:-1]:
                    Yp0 = odefctn_CN( YY, xpos+0*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Yp1 = odefctn_CN( YY+1*ssrk*Yp0, xpos+1*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Yp2 = odefctn_CN( YY+1*ssrk*Yp1, xpos+1*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)
                    Ypb = odefctn_CN( YY+2*ssrk*Yp2, xpos+2*ssrk ,      abb0,abb1,abb7,abb8z,lex,sgbeta,y0,c1)

                    YY =  YY  + ministep*(  Yp0 + 2*(Yp1+Yp2) + Ypb  )/6.0

            refl1 = YY.real
            refl2 = YY.imag
            refl  = refl1*refl1+refl2*refl2


    de = dev * e * 1.0e6 /np.tan(thetab)

    lam    = lam *1.0e+09        
    dsp    = dsp*1.0e+09        
    alpha  = alpha/np.pi*180.0        
    ystart = ystart/4.848136811e-06
    yend   = yend/4.848136811e-06   
    ystep  = ystep/4.848136811e-06
    dev    = dev/4.848136811e-06 # dev in arcsecs
    
    dev = dev/3600.0 # in degrees
    thb = th
    th  = thb + dev
    e0  = e
    e   = energy(dspace(hkl,crystals),th)-e0
    e = e*1e6

    dev = dev*3600.0 # back to arcsecs

    return refl,e,dev,e0, refl1, refl2


def hlike_Rwfn(n,l,r,Z):
    """ **hlike_Rwfn**
    Returns an array with the radial part of a hydrogen-like 
    wave function.

    Args:
      * n (integer): main quantum number n
      * l (integer): orbitalquantum number l
      * r (array): vector of radii on which the function should be evaluated
      * Z (float): effective nuclear charge
    """
    import math
    from scipy import special
    a0 = 0.52917721092
    factor1 = np.sqrt(   (2.0*Z/(n*a0))**3 * math.factorial((n-l-1.0))/(2.0*n*math.factorial((n+1.0)) )   )
    factor2 = (2.0*Z*r/(n*a0))**(l)
    factor3 = np.exp(-(Z*r/(n*a0)))
    lag = special.eval_genlaguerre(n-l-1.0,2.0*l+1.0,2.0*Z*r/(n*a0))
    return factor1*factor2*factor3*lag#*np.sqrt(n+1.0)

def compute_matrix_elements(R1,R2,k,r):
    # for ii=1:length(q);
    #     fun=y3d.^2.*besselj(4,q(ii)*r);
    #     int4(ii)=simpson(r,fun);
    # end    
    from scipy import special
    q        = np.linspace(0,30,len(r))
    r2RsphBR = np.linspace(0,30,len(r))
    for ii in range(len(q)):
        sphB = np.zeros_like(q)
        for jj in range(len(sphB)):
            # special.sph_jn returns the function in [0] and the derivative in [1] 
            sphB[jj] = special.sph_jn(k,q[ii]*r[jj])[0][-1]
        fun  = r**2*R1*sphB*R2
        r2RsphBR[ii] = np.trapezoid(fun,r)
    return q,r2RsphBR

def read_dft_wfn(element, n, l, spin=None, directory=data_installation_dir):
    """ **read_dft_wfn**
    Parses radial parts of wavefunctions.

    Args:
      * element (str): Element symbol.
      * n (int): Main quantum number.
      * l (int): Orbital quantum number.
      * spin (str): Which spin channel, default is average over up and down.
      * directory (str): Path to directory where the wavefunctions can be found.

    Returns:
      * r (np.array): radius
      *  wfn (np.array):
    """
    element_name = element.lower()
    subfolder    = 'ae'
    l_name       = ['s', 'p', 'd', 'f', 'g'][l]
    prefix       = 'wf-%d%s_'%(n, l_name)
    postfix1     = 'up'
    postfix2     = 'dn'
    path1 = os.path.join(directory, 'wave_functions', element_name, subfolder, prefix+postfix1)
    path2 = os.path.join(directory, 'wave_functions', element_name, subfolder, prefix+postfix2)

    # load
    au2a = constants.physical_constants['atomic unit of length'][0]*1e10
    wfn1 = np.loadtxt(path1)
    wfn2 = np.loadtxt(path2)
    r    = wfn1[:,0]
    wfn  = np.zeros_like(r)
    if not spin:
        wfn = wfn1[:,1]/2. + wfn2[:,1]/2.
    elif spin == 'up':
        wfn = wfn1[:,1]
    elif spin == 'dn' or spin == 'down':
        wfn = wfn2[:,1]
    else:
        print('unknown keyword for spin') # raise proper error
        return
    # normalize
    norm =  np.trapezoid( r**2 * wfn*np.conj(wfn), r )
    wfn /= norm
    return r, wfn




def readfio(prefix, scannumber, repnumber=0):
    """
    if repnumber = 0:
    reads a spectra-file (name: prefix_scannumber.fio)
    if repnumber > 1:
    reads a spectra-file (name: prefix_scannumber_rrepnumber.fio)
    """
    suffix = '.fio'
    filename = prefix + '%05d' % scannumber + suffix
    if repnumber > 0:
        filename = prefix + '%05d' % scannumber + 'r' + '%d' % repnumber + suffix

    # analyze structure of file    
    fid     = open(filename,'r')

    colnameflag  = ' Col'
    colstartflag = '%d'    
    colnames = []

    linenum  = 0
    for line in fid:
        linenum +=1
        if colnameflag in line: 
            colnames.append(line.strip())
        if colnameflag in line: 
            startline = linenum
    fid.close()
    thefile = open(filename,'r').readlines()
    data = []
    for line in thefile[(len(colnames)+startline):]:
        data.append([float(x) for x in line.strip().split()])

    return np.array(data), colnames

def energy_monoangle(angle,d=5.4307/np.sqrt(11)):
    """
    % ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
    %         function e=energy(dspace,bragg_angle)
    %
    %         dspace for reflection (defaulf for Si(311) reflection)
    %         bragg_angle in DEG
    %
    %         KH 28.09.93
    %
    """
    hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
    e  = (2.0*d*math.sin(angle/180.0*np.pi)/hc)**(-1.0)
    return e





def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.

    Parameters: 
      * y : array_like, shape (N,)
        the values of the time history of the signal.
      * window_size : int
        the length of the window. Must be an odd integer number.
      * order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
      * deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)

    Returns
      * ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).

    Notes:
      The Savitzky-Golay is a type of low-pass filter, particularly
      suited for smoothing noisy data. The main idea behind this
      approach is to make for each point a least-square fit with a
      polynomial of high order over a odd-sized window centered at
      the point.

    Examples ::

      t = np.linspace(-4, 4, 500)
      y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
      ysg = savitzky_golay(y, window_size=31, order=4)
      import matplotlib.pyplot as plt
      plt.plot(t, y, label='Noisy signal')
      plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
      plt.plot(t, ysg, 'r', label='Filtered signal')
      plt.legend()
      plt.show()

    References ::
       .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
          Data by Simplified Least Squares Procedures. Analytical
          Chemistry, 1964, 36 (8), pp 1627-1639.
       .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
          W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
          Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError :
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid'), signal.fftconvolve(Z, -c, mode='valid')






def readp01image(filename):
    """
    reads a detector file from PetraIII beamline P01
    """
    dim = 256
    f = open(filename,'rb')
    data = np.fromfile(f,np.int32)
    image = np.reshape(data,(dim,dim))
    f.close()
    return image

def readp01scan(prefix,scannumber):
    """
    reads a whole scan from PetraIII beamline P01 (experimental)
    """
    print ("parsing files of scan No. %s" % scannumber)
    #fioname = prefix + 'online/hasylab_' + "%05d" % scannumber + '.fio'
    fioprefix = prefix + 'online/ixs_scan_'
    fiodata = readfio(fioprefix,scannumber)[0]

    mats1 = np.zeros((np.shape(fiodata)[0],256,256))
    mats2 = np.zeros((np.shape(fiodata)[0],256,256))
    mats  = np.zeros((np.shape(fiodata)[0],256,256*2))

    for n in range(np.shape(fiodata)[0]):
        matname1 = prefix + 'ixs_scan_' + "%05d" % scannumber + '/mdpxa/ixs_scan_' + "%05d" % scannumber + '_a_' + "%05d" % (n+1)
        matname2 = prefix + 'ixs_scan_' + "%05d" % scannumber + '/mdpxa/ixs_scan_' + "%05d" % scannumber + '_b_' + "%05d" % (n+1)

        mats1[n,:,:] = readp01image(matname1)
        mats2[n,:,:] = readp01image(matname2)
        mats[n,:,0:256] = mats1[n,:,:]
        mats[n,:,256:]  = mats2[n,:,:]
    return fiodata, mats, mats1, mats2

def readp01scan_rep(prefix,scannumber,repetition):
    """
    reads a whole scan with repititions from PetraIII beamline P01 (experimental)
    """
    print ("parsing files of scan No. %s" % scannumber)
    #fioname = prefix + 'online/hasylab_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '.fio'
    fioprefix = prefix + 'online/ixs_scan_'
    fiodata = readfio(fioprefix,scannumber,repetition)[0]

    mats1 = np.zeros((np.shape(fiodata)[0],256,256))
    mats2 = np.zeros((np.shape(fiodata)[0],256,256))
    mats  = np.zeros((np.shape(fiodata)[0],256,256*2))

    for n in range(np.shape(fiodata)[0]):
        matname1 = prefix + 'ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition  + '/mdpxa/ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '_a_' + "%05d" % (n+1)
        matname2 = prefix + 'ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition  + '/mdpxa/ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '_b_' + "%05d" % (n+1)

        mats1[n,:,:] = readp01image(matname1)
        mats2[n,:,:] = readp01image(matname2)
        mats[n,:,0:256] = mats1[n,:,:]
        mats[n,:,256:]  = mats2[n,:,:]
    return fiodata, mats, mats1, mats2


def split_hdf5_address(dataadress):
    pos = dataadress.rfind(":")
    if ( pos==-1):
        raise Exception( """
roiaddress   must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
""")
    filename, groupname = dataadress[:pos], dataadress[pos+1:]
    return filename, groupname





