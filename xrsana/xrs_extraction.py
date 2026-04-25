#!/usr/bin/python
# Filename: xrs_extraction.py

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

import os
import copy
import numpy as np
from . import xrs_utilities, math_functions, xrs_ComptonProfiles
import matplotlib.pyplot as plt

from scipy import optimize
from .math_functions import *

installation_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources" )

debug = 0

if not debug:
    # when testing it from the source you can create a link in ./ to ../data/
    HFCP_PATH = os.path.join(installation_dir,'data/ComptonProfiles.dat')
    LOGTABLE_PATH = os.path.join(installation_dir,'data/logtable.dat')
else:
    HFCP_PATH     = 'data/ComptonProfiles.dat' #'/home/christoph/sources/XRStools/data/ComptonProfiles.dat'
    LOGTABLE_PATH = 'data/logtables.dat' #'/home/christoph/sources/XRStools/data/logtable.dat'


def map_chamber_names(name):
    """
    **map_chamber_names**
    Maps names of chambers to range of ROI numbers.
    """
    chamberNames = {'VD':list(range(0,12)),
                    'VU':list(range(12,24)),
                    'VB':list(range(24,36)),
                    'HR':list(range(36,48)),
                    'HL':list(range(48,60)),
                    'HB':list(range(60,72))}
    return chamberNames[name.upper()]


class HF_dataset:
    """
    **dataset**
    A class to hold all information from HF Compton profiles necessary to subtract background from the experiment.
    """
    def __init__(self, data, formulas, stoich_weights, edges):
        self.formulas       = formulas
        self.stoich_weights = stoich_weights
        self.edges          = edges #e.g. {'Li':['K','L23'], 'O':'K'}
        self.E0             = data.E0
        self.cenom          = data.cenom
        self.tth            = data.tth
        self.eloss          = data.eloss
        self.HFProfile      = xrs_ComptonProfiles.HFProfile(formulas, stoich_weights, HFCP_PATH)
        self.HFProfile.get_elossProfiles(self.E0,self.tth)

        # interpolate total HF profiles onto experimental eloss scale
        self.J_total   = np.zeros((len(self.eloss),len(self.tth)))
        self.C_total   = np.zeros((len(self.eloss),len(self.tth)))
        self.V_total   = np.zeros((len(self.eloss),len(self.tth)))
        self.q_vals    = np.zeros((len(self.eloss),len(self.tth)))
        for ii in range(len(self.tth)):
            self.J_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.J_total[:,ii])
            self.C_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.C_total[:,ii])
            self.V_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.V_total[:,ii])
            self.q_vals[:,ii]  = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.q_vals[:,ii])

        # initialize double dicts for {'element1':{'edge1','edge2',...}, 'element2':{'edge1','edge2',...} }
        self.C_edges = {}
        for key in self.edges:
            self.C_edges[key] = {}
            for edge in self.edges[key]:
                self.C_edges[key][edge] = np.zeros((len(self.eloss),len(self.tth)))

        # interpolate the core profiles for the desired elements
        for key in self.edges:
            for edge in self.edges[key]:
                edge_keyword = xrs_ComptonProfiles.mapShellNames(edge,xrs_utilities.element(key))
                for formula in self.formulas:
                    if key in formula:
                        # cp core-edge profile
                        atom_profile = self.HFProfile.FormulaProfiles[formula].AtomProfiles[key]
                        for jj in range(len(self.tth)):
                            self.C_edges[key][edge][:,jj] = np.interp(self.eloss,atom_profile.eloss,atom_profile.CperShell[edge_keyword][:,jj])
                    else:
                        print('Could not find ' + key + ' in any of the provided formulas.')

    def get_J_total_av(self,columns):
        return np.mean(self.J_total[:,columns])

    def get_C_total(self,columns):
        return np.mean(self.J_total[:,columns])

    def get_C_edges_av(self,element,edge,columns):
        return np.mean(self.C_edges[element][edge][:,columns])

class valence_CP:
    """
    **valence_CP**
    Class to organize information about extracted experimental valence Compton profiles.
    """
    def __init__(self):
        self.pzscale        = np.flipud(np.arange(-10,10,0.05)) # definition according to Huotari et al, JPCS 62 (2001) 2205
        self.valencepz      = np.array([])
        self.valasymmetrypz = np.array([])
        self.valence        = np.array([])
        self.valasymmetry   = np.array([])

    def get_asymmetry(self):
        pass

    def get_pzscale(self):
        pass

class functorObjectV:

    def __init__( self, y, eloss, hfcore, lam ):
        self.y      =  y
        self.eloss  =  eloss
        self.hfcore =  hfcore
        self.lam    =  lam

    def funct( self, a, eloss ):
        pear = math_functions.pearson7_zeroback( eloss, a )
        poly = np.polyval( a[4:6], eloss )
        tot  = pear+poly
        return tot

    def __call__( self, x ):
        pea = math_functions.pearson7_zeroback( self.eloss, x[0:4] )
        if len(x)==7:
            pol = np.polyval(x[4:6],self.eloss)
            hf  = self.hfcore
        else:
            pol = 0
            hf  = 0

        self.hf_fit = hf
        self.fit    = pea+pol+hf
        self.peapol = pea+pol

        diff = self.y*x[6]-self.fit
        # with extra cost for large linear background
        res  = diff/np.sqrt(self.y*x[6]) + self.lam*(x[4]**2+x[5]**2)

        return res

class edge_extraction:
    """
    **edge_extraction**
    Class to destill core edge spectra from x-ray Raman scattering experiments.
    """
    def __init__(self,exp_data, formulas, stoich_weights, edges ,prenormrange=[5,np.inf]):
        # input
        self.eloss   = copy.deepcopy( exp_data.eloss )
        self.signals = copy.deepcopy( exp_data.signals )
        self.errors  = copy.deepcopy( exp_data.errors )
        self.E0      = exp_data.E0
        self.tth     = exp_data.tth
        self.prenormrange = prenormrange
        self.HF_dataset = HF_dataset(exp_data, formulas, stoich_weights, edges)

        # output
        self.background     = np.zeros(np.shape(exp_data.signals))
        self.sqw            = np.zeros(np.shape(exp_data.signals))
        self.sqwav          = np.zeros(np.shape(exp_data.eloss))
        self.sqwaverr       = np.zeros(np.shape(exp_data.eloss))
        self.valence_CP     = valence_CP()

        # some variables for averaging rawdata over analyzers/whole chambers
        self.avsignals  = np.array([])
        self.averrors   = np.array([])
        self.av_C       = {}
        self.av_J       = np.array([])
        self.avqvals    = np.array([])

        # rough normalization over range given by prenormrange
        scales = []
        if prenormrange:
            for n in range(self.signals.shape[1]):
                HFnorm  = np.trapezoid(self.HF_dataset.J_total[:,n], self.eloss)
                inds    = np.where(np.logical_and(self.eloss >= prenormrange[0], self.eloss <= prenormrange[1]))[0]
                EXPnorm = np.trapezoid(self.signals[inds,n], self.eloss[inds])
                scales.append(EXPnorm*HFnorm)
                #print 'HFnorm = ', HFnorm, ', EXPnorm = ', EXPnorm
                #self.signals[:,n] /= EXPnorm
                #self.signals[:,n] *= HFnorm
                #self.errors[:,n]  /= EXPnorm
                #self.errors[:,n]  *= HFnorm
            scale = np.mean(scales)
            for n in range(self.signals.shape[1]):
                self.signals[:,n] *= scale
                self.errors[:,n]  *= scale

    def analyzerAverage(self,roi_numbers,errorweighing=True):
        """
        **analyzerAverage**
        Averages signals from several crystals before background subtraction.

        Args:

          * roi_numbers : list, str
              list of ROI numbers to average over of keyword for analyzer chamber
              (e.g. 'VD','VU','VB','HR','HL','HB')

          *  errorweighing : boolean (True by default)
                  keyword if error weighing should be used for the averaging or not

        """    
        if isinstance(roi_numbers,list):
            columns = roi_numbers
        elif isinstance(roi_numbers,str): 
            columns = map_chamber_names(roi_numbers)
        else:
            print('Unsupported type for keyword \'roi_numbers\'.')
            return

        self.avsignals = np.zeros_like(self.eloss)
        self.averror   = np.zeros_like(self.eloss)
        self.av_J      = np.zeros_like(self.eloss)
        self.avqvals   = np.zeros_like(self.eloss)

        # build matricies to sum over
        av      = np.zeros((len(self.eloss),len(columns)))
        averr   = np.zeros((len(self.eloss),len(columns)))
        avqvals = np.zeros((len(self.eloss),len(columns)))
        avJ     = np.zeros((len(self.eloss),len(columns)))
        avcvals = {}
        for key in self.HF_dataset.edges:
            avcvals[key] = {}
            for edge in self.HF_dataset.edges[key]:
                avcvals[key][edge] = np.zeros((len(self.eloss),len(columns)))

        # assign the signals to sum over
        for column,ii in zip(columns,list(range(len(columns)))):
            # find data points with error = 0.0 and replace error by 1.0
            inds = np.where(self.errors[:,column] == 0.0)[0]
            self.errors[inds,column] = 1.0
            av[:,ii]      = self.signals[:,column]
            averr[:,ii]   = self.errors[:,column]
            avqvals[:,ii] = self.HF_dataset.q_vals[:,column]
            for key in self.HF_dataset.edges:
                for edge in self.HF_dataset.edges[key]:
                    avcvals[key][edge][:,ii] = self.HF_dataset.C_edges[key][edge][:,column]

        # sum things up
        if errorweighing:
            self.avsignals = np.sum( av/(averr**2.0) ,axis=1)/( np.sum(1.0/(averr**2.0),axis=1))
            self.averrors  = np.sqrt( 1.0/np.sum(1.0/(averr**2.0),axis=1) )
        else: 
            self.avsignals = np.mean(av,axis=1)
            self.averrors  = np.sqrt(np.sum(np.absolute(averr)**2.0,axis=1))/np.sqrt(averr.shape[1]*(averr.shape[1]-1)) # check this again

        # average over HF core profiles 
        self.avqvals = np.mean(avqvals,axis=1)
        self.av_J    = np.mean(self.HF_dataset.J_total[:,columns],axis=1)
        for key in self.HF_dataset.edges:
            self.av_C[key] = {}
            for edge in self.HF_dataset.edges[key]:
                self.av_C[key][edge] = np.mean(avcvals[key][edge],axis=1)

    def removePolyCoreAv(self,element,edge,range1,range2,weights=[1,1],guess=[1.0,0.0,0.0],ewindow=100.0):
        """
        **removePolyCoreAv**
        Subtract a polynomial from averaged data guided by the HF core Compton profile.

        Args

          * element : str
            String (e.g. 'Si') for the element you want to work on.
          * edge: str
            String (e.g. 'K' or 'L23') for the edge to extract.
          * range1 : list
            List with start and end value for fit-region 1.
          * range2 : list
            List with start and end value for fit-region 2.
          * weigths : list of ints
            List with weights for the respective fit-regions 1 and 2. Default is [1,1].
          * guess : list
            List of starting values for the fit. Default is [1.0,0.0,0.0] (i.e. a quadratic
            function. Change the number of guess values to get other degrees of polynomials
            (i.e. [1.0, 0.0] for a constant, [1.0,0.0,0.0,0.0] for a cubic, etc.).
            The first guess value passed is for scaling of the experimental data to the HF
            core Compton profile.
          * ewindow: float
            Width of energy window used in the plot. Default is 100.0.

        """
        # check that there are averaged signals available
        if not np.any(self.avsignals):
            print('Found no averaged signals. Use \'analyzerAverage\'-method first to create some averages.')
            return

        # check that desired edge is available
        if not element in list(self.HF_dataset.edges.keys()):
            print('Cannot find HF profiles for desired atom.')
            return
        if not edge in self.HF_dataset.edges[element]:
            print('Cannot find HF core profiles for desired edge.' )
            return

        # define fitting ranges
        region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))
        region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # prepare plotting window
        plt.ion()
        plt.cla()

        # get the HF core spectrum
        HF_core = self.av_C[element][edge]

        # estimate start value for scaling parameter
        HF_core_norm = np.trapezoid(HF_core[region2],self.eloss[region2])
        exp_norm     = np.trapezoid(self.avsignals[region2],self.eloss[region2])
        self.avsignals *= HF_core_norm/exp_norm

        # define fit-function, boundaries, and constraints
        cons    =  ({'type': 'eq',   'fun': lambda x: np.trapezoid(HF_core[region2],self.eloss[region2]) - 
                                                        np.trapezoid(x[0]*self.avsignals[region2]-HF_core[region2] - 
                                                                np.polyval(x[1::],self.eloss[region2]),self.eloss[region2] )  },
                    {'type': 'eq',   'fun': lambda x: np.trapezoid( np.abs( x[0]*self.avsignals[region1] - HF_core[region1] - 
                                                                np.polyval(x[1::],self.eloss[region1]),self.eloss[region1] )) },
                    {'type': 'ineq', 'fun': lambda x: x[0]})

        fitfct = lambda a: np.sum( (a[0]*self.avsignals[region] - HF_core[region] - np.polyval(a[1::],self.eloss[region]) )**2.0 )
        res    = optimize.minimize(fitfct, guess, method='SLSQP',constraints=cons).x
        print( 'The fit parameters are: ', res)

        yres = np.polyval(res[1::], self.eloss)
        plt.plot(self.eloss,HF_core)    
        plt.plot(self.eloss,self.avsignals*res[0], self.eloss,yres+HF_core, self.eloss,self.avsignals*res[0]-yres)
        plt.legend(['HF core','scaled signal','poly-fit + core','scaled signal - poly'])
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(range1[0]-ewindow,range2[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav    = self.avsignals*res[0] - yres
        self.sqwaverr = self.averrors*res[0]

    def removeCorePearsonAv(self,element,edge,range1,range2,weights=[2,1],HFcore_shift=0.0,guess=None,scaling=None,return_background=False):
        """
        **removeCorePearsonAv**
        """
        # check that there are averaged signals available
        if not np.any(self.avsignals):
            print('Found no averaged signals. Use \'analyzerAverage\'-method first to create some averages.')
            return

        # check that desired edge is available
        if not element in list(self.HF_dataset.edges.keys()):
            print('Cannot find HF profiles for desired atom.') # @@@@@@@@@@@@@@@@@@@ raise
            return
        if not edge in self.HF_dataset.edges[element]:
            print('Cannot find HF core profiles for desired edge.' ) # @@@@@@@@@@@@ raise
            return

        # define fitting ranges
        region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))
        region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # find indices for guessing start values from HF J_total
        fitfct = lambda a: np.sum( (self.av_J[region] - pearson7_zeroback(self.eloss,a)[region] - 
                                    np.polyval(a[4:6],self.eloss[region]) )**2.0 )
        guess1 = optimize.minimize(fitfct, [1.0,1.0,1.0,1.0], method='SLSQP').x
        #guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]
        #if not guess: 
        #    guess = []
        #    ind   = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
        #    guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
        #    guess.append(guess[0]*1.0) # once the position of the peason maximum
        #    guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
        #    guess.append(self.avsignals[guessregion][ind]) # Peak intensity
        #    guess.append(0.0) # linear slope
        #    guess.append(0.0) # no background
        #    guess.append(1.0) # scaling factor for exp. data
        if not guess:
            guess = guess1
            guess = np.append(guess,[0.0,0.0,1.0]) # append starting values for linear and scaling

        # manage some plotting things
        plt.ion()
        plt.cla()

        # get the HF core spectrum
        HF_core = np.interp(self.eloss,self.eloss+HFcore_shift,self.av_C[element][edge])

        # approximately scale data in post-edge region
        HF_core_norm = np.trapezoid(HF_core[region2],self.eloss[region2])
        exp_norm     = np.trapezoid(self.avsignals[region2],self.eloss[region2])
        the_signals  = self.avsignals*HF_core_norm/exp_norm

        if scaling:
            bnds = ((None, None), (0, None), (0, None), (0, None), (0, None), (0, None), (scaling-scaling*0.001, scaling+scaling*0.001))
        else:
            scaling = 1.0
            bnds = ((None, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None))

        cons    =  (#{'type': 'ineq', 'fun': lambda x:  x[2]  },
                    #{'type': 'ineq', 'fun': lambda x:  x[3]  },
                    #{'type': 'ineq', 'fun': lambda x:  x[6]  },
                    {'type': 'eq',   'fun': lambda x:  np.trapezoid(np.abs(scaling*the_signals[region1] - 
                                                                pearson7_zeroback(self.eloss[region1],x[0:4]) - 
                                                                np.polyval(x[4:6],self.eloss[region1]) - 
                                                                HF_core[region1] ),
                                                                self.eloss[region1]  ) },
                    {'type': 'eq',   'fun': lambda x:  np.trapezoid(scaling*the_signals[region2] - 
                                                                pearson7_zeroback(self.eloss[region2],x[0:4]) - 
                                                                np.polyval(x[4:6],self.eloss[region2])-HF_core[region2],
                                                                self.eloss[region2])})

        fitfct = lambda a: np.sum( (scaling*the_signals[region] - 
                                    pearson7_zeroback(self.eloss[region],a[0:4]) - 
                                    np.polyval(a[4:6],self.eloss[region]) - 
                                    HF_core[region] )**2.0 )

        res    = optimize.minimize(fitfct, guess, method='SLSQP', bounds=bnds, constraints=cons).x
        print( 'The fit parameters are: ', res)

        yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)
        plt.plot(self.eloss,the_signals*scaling,self.eloss,yres+HF_core,self.eloss,the_signals*scaling-yres,self.eloss,HF_core)
        plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear)','core'))
        plt.draw()

        self.sqwav = the_signals*scaling - yres
        self.sqwaverr = self.averrors*scaling*HF_core_norm/exp_norm

        if return_background:
            return self.eloss, yres, HF_core

    def removePearsonAv(self,element,edge,range1,range2=None,weights=[2,1],guess=None,scale=1.0,HFcore_shift=0.0):
        """
        **removePearsonAv**
        """
        # check that there are averaged signals available
        if not np.any(self.avsignals):
            print('Found no averaged signals. Use \'analyzerAverage\'-method first to create some averages.')
            return

        # define fitting ranges
        region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))
        if range2:
            region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))
            region  = np.append(region1*weights[0],region2*weights[1])
        else:
            region  = region1

        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]
        if not guess: 
            guess = []
            ind   = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*2.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[guessregion][ind]) # Peak intensity
            guess.append(0.0) # no background

        # scale data by hand
        thespec = self.avsignals * scale

        # get the HF core spectrum
        HF_core = np.interp(self.eloss,self.eloss+HFcore_shift,self.av_C[element][edge])

        # define fitfunction
        fitfct  = lambda a: np.sum( (thespec[region] - pearson7_zeroback(self.eloss[region],a[0:4]) - np.polyval(a[4:6],self.eloss[region]) - HF_core[region])**2.0 )

        res = optimize.minimize(fitfct,guess).x
        print( 'the fitting results are: ', res)

        yres = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)
        plt.cla()
        plt.plot(self.eloss,thespec,self.eloss,yres,self.eloss,thespec-yres,self.eloss,HF_core)
        plt.legend(('data','pearson fit','data - pearson','core'))
        plt.draw()

        self.sqwav = thespec-yres
        self.sqwaverr = self.averrors * scale





        

    def save_average_Sqw(self,filename, emin=None, emax=None, normrange=None):
        """
        **save_average_Sqw**
        Save the S(q,w) into a ascii file (energy loss, S(q,w), Poisson errors).

        Args:
          * filename : str
            Filename for the ascii file.
          * emin : float
            Use this to save only part of the spectrum.
          * emax : float
            Use this to save only part of the spectrum.
          * normrange : list of floats
            E_start and E_end for possible area-normalization before saving.

        """
        # check that there are is an extracted S(q,w) available
        if not np.any(self.sqwav):
            print('Found no extracted S(q,w).')
            return

        if emin and emax: 
            inds = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))[0]
            data = np.zeros((len(inds),3))
            data[:,0] = self.eloss[inds]
            data[:,1] = self.sqwav[inds]
            data[:,2] = self.sqwaverr[inds]
        else:
            data = np.zeros((len(self.eloss),3))
            data[:,0] = self.eloss
            data[:,1] = self.sqwav
            data[:,2] = self.sqwaverr
        if normrange:
            #assert type(normrange) is list and len(normrange) is 2, "normrange has to be a list of length two!"
            assert isinstance(normrange, list) and len(normrange) == 2, "normrange has to be a list of length two!"
            inds = np.where(np.logical_and(data[:,0]>=normrange[0],data[:,0]<=normrange[1]))[0]
            norm = np.trapezoid(data[inds,1],data[inds,0])
            data[:,1] /= norm
            data[:,2] /= norm
        np.savetxt(filename,data)



    def removeCorePearsonAv_new( self, element, edge, range1, range2,
                                     HFcore_shift=0.0, guess=None, scaling=None,
                                     return_background=False, reg_lam=10):
        """
        **removeCorePearsonAv_new**
        """
        # check that there are averaged signals available
        if not np.any(self.avsignals):
            raise ValueError('Found no averaged signals. Use \'analyzerAverage\'-method first to create some averages.')

        # check that desired edge is available
        if not element in list(self.HF_dataset.edges.keys()):
            raise ValueError('Cannot find HF profiles for desired atom.')

        if not edge in self.HF_dataset.edges[element]:
            raise ValueError('Cannot find HF core profiles for desired edge.')

        # define fitting ranges
        region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))
        region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))
        region  = np.append(region1, region2)

        # get the HF core spectrum
        HF_core = np.interp( self.eloss, self.eloss+HFcore_shift, self.av_C[element][edge] )
        
        y_reg1       = self.avsignals[region1]
        eloss_reg1   = self.eloss[region1]
        HF_core_reg1 = HF_core[region1]
        
        y_reg2       = self.avsignals[region2]
        eloss_reg2   = self.eloss[region2]
        HF_core_reg2 = HF_core[region2]
        
        y_reg       = self.avsignals[region]
        eloss_reg   = self.eloss[region]
        HF_core_reg = HF_core[region]

        HF_reg_max = HF_core[region].max()
        y_reg_max  = y_reg.max()
        fact       =  HF_reg_max/y_reg_max

        if not guess:
            fitfct   = functorObjectV( y_reg, eloss_reg,  HF_core_reg , reg_lam )
            bndsa    = [ -np.inf]+ [ 0  for tmp in range(6)]
            bndsa[5] = -np.inf
            bndsb    = [ np.inf for tmp in range(7)]
            solution = optimize.least_squares(fitfct, np.random.rand(7),method='trf', bounds=[bndsa,bndsb] )
            guess = solution.x
            # print ( "Using following guess parameters: ", solution.x)

        # manage some plotting things
        plt.ion()
        plt.cla()

        # do the actual fit using the guess parameters
        fitfct_res = functorObjectV( y_reg, eloss_reg,  HF_core_reg , reg_lam )  
        bndsa      = [ -np.inf]+ [ 0  for tmp in range(6)]
        bndsa[5]   = -np.inf
        bndsb      = [ np.inf for tmp in range(7)]
        solution   = optimize.least_squares(fitfct_res, guess, method='trf', bounds=[bndsa,bndsb] )
        print ( "Best fit parameters: ", solution.x)

        # plot results
        fitfct_res = functorObjectV( self.avsignals, self.eloss,  HF_core , reg_lam  )
        res = fitfct_res(solution.x)

        x  = self.eloss                   # eloss
        y1 = self.avsignals*solution.x[6] # scaled data
        y2 = fitfct_res.fit               # pearson + linear + core
        y3 = y1 - fitfct_res.peapol       # data - (pearson + linear)
        y4 = HF_core                      # core

        plt.plot( x, y1, x, y2, x, y3, x, y4 )
        plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear)','core'))
        plt.draw()

        self.sqwav = y3
        self.sqwaverr = self.averrors*solution.x[6]

        if return_background:
            yres = pearson7_zeroback(self.eloss, solution.x[0:4]) + np.polyval(solution.x[4:6], self.eloss)
            return self.eloss, yres, HF_core





