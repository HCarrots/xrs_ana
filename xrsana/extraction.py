#!/usr/bin/python
# Filename: extraction.py

from .math_functions import *
from .xrs_utilities import *

import numpy as np
import pylab
import math
from scipy import interpolate, signal, integrate, constants, optimize, ndimage
import matplotlib.pyplot as plt

class extraction:
    """
    class for extraction of S(q,w) from an instance of the id20 class "data" and the predictthings class "theory"
    """
    def __init__(self,data,theory,prenormrange=[5,np.inf]):
        self.data=data
        # the data    
        self.eloss   = data.eloss
        self.signals = data.signals
        self.errors  = data.errors
        self.E0      = data.E0
        self.tth     = data.tth
        self.prenormrange = prenormrange
        # the theory 
        self.J        = theory.J
        self.C        = theory.C
        self.V        = theory.V
        self.qvals    = theory.q
        self.formulas = theory.formulas
        self.concentrations = theory.concentrations
        # output
        self.background = np.zeros(np.shape(data.signals))
        self.sqw        = np.zeros(np.shape(data.signals))
        self.pzscale      = np.flipud(np.arange(-10,10,0.05)) # definition according to Huotari et al, JPCS 62 (2001) 2205
        self.valencepz      = np.zeros((len(self.pzscale),len(self.signals[0,:])))
        self.valasymmetrypz = np.zeros((len(self.pzscale),len(self.signals[0,:])))
        self.valence      = np.zeros((len(self.eloss),len(self.signals[0,:])))
        self.valasymmetry = np.zeros((len(self.eloss),len(self.signals[0,:])))
        self.sqwav      = np.zeros(np.shape(data.eloss))
        self.sqwaverr   = np.zeros(np.shape(data.eloss))
        # some variables for averaging over analyzers/whole chambers
        self.avsignals  = np.array([])
        self.averrors   = np.array([])
        self.avC        = np.array([])
        self.avqvals    = np.array([])

        # rough normalization over range given by prenormrange
        for n in [ nn for nn in  range(len(self.signals[0,:])) if   ((self.signals[:,nn]).sum()>0)  ]:
        # for n in range(len(self.signals[0,:])):
            HFnorm = np.trapezoid(self.J[:,n],self.eloss)
            inds   = np.where(np.logical_and(self.eloss>=prenormrange[0],self.eloss<=prenormrange[1]))[0]
            EXPnorm = np.trapezoid(self.signals[inds,n],self.eloss[inds])
            self.signals[:,n] = self.signals[:,n]/EXPnorm*HFnorm
            self.errors[:,n]  = self.errors[:,n]/EXPnorm*HFnorm

    def areanorm(self,whichq,emin=None,emax=None):
        """
        normalizes self.signals to area in between emin and emax, 
        default values cover the whole self.eloss axis
        """        
        cols = []        
        if not isinstance(whichq,list):
            cols.append(whichq)
        else:
            cols = whichq

        if not emin:
            emin = self.eloss[0]
        if not emax:
            emax = self.eloss[-1]

        for col in cols:
            inds = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))
            self.signals[:,col] = self.signals[:,col]/np.trapezoid(self.signals[inds,col],self.eloss[inds])

    def analyzerAverage(self,whichq,errorweighing=True):
        """
        average signals over several analyzers before background subtraction, either with (default) or without error weighing whichq = either list of analyzer numbers to be averaged over or keywords ('VD','VB', ... ) to average over chambers
        """    
        analyzerNames = ['VD','VU','VB','HR','HL','HB']
        if type(whichq) == str:
            if whichq.upper() == 'VD':
                columns = list(range(0,12))
            elif whichq.upper() == 'VU':
                columns = list(range(12,24))
            elif whichq.upper() == 'VB':
                columns = list(range(24,36))
            elif whichq.upper() == 'HL':
                columns = list(range(36,48))
            elif whichq.upper() == 'HR':
                columns = list(range(48,60))
            elif whichq.upper() == 'HB':
                columns = list(range(60,72))
            elif whichq not in analyzerNames:
                print( 'Unknown keyword ' + '\'' + whichq + '\'' + '! Try one of these: \'VD\', \'VU\', \'VB\', \'HR\', \'HL\', \'HB\'.')
                return
        else:
            if not isinstance(whichq,list):
                columns = []
                columns.append(whichq)
            else: 
                columns = whichq

        # build the matricies
        av    = np.zeros((len(self.eloss),len(columns)))
        averr = np.zeros((len(self.eloss),len(columns)))
        avC   = np.zeros((len(self.eloss),len(columns)))
        avqvals = np.zeros((len(self.eloss),len(columns)))
        for n in [ nn for nn in  range(len(columns)) if  ((self.signals[:,nn]).sum()>0)   ]:
            # find data points with error = 0.0 and replace by 1.0
            inds = np.where(self.errors[:,columns[n]] == 0.0)[0]
            for ind in inds:
                self.errors[ind,columns[n]] = 1.0
            # arrange the desired columns into a matrix
            av[:,n]    = self.signals[:,columns[n]]
            averr[:,n] = self.errors[:,columns[n]]
            avC[:,n]   = self.C[:,columns[n]]
            avqvals[:,n] = self.qvals[:,columns[n]]
        # sum things up
        if errorweighing:
            self.avsignals = np.sum( av/averr**2.0 ,axis=1)/( np.sum(1.0/averr**2.0,axis=1))
            self.averrors  = np.sqrt( 1.0/np.sum(1.0/(averr)**2.0,axis=1) )
        else: 
            self.avsignals = np.sum(av,axis=1)
            self.averrors  = np.sqrt(np.sum(np.absolute(averr)**2.0,axis=1)) # check this again

        self.avC = np.mean(avC,axis=1)
        self.avqvals = np.mean(avqvals,axis=1)

    def energycorrect(self,whichq,alpha,densities,samthickness):
        """
        apply energy dependent corrections to the measured data based on scattering angles, sample material, ... 
        whichq       = single value or list of columns to apply the corrections to (index starts at zero)
        alpha        = incident beam angle (relative to sample surface normal), need to be negative for transmission geometry
        density      = single value (just one compound) or list of values (mixture of several compounds) 
        samthickness = sample thickness in [cm]        
        """
        # make things iterable
        cols = []        
        if not isinstance(whichq,list):
            cols.append(whichq)
        else:
            cols = whichq
        denses = []        
        if not isinstance(densities,list):
            denses.append(densities)
        else:
            denses = densities

        # calculate beta (exit angle) from alpha and tth for all columns in whichq
        beta = []
        
        if alpha >=0: # reflection geometry (check this again)
            for n in range(len(cols)):
                beta.append(np.absolute(180.0 - alpha - self.tth[cols[n]]))
        else: # transmission geometry (check this again)
            for n in range(len(cols)):
                beta.append(-1.0*(np.absolute(np.absolute(alpha) - self.tth[cols[n]])))

        # absorption and self-absorption
        # mu_in and mu_out from log-log table
        mu_in, mu_out = mpr_compds(self.eloss/1e3+self.E0,self.formulas,self.concentrations,self.E0,denses)
        ac = np.zeros((len(self.eloss),len(cols)))
        for n in range(len(cols)):
            ac[:,n] = abscorr2(mu_in,mu_out,alpha,beta[n],samthickness)

        # cross section correction (use cf output from e2pz routine)
        pz = np.zeros((len(self.eloss),len(cols)))
        cf = np.zeros((len(self.eloss),len(cols)))
        for n in range(len(cols)):
            pz[:,n], cf[:,n] = e2pz(self.E0+self.eloss/1.0e3,self.E0,self.tth[n])

        for col in cols:
            self.signals[:,col] = self.signals[:,col]*ac[:,n]*cf[:,n]
            self.errors[:,col]  = self.errors[:,col]*ac[:,n]*cf[:,n]
            # normalize back to HF profiles (because we know those are on eV / [1/eV] scale (is this correct)
            HFnorm = np.trapezoid(self.J[:,col],self.eloss)
            inds   = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]
            EXPnorm = np.trapezoid(self.signals[inds,col],self.eloss[inds])
            self.signals[:,col] = self.signals[:,col]/EXPnorm*HFnorm
            self.errors[:,col]  = self.errors[:,col]/EXPnorm*HFnorm

    def removeelastic(self,whichq,range1,range2,guess=None,stoploop=True,overwrite=False):
        """
        subtract a pearson function before starting extraction procedure, e.g. for subtracting the elastic peak tail
        guess values: 
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        a[4] = background
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        region1 = np.where(np.logical_and(self.eloss>=range1[0],self.eloss<=range1[1]))[0]
        region2 = np.where(np.logical_and(self.eloss>=range2[0],self.eloss<=range2[1]))[0]
        region  = np.append(region1,region2)

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = self.signals[:,col].argmax(axis=0) 
                guess.append(self.eloss[ind]) # max of signal (in range of prenorm from __init__)
                guess.append(1.0) # twice the position of the peason maximum
                guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(1e2) # Peak intensity
                guess.append(0.0) # background
            
            fitfct  = lambda a: self.signals[region,col] - pearson7(self.eloss[region],a)
            res     = optimize.leastsq(fitfct,guess)
            yres    = pearson7(self.eloss,res[0])

            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('data','pearson fit','data - pearson'))
            plt.draw()

            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
            if overwrite: 
                self.signals[:,col] = self.signals[:,col] - yres
        plt.ioff()

    def removeconst(self,whichq,emin,emax,ewindow=100.0,stoploop=True):
        """
        fits a constant as background in the range emin-emax and 
        saves the constant in self.back and the background subtracted self.signals in self.sqw
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        plt.ion()
        for col in columns:
            inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
            res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), 0)
            yres = np.polyval(res, self.eloss)

            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('signal','constant fit','signal - constant'))
            plt.title('Hit [enter] in the python shell to continue')            
            plt.xlabel('energy loss [eV]')
            plt.ylabel('signal [a.u.]')
            plt.xlim(emin-ewindow,emax+ewindow)
            plt.autoscale(enable=True, axis='y', tight=False)
            
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
            # close the figure to show the next one
        plt.ioff()

    def removeconstav(self,emin,emax,ewindow=100.0):
        """
        fits a constant as background in the range emin-emax and 
        saves the constant in self.back and the background subtracted self.signals in self.sqw
        """

        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
        res  = np.polyfit(self.eloss[inds],self.avsignals[inds], 0)
        yres = np.polyval(res, self.eloss)

        plt.plot(self.eloss,self.avsignals,self.eloss,yres,self.eloss,self.avsignals-yres)
        plt.legend(('signal','constant fit','signal - constant'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(emin-ewindow,emax+ewindow)
        plt.autoscale(enable=True, axis='y', tight=False)
        plt.draw()

        self.sqwav    = self.avsignals - yres
        self.sqwaverr = self.averrors

    def removelinear(self,whichq,emin,emax,ewindow=100.0,stoploop=True):
        """
        fits a linear function as background in the range emin-emax and 
        saves the linear in self.back and the background subtracted self.signals in self.sqw
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        plt.ion()
        for col in columns:
            inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
            res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), 1)
            yres = np.polyval(res, self.eloss)
            
            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('signal','linear fit','signal - linear'))
            plt.title('Hit [enter] in the python shell to continue')            
            plt.xlabel('energy loss [eV]')
            plt.ylabel('signal [a.u.]')
            plt.grid(False)
            plt.xlim(emin-ewindow,emax+ewindow) 
            plt.autoscale(enable=True, axis='y')
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removelinearav(self,emin,emax,ewindow=100.0):
        """
        fits a linear function as background in the range emin-emax from averaged data and saves the result in self.sqwav and self.sqwerrav
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
        res  = np.polyfit(self.eloss[inds],self.avsignals[inds], 1)
        yres = np.polyval(res, self.eloss)
            
        plt.plot(self.eloss,self.avsignals,self.eloss,yres,self.eloss,self.avsignals-yres)
        plt.legend(('signal','linear fit','signal - linear'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.grid(False)
        plt.xlim(emin-ewindow,emax+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav    = self.avsignals - yres
        self.sqwaverr = self.averrors

    def removeLinearAv(self,region1,region2=None,ewindow=100.0,scale=1, view=False):
        """
        fits a linear function as background in the range emin-emax from averaged data and saves the result in self.sqwav and self.sqwerrav
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        range1 = np.where(np.logical_and(self.eloss >= region1[0], self.eloss <= region1[1]))
        if region2:
            range2 = np.where(np.logical_and(self.eloss >= region2[0], self.eloss <= region2[1]))
            region = np.append(range1,range2)
        else:
            region = range1


        res  = np.polyfit(self.eloss[region],self.avsignals[region], 1)
        yres = np.polyval(res, self.eloss)

        newspec = (self.avsignals-yres)*scale
        newerrs = self.averrors*scale

        self.sqwav    = newspec
        self.sqwaverr = newerrs

        self.yres = yres
        self.newspec = newspec

        if view:
                plt.ion()
                plt.cla()
                plt.plot(self.eloss,self.avsignals,self.eloss,yres,self.eloss,newspec,self.eloss,self.avC)
                plt.legend(('signal','linear fit','signal - linear','av. core Compton'))
                plt.xlabel('energy loss [eV]')
                plt.ylabel('signal [a.u.]')
                plt.grid(False)
                plt.xlim(region1[0]-ewindow,region1[-1]+ewindow) 
                plt.autoscale(enable=True, axis='y')
                plt.draw()
                input()


    def removepoly(self,whichq,emin,emax,polyorder=2.0,ewindow=100.0):
        """
        fits a polynomial of order "polyorder" (default is quadratic) as background 
        in the range emin-emax and saves the polynomial in self.back and the background 
        subtracted self.signals in self.sqw
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        plt.ion()
        for col in columns:
            inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
            res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), polyorder)
            yres = np.polyval(res, self.eloss)
            
            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('signal','poly fit','signal - poly'))
            plt.title('Hit [enter] in the python shell to continue')            
            plt.xlabel('energy loss [eV]')
            plt.ylabel('signal [a.u.]')
            plt.xlim(emin-ewindow,emax+ewindow) 
            plt.autoscale(enable=True, axis='y')
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removepolyav(self,polyregion,coreregion,weights=[1,1],scale=1.0,polyorder=2.0,ewindow=100.0,hfcoreshift=0.0):
        """
        fits a polynomial of order "polyorder" (default is quadratic) as background 
        in the range emin-emax to averaged signals, save it to self.sqwav and self.sqwerrav
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        region1 = np.where(np.logical_and(self.eloss >= polyregion[0], self.eloss <= polyregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # shift the HF core edge onset by hfcoreshift
        thecore = np.interp(self.eloss,self.eloss+hfcoreshift,self.avC)

        plt.ion()
        plt.cla()
        print( scale)
        newspec = (self.avsignals)*scale

        #inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
        res  = np.polyfit(self.eloss[region],newspec[region]-thecore[region], polyorder)
        yres = np.polyval(res, self.eloss)

        resspec = (newspec-yres)
        plt.plot(self.eloss,newspec,self.eloss,yres,self.eloss,resspec,self.eloss,thecore)
        plt.legend(('signal','poly fit','signal - poly (scaled)','av. C'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(polyregion[0]-ewindow,coreregion[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav    = resspec
        self.sqwaverr = self.averrors * scale

    def removepolyav1(self,polyregion1,polyregion2=None,polyorder=2.0,weights=[1,1]):
        """
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        region1 = np.where(np.logical_and(self.eloss >= polyregion1[0], self.eloss <= polyregion1[1]))[0]
        if polyregion2:
            region2 = np.where(np.logical_and(self.eloss >= polyregion2[0], self.eloss <= polyregion2[1]))[0]
            region  = np.append(region1*weights[0],region2*weights[1])
        else:
            region = region1

        plt.ion()
        plt.cla()
        newspec = (self.avsignals)

        #inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
        res  = np.polyfit(self.eloss[region],newspec[region], polyorder)
        yres = np.polyval(res, self.eloss)

        resspec = (newspec-yres)
        plt.plot(self.eloss,newspec,self.eloss,yres,self.eloss,resspec)
        plt.legend(('signal','poly fit','signal - poly (scaled)'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.draw()

        self.sqwav    = resspec
        self.sqwaverr = self.averrors


    def removePolyCoreAv(self,polyregion,coreregion,weights=[1,1],guess=[1.0,0.0,0.0],ewindow=100.0):
        """
        fits a polynomial of order "polyorder" (default is quadratic) as background 
        in the range emin-emax to averaged signals, save it to self.sqwav and self.sqwerrav
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        region1 = np.where(np.logical_and(self.eloss >= polyregion[0], self.eloss <= polyregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        funct = lambda a: np.sum( (a[0]*self.avsignals[region] - self.avC[region] - np.polyval(a[1::],self.eloss[region]) )**2.0 )

        res   = optimize.minimize(funct,guess).x
        print( 'the fit results are: ', res)

        yres = np.polyval(res[1::], self.eloss)
        plt.plot(self.eloss,self.avC)    
        plt.plot(self.eloss,self.avsignals*res[0],self.eloss,yres+self.avC,self.eloss,self.avsignals*res[0]-yres)
        plt.legend(('scaled signal','poly fit + core','scaled signal - poly'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(polyregion[0]-ewindow,polyregion[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav    = self.avsignals*res[0] - yres
        self.sqwaverr = self.averrors*res[0]

    def removePolyCoreAv2(self,polyregion,coreregion,weights=[1,1],guess=[1.0,0.0],ewindow=100.0,scale=1.0,hfcoreshift=0.0):
        """
        fits a polynomial of order "polyorder" (default is quadratic) as background 
        in the range emin-emax to averaged signals, save it to self.sqwav and self.sqwerrav
        """
        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        region1 = np.where(np.logical_and(self.eloss >= polyregion[0], self.eloss <= polyregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # try scaling data before the fit (so scale does not have to be a parameter in the fit)
        thespec = self.avsignals * scale

        # shift the HF core edge onset by hfcoreshift
        thecore = np.interp(self.eloss,self.eloss+hfcoreshift,self.avC)

        funct = lambda a: np.sum( (thespec[region] - thecore[region] - np.polyval(a,self.eloss[region]) )**2.0 )

        res   = optimize.minimize(funct,guess).x
        print( 'the fit results are: ', res)

        yres = np.polyval(res, self.eloss)
        plt.plot(self.eloss,thespec,self.eloss,yres+thecore,self.eloss,thespec-yres,self.eloss,thecore)
        plt.legend(('scaled signal','poly fit + core','scaled signal - poly','core profile'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(polyregion[0]-ewindow,polyregion[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav    = thespec - yres
        self.sqwaverr = self.averrors*scale

    def removeconstpcore(self,whichq,constregion,coreregion,weights=[5,1],guess=[0.0, 1.0],ewindow=100.0,stoploop=True):
        """
        fit a const to the preedge and scale data to postedge
        matches the theory profiles
        fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq        

        region1 = np.where(np.logical_and(self.eloss >= constregion[0], self.eloss <= constregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        plt.ion()
        for col in columns:
            # first scale data to same area as core profile in region2
            corenorm = np.trapezoid(self.C[region2,col],self.eloss[region2])
            self.signals[:,col] = self.signals[:,col]/np.trapezoid(self.signals[region2,col],self.eloss[region2])*corenorm

            #fitfct  = lambda a: a[1]*self.signals[region,col] - (np.polyval([a[0]],self.eloss[region])+self.C[region,col])
            #res = optimize.leastsq(fitfct,guess)[0]
            #yres = np.polyval([res[0]],self.eloss)

            c1 = lambda a: -np.sum((a[1]*self.signals[region2,col] - (np.polyval([a[0]],self.eloss[region2])+self.C[region2,col]))**2.0) # post edge should oscillate around HF core profile
            
            fitfct  = lambda a: np.sum( (a[1]*self.signals[region1,col] - (np.polyval([a[0]],self.eloss[region1])+self.C[region1,col])) )
            cons = [c1] #, c2, c3, c4, c5, c6
            res     = optimize.fmin_cobyla(fitfct,guess,cons,maxfun=10000)
            yres = np.polyval([res[0]],self.eloss)

            plt.plot(self.eloss,self.signals[:,col]*res[1],self.eloss,yres+self.C[:,col],self.eloss,self.signals[:,col]*res[1]-yres,self.eloss,self.C[:,col])
            plt.legend(('data','fit','data - constant','core profile'))
            plt.title('Hit [enter] in the python shell to continue')            
            plt.xlabel('energy loss [eV]')
            plt.ylabel('signal [a.u.]')
            plt.xlim(constregion[0]-ewindow,coreregion[1]+ewindow) 
            plt.autoscale(enable=True, axis='y')
            plt.draw()
            # save some data for paper-plot
            #thedata = np.zeros((len(self.eloss),5))
            #thedata[:,0] = self.eloss
            #thedata[:,1] = self.signals[:,col]*res[1]
            #thedata[:,2] = yres+self.C[:,col]
            #thedata[:,3] = self.signals[:,col]*res[1]-yres
            #thedata[:,4] = self.C[:,col]
            #filename = '/home/csahle/Dropbox/tool_paper/figures/analysis/licl_linpcore_background_det' + '%s' % col + '.dat'
            #np.savetxt(filename,thedata)
            self.background[:,col] = yres
            self.sqw[:,col]        = res[1]*self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removelinpcore(self,whichq,linregion,coreregion,weights=[5,1],guess=[0.0, 0.0, 1.0],ewindow=100.0,stoploop=True):
        """
        fit a linear to the preedge and scale data so postedge
        matches the theory profiles
        fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq        

        region1 = np.where(np.logical_and(self.eloss >= linregion[0], self.eloss <= linregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        plt.ion()
        for col in columns:
            # first scale data to same area as core profile in region2
            corenorm = np.trapezoid(self.C[region2,col],self.eloss[region2])
            self.signals[:,col] = self.signals[:,col]/np.trapezoid(self.signals[region2,col],self.eloss[region2])*corenorm

            # then try a minimization
            #fitfct  = lambda a: (a[2]*self.signals[region,col] - (np.polyval(a[0:2],self.eloss[region])+self.C[region,col]))
            #res = optimize.leastsq(fitfct,guess)[0]
            #yres = np.polyval(res[0:2],self.eloss)

            c1 = lambda a: -np.sum((a[2]*self.signals[region2,col] - (np.polyval(a[0:2],self.eloss[region2])+self.C[region2,col]))**2.0) # post edge should oscillate around HF core profile
            c2 = lambda a: a[2] # scaling should not be negative
            fitfct  = lambda a: np.sum( (a[2]*self.signals[region1,col] - (np.polyval(a[0:2],self.eloss[region1])+self.C[region1,col])) )
            cons = [c1,c2] #, c2, c3, c4, c5, c6
            res     = optimize.fmin_cobyla(fitfct,guess,cons)
            yres = np.polyval(res[0:2],self.eloss)

            plt.plot(self.eloss,res[2]*self.signals[:,col],self.eloss,yres+self.C[:,col],self.eloss,res[2]*self.signals[:,col]-yres,self.eloss,self.C[:,col])
            plt.legend(('data','fit','data - linear','core profile'))
            plt.title('Hit [enter] in the python shell to continue')            
            plt.xlabel('energy loss [eV]')
            plt.ylabel('signal [a.u.]')
            plt.xlim(linregion[0]-ewindow,coreregion[1]+ewindow) 
            plt.autoscale(enable=True, axis='y')
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = res[2]*self.signals[:,col] - yres
            if stoploop:            
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removelinpcoreav(self,linregion,coreregion,weights=[5,1],guess=[0.0, 0.0, 1.0],ewindow=100.0):
        """
        fit a linear to the preedge and scale data so postedge
        matches the theory profiles
        fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        """

        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        region1 = np.where(np.logical_and(self.eloss >= linregion[0], self.eloss <= linregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # remove a constant from before the edge
        linguess = guess[0:2]
        fitfct  = lambda a: self.avsignals[region1] - np.polyval(a,self.eloss[region1])    
        res1    = optimize.leastsq(fitfct,linguess)[0]
        back    = np.polyval(res1,self.eloss)

        newspec = self.avsignals - back

        # scale results in a way that it oscillates around the core profile
        coreguess = guess[2]
        fitfct  = lambda a: (a*newspec[region2] - self.avC[region2])
        res2    = optimize.leastsq(fitfct,coreguess)[0]

        # first scale data to same area as core profile in region2
        #corenorm = np.trapezoid(self.avC[region2],self.eloss[region2])
        #self.signalsav = self.avsignals/np.trapezoid(self.avsignals[region2],self.eloss[region2])*corenorm

        #c1 = lambda a: -np.sum((a[2]*self.avsignals[region2] - (np.polyval(a[0:2],self.eloss[region2])+self.avC[region2]))**2.0) # post edge should oscillate around HF core profile
        #c2 = lambda a: a[2] # scaling should not be negative
        #fitfct  = lambda a: np.sum( (a[2]*self.avsignals[region] - (np.polyval(a[0:2],self.eloss[region])+self.avC[region]) )**2.0 )
        #cons = [c1,c2] #, c2, c3, c4, c5, c6
        #res  = optimize.fmin_cobyla(fitfct,guess,cons)
        #yres = np.polyval(res[0:2],self.eloss)

        plt.plot(self.eloss,res2*self.avsignals,self.eloss,self.avsignals,self.eloss,back+self.avC,self.eloss,res2*newspec,self.eloss,self.avC)
        plt.legend(('scaled data','data','fit','data - linear','core profile'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(linregion[0]-ewindow,coreregion[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav  = res2*newspec
        self.sqwaverr = self.averrors*res2

    def removeLinCoreAv(self,linregion,coreregion,weights=[5,1],guess=[0.0, 0.0, 1.0],ewindow=100.0):
        """
        fit a linear to the preedge and scale data so postedge
        matches the theory profiles
        fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        """

        if not np.any(self.avsignals):
            print( 'use averageAnalyzers first to create some averages')
            return

        plt.ion()
        plt.cla()

        region1 = np.where(np.logical_and(self.eloss >= linregion[0], self.eloss <= linregion[1]))
        region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
        region  = np.append(region1*weights[0],region2*weights[1])

        # remove a constant from before the edge
        fitfct  = lambda a: np.sum( (a[2]*self.avsignals[region] - (a[0]*self.eloss[region] + a[1] + self.avC[region]) )**2.0)
        res1    = optimize.minimize(fitfct,guess).x
        back    = np.polyval(res1[0:2],self.eloss)
        print( 'the result of the fit is: ', res1)
        newspec = res1[2]*self.avsignals - back

        plt.plot(self.eloss,res1[2]*self.avsignals,self.eloss,self.avsignals,self.eloss,back+self.avC,self.eloss,newspec,self.eloss,self.avC)
        plt.legend(('scaled data','data','fit','data - linear','core profile'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('signal [a.u.]')
        plt.xlim(linregion[0]-ewindow,coreregion[1]+ewindow) 
        plt.autoscale(enable=True, axis='y')
        plt.draw()

        self.sqwav  = res1[2]*newspec
        self.sqwaverr = self.averrors*res1[2]

    def removepearson(self,whichq,emin,emax,guess=None,stoploop=True):
        """
        guess values: 
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        a[4] = background
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        region = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))[0]
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = np.where( self.signals[guessregion,col] == np.max(self.signals[guessregion,col]) )[0][0]
                guess.append(self.eloss[ind]) # max of signal (in range of prenorm from __init__)
                guess.append(guess[0]*2.0) # twice the position of the peason maximum
                guess.append(1000.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(1.0) # Peak intensity
                guess.append(0.0) # background
            
            fitfct  = lambda a: self.signals[region,col] - pearson7(self.eloss[region],a)
            res     = optimize.leastsq(fitfct,guess)
            yres    = pearson7(self.eloss,res[0])

            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('data','pearson fit','data - pearson'))
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removepearson2(self,whichq,emin,emax,guess=None,stoploop=True):
        """
        guess values: 
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        a[4] = background
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        region = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))[0]
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = np.where( self.signals[guessregion,col] == np.max(self.signals[guessregion,col]) )[0][0]
                guess.append(self.eloss[ind]) # max of signal (in range of prenorm from __init__)
                guess.append(guess[0]*2.0) # twice the position of the peason maximum
                guess.append(1000.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(1.0) # Peak intensity
                guess.append(0.0) # background
                print( guess)
            
            popt, pcov = optimize.curve_fit(pearson7_forcurvefit, self.eloss[region], self.signals[region,col],p0=guess)            
            yres    = pearson7(self.eloss,popt)

            plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
            plt.legend(('data','pearson fit','data - pearson'))
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removePearsonAv(self,region1,region2=None,guess=None,scale=1.0):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        range1 = np.where(np.logical_and(self.eloss >= region1[0], self.eloss <= region1[1]))
        if region2:
            range2 = np.where(np.logical_and(self.eloss >= region2[0], self.eloss <= region2[1]))
            region = np.append(range1,range2)
        else:
            region = range1

        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        plt.cla()

        if not guess: 
            guess = []
            ind = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*1.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[guessregion][ind]) # Peak intensity
            guess.append(0.0) # linear slope
            guess.append(0.0) # linear background

        # fit a pearson to the whole region
        res1  = optimize.curve_fit(pearson7_linear_forcurvefit, self.eloss[region], self.avsignals[region],p0=guess)[0]
        yres1 = pearson7_zeroback(self.eloss,res1[0:4]) + np.polyval(res1[4:6],self.eloss)
        print( 'the fitting results are: ', res1)

        newspec = (self.avsignals-yres1)*scale
        plt.plot(self.eloss,self.avsignals,self.eloss,yres1,self.eloss,newspec,self.eloss,self.avC)
        plt.legend(('data','pearson fit','data - pearson'))
        plt.draw()

        self.sqwav = newspec
        self.sqwaverr = self.averrors * scale

    def removePearsonAv2(self,region1,region2=None,weights=[1,1],guess=None,ewindow=100.0,scale=1.0,hfcoreshift=0.0):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        range1 = np.where(np.logical_and(self.eloss >= region1[0], self.eloss <= region1[1]))
        if region2:
            range2 = np.where(np.logical_and(self.eloss >= region2[0], self.eloss <= region2[1]))
            region = np.append(range1*weights[0],range2*weights[1])
        else:
            region = range1

        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        plt.cla()

        if not guess: 
            guess = []
            ind = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*1.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[guessregion][ind]) # Peak intensity
            guess.append(0.0) # linear slope
            guess.append(0.0) # linear background

        # try scaling data before the fit (so scale does not have to be a parameter in the fit)
        thespec = self.avsignals * scale

        # shift the HF core edge onset by hfcoreshift
        thecore = np.interp(self.eloss,self.eloss+hfcoreshift,self.avC)

        # fit a pearson to the whole region
        fitfct  = lambda a: np.sum( (thespec[region] - pearson7_zeroback(self.eloss[region],a[0:4]) - np.polyval(a[4:6],self.eloss[region]) - thecore[region])**2.0 )
        res = optimize.minimize(fitfct,guess).x
        #res1  = optimize.curve_fit(pearson7_linear_forcurvefit, self.eloss[region], thespec[region],p0=guess)[0]
        #yres1 = pearson7_zeroback(self.eloss,res1[0:4]) + np.polyval(res1[4:6],self.eloss)
        print( 'the fitting results are: ', res)

        yres = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)

        plt.plot(self.eloss,thespec,self.eloss,yres,self.eloss,thespec-yres,self.eloss,thecore)
        plt.legend(('data','pearson fit','data - pearson','core'))
        plt.draw()

        self.sqwav = thespec-yres
        self.sqwaverr = self.averrors * scale

    def removecoreppearson(self,whichq,pearsonrange,postrange,weights=[2,1],guess=None,stoploop=True):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq        

        region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))
        region2 = np.where(np.logical_and(self.eloss >= postrange[0], self.eloss <= postrange[1]))
        region  = np.append(region1*weights[0],region2*weights[1])
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        print( len(self.eloss[region]))

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = self.signals[guessregion,col].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
                guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
                guess.append(guess[0]*1.0) # once the position of the peason maximum
                guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(self.signals[guessregion,col][ind]) # Peak intensity
                guess.append(0.0) # linear slope
                guess.append(0.0) # linear background
                guess.append(1.0) # scaling factor for exp. data

            # some sensible boundary conditions for the fit:
            c1 = lambda a: a[1]*np.absolute(2e2  - a[1]) # FWHM should not be bigger than 200 eV and positive
            c2 = lambda a: a[2] # shape should not be negative
            c3 = lambda a: a[3] # peak intensity should not be negative
            c4 = lambda a: np.absolute(5e-1 - a[4]) # slope for linear background should be small
            c5 = lambda a: a[3] - a[5] # offset for linear should be smaller than maximum of pearson
            c6 = lambda a: a[6]*np.absolute(1e10 - a[6]) # scaling factor for the data should not be negative
            c7 = lambda a: np.sum( (a[5]*self.signals[region2,col] - pearson7_zeroback(self.eloss[region2],a[0:5]) - self.C[region2,col])**2.0 )

            fitfct  = lambda a: np.sum( (a[6]*self.signals[region,col] - pearson7_zeroback(self.eloss[region],a[0:4]) - np.polyval(a[4:6],self.eloss[region]) - self.C[region,col])**2.0 )
            cons = [c7] #[c1, c2, c3, c4, c5, c6, c7]
            res     = optimize.fmin_cobyla(fitfct,guess,cons)
            print( res)
            yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss) + self.C[:,col]

            plt.plot(self.eloss,self.signals[:,col]*res[6],self.eloss,yres,self.eloss,self.signals[:,col]*res[6]-yres,self.eloss,self.C[:,col])
            plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear +core)','core'))
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removeCorePearsonAv(self,fitrange,guess=None):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        inds = np.where(np.logical_and(self.eloss >= fitrange[0], self.eloss <= fitrange[1]))
        plt.ion()
        plt.cla()

        if not guess: 
            guess = []
            ind = self.avsignals[inds].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[inds][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*1.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[inds][ind]) # Peak intensity
            guess.append(0.0) # linear slope
            guess.append(0.0) # linear background
            guess.append(1.0) # scaling factor for exp. data

        fitfct  = lambda a: a[6]*self.avsignals[inds] - (pearson7_zeroback(self.eloss[inds],a[0:4]) + np.polyval(a[4:6],self.eloss[inds]) + self.avC[inds])
        res     = optimize.leastsq(fitfct,guess)[0]
        yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)

        print( 'the current fit-parameters are: ' + str(res) + ', try using these as guess parameters in a more refined fit!')

        plt.plot(self.eloss,self.avsignals*res[6],self.eloss,yres+self.avC,self.eloss,self.avsignals*res[6]-yres,self.eloss,self.avC)
        plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear +core)','core'))
        plt.draw()

        self.sqwav = self.avsignals*res[6] - yres
        self.sqwaverr = self.averrors*res[6]

    def removeCorePearsonAv2(self,pearsonrange,corerange,weights=[2,1],guess=None):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))
        region2 = np.where(np.logical_and(self.eloss >= corerange[0], self.eloss <= corerange[1]))
        region  = np.append(region1*weights[0],region2*weights[1])
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        plt.cla()

        # scale the data approximately to get a better match in edge jumb btw. data and compton profiles
        #try:
        #theojump = np.abs(self.avC[region2[0][0]] - self.avC[region1[0][-1]])
        #datajump = np.abs(self.avsignals[region2[0][0]] - self.avsignals[region1[0][-1]])
        #prescaling = datajump/theojump
        #except:
        #    prescaling = 1.0
        #self.avsignals *= prescaling

        if not guess: 
            guess = []
            ind = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*1.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[guessregion][ind]) # Peak intensity
            guess.append(0.0) # linear slope
            guess.append(0.0) # linear background
            guess.append(1.0) # scaling factor for exp. data

        # some sensible boundary conditions for the fit:
        c1 = lambda a: a[1]*np.absolute(2e2 - a[1]) # FWHM should not be bigger than 200 eV and positive
        c2 = lambda a: a[2] # shape should not be negative
        c3 = lambda a: a[3] # peak intensity should not be negative
        c4 = lambda a: np.absolute(5e-1 - a[4]) # slope for linear background should be small
        c5 = lambda a: a[3] - a[5] # offset for linear should be smaller than maximum of pearson
        c6 = lambda a: a[6] # scaling factor for the data should not be negative
        c7 = lambda a: np.sum( (a[6]*self.avsignals[region2] - pearson7_zeroback(self.eloss[region2],a[0:4]) - self.avC[region2] - np.polyval(a[4:6],self.eloss[region2]))**2.0 )

        fitfct  = lambda a: np.sum( (a[6]*self.avsignals[region1] - pearson7_zeroback(self.eloss[region1],a[0:4]) - np.polyval(a[4:6],self.eloss[region1]) - self.avC[region1])**2.0 ) + np.sum( (a[6]*self.avsignals[region2] - pearson7_zeroback(self.eloss[region2],a[0:4]) - np.polyval(a[4:6],self.eloss[region2]) - self.avC[region2])**2.0 )
        cons = []#[c1, c2, c3, c4, c5, c6, c7]
        res     = optimize.minimize(fitfct,guess)
        yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)

        plt.plot(self.eloss,self.avsignals*res[6],self.eloss,yres+self.avC,self.eloss,self.avsignals*res[6]-yres,self.eloss,self.avC)
        plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear)','core'))
        plt.draw()

        self.sqwav = self.avsignals*res[6] - yres
        self.sqwaverr = self.averrors*res[6]

    def removeCorePearsonAv3(self,pearsonrange,corerange,weights=[1,1],guess=None):
        """
        weights must be integers!
        
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))
        region2 = np.where(np.logical_and(self.eloss >= corerange[0], self.eloss <= corerange[1]))
        region  = np.append(region1*weights[0],region2*weights[1])
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        plt.cla()

        if not guess: 
            guess = []
            ind = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
            guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
            guess.append(guess[0]*1.0) # once the position of the peason maximum
            guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
            guess.append(self.avsignals[guessregion][ind]) # Peak intensity
            guess.append(0.0) # linear slope
            guess.append(0.0) # linear background

        # fit a pearson to the whole region
        res1  = optimize.curve_fit(pearson7_linear_forcurvefit, self.eloss[region], self.avsignals[region],p0=guess)[0]
        yres1 = pearson7_zeroback(self.eloss,res1[0:4]) + np.polyval(res1[4:6],self.eloss)

        plt.plot(self.eloss,self.avsignals,self.eloss,yres1+self.avC,self.eloss,self.avsignals-yres1,self.eloss,self.avC)

        # estimate a scaling factor from area in the coreregion
        scale = np.trapezoid(self.avC[region2],self.eloss[region2])/np.abs(np.trapezoid(self.avsignals[region2]-yres1[region2],self.eloss[region2]))
        print( 'trapz of core is ',np.trapezoid(self.avC[region2],self.eloss[region2]))
        print( 'trapz of data is ',np.trapezoid(self.avsignals[region2]-yres1[region2],self.eloss[region2]))
        print( 'scale is ', scale)
        newspec = self.avsignals * scale
        newerrs = self.averrors * scale
        # fit scaled averaged data again (fit scaling with this too)
        guess = res1
        guess = np.append(guess,1.0) # scaling factor for exp. data

        res2  = optimize.curve_fit(pearson7_linear_scaling_forcurvefit, self.eloss[region], newspec[region],p0=guess)[0]
        yres2 = pearson7_zeroback(self.eloss,res2[0:4]) + np.polyval(res2[4:6],self.eloss)

        print(  res2)
        plt.figure()
        plt.plot(self.eloss,newspec*res2[6],self.eloss,yres2+self.avC,self.eloss,newspec*res2[6]-yres2,self.eloss,self.avC)
        plt.legend(('scaled data','pearson + linear + core','scaled data - (pearson + linear)','core'))
        plt.draw()

        self.sqwav = newspec*res2[6] - yres2
        self.sqwaverr = newerrs*res2[6]

    def removecoreppearson2(self,whichq,pearsonrange,postrange,guess=None,stoploop=True):
        """
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        linear:
        a[4] = linear slope
        a[5] = linear background/offset
        data: 
        a[6] = scaling factor
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq        

        region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))[0]
        region2 = np.where(np.logical_and(self.eloss >= postrange[0], self.eloss <= postrange[1]))[0]
        region = np.append(region1,region2)
        
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = self.signals[guessregion,col].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
                guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
                guess.append(guess[0]*1.0) # once the position of the peason maximum
                guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(self.signals[guessregion,col][ind]) # Peak intensity
                guess.append(0.0) # ax
                guess.append(0.0) # b
                guess.append(1.0) # scaling factor for exp. data
                print( guess)

            # boundary conditions for the fit: let the post-edge region oscilate around the HF core profile
            c1 = lambda a: np.sum( (a[6]*self.signals[region2,col] - pearson7_zeroback(self.eloss[region2],a[0:4]) - self.C[region2,col])**2.0 )
            c2 = lambda a: a[3]
            c3 = lambda a: a[1]
            c4 = lambda a: a[5]

            #fitfct  = lambda a: np.sum( (a[5]*self.signals[region1,col] - pearson7_zeroback(self.eloss[region1],a[0:4]) - self.C[region1,col])**2.0 )
            fitfct  = lambda a: np.sum( ( a[6]*self.signals[region,col] - pearson7_zeroback(self.eloss[region],a[0:4]) - np.polyval(a[4:6],self.eloss[region]) - self.C[region,col])**2.0 )
            cons = [c2,c3,c4] #, c2, c3, c4, c5, c6
            res     = optimize.fmin_cobyla(fitfct,guess,cons)
            #res     = optimize.fmin(fitfct,guess)

            yres    = pearson7_zeroback(self.eloss,res[0:4])

            plt.plot(self.eloss,res[6]*self.signals[:,col], self.eloss,pearson7_zeroback(self.eloss[:],res[0:4]), self.eloss,np.polyval(res[4:6],self.eloss[:]) )

            #plt.plot(self.eloss,self.signals[:,col]*res[6],self.eloss,yres,self.eloss,self.signals[:,col]*res[6]-yres ,self.eloss,self.C[:,col])

            #plt.legend(('scaled data','pearson','data - pearson','core'))
            plt.draw()

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def removecoreppearson3(self,whichq,pearsonrange,postrange,guess=None,stoploop=True):
        """
        guess values: 
        pearson (always zero background):
        a[0] = Peak position
        a[1] = FWHM
        a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
        a[3] = Peak intensity
        a[4] = pearson offset
        linear:
        a[5] = linear slope
        a[6] = linear background/offset
        data: 
        a[7] = scaling factor
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq        

        region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))[0]
        region2 = np.where(np.logical_and(self.eloss >= postrange[0], self.eloss <= postrange[1]))[0]
        
        guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

        plt.ion()
        for col in columns:
            if not guess: 
                guess = []
                ind = self.signals[guessregion,col].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
                guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
                guess.append(guess[0]*1.0) # once the position of the peason maximum
                guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
                guess.append(self.signals[guessregion,col][ind]) # Peak intensity
                guess.append(0.0) # pearson offset
                guess.append(0.0) # scaling factor for exp. data
                guess.append(0.0)
                guess.append(1.0)
                print( 'guessing start values: ', guess)

            # formulate some contraints:
            def constr1(a):
                return a[1] # let FWHM be positive
            def constr2(a):
                return a[3] # let the peak intensity be positive
            def constr3(a):
                return a[7] # scaling should be positive

            # initial fit:
            fitfct  = lambda a: np.sum( (a[5]*self.signals[region1,col] - pearson7_zeroback(self.eloss[region1],a[0:5]))**2.0 )
            res     = optimize.fmin_cobyla(fitfct,guess,[constr1, constr2, constr3])

            # try again and again to optimize until post-edge region oscilates around the HF core profile
            etol = 1e-4
            while np.sum( (res[5]*self.signals[region2,col] - pearson7_zeroback(self.eloss[region2],res[0:5]) - self.C[region2,col])**2.0 )>etol:
                res = optimize.fmin_cobyla(fitfct,res,[constr1, constr2, constr3])

            yres    = pearson7_zeroback(self.eloss,res[0:5])
            if stoploop:
                plt.plot(self.eloss,self.signals[:,col]*res[7],self.eloss,yres,self.eloss,self.signals[:,col]*res[7]-yres ,self.eloss,self.C[:,col])
                plt.legend(('scaled data','pearson','data - pearson','core'))
                plt.draw()
                _ = input("Press [enter] to continue.") # wait for input from the user
                plt.close()    # close the figure to show the next one

            self.background[:,col] = yres
            self.sqw[:,col]        = self.signals[:,col] - yres

    def remquickval(self,whichq,corefitrange,interpolrange,convwidth,stoploop=True):
        """
        quick and dirty way of valence profile extraction from a single spectrum. 
        works if the edge rides on the tail of the valence profile at high q. the HF core
        profile is fitted to the spectrum in the 'corefitrange' and the resulting valence 
        profile is cut out and interpolated over in the 'interpolrange'. finally the valence 
        profile is smoothed by convolution with a gaussian of FWHM 'convwidth' and subtracted
        from the original data such that the resulting S(q,w) oscillates around the HF
        core profile.
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        fitrange = np.where(np.logical_and(self.eloss >= corefitrange[0], self.eloss <= corefitrange[1]))[0]

        interprange1 = np.where(self.eloss<=interpolrange[0])[0]
        interprange2 = np.where(self.eloss>=interpolrange[1])[0]
        interprange  = np.append(interprange1,interprange2)

        plt.ion()
        for col in columns:
            fitfct   = lambda a: np.sum((self.signals[fitrange,col] - a*self.C[fitrange,col])**2.0)
            constr   = lambda a: a # scaling factor for the HF core profile should not be negative

            res = optimize.fmin_cobyla(fitfct,[1.0],cons=constr)[0]

            # subtract the HF core compton profile and interpolate through the edge
            f       = interpolate.interp1d(self.eloss[interprange],self.signals[interprange,col]-res*self.C[interprange,col], bounds_error=False, fill_value=0.0)
            valdata = f(self.eloss)
            valdata = convg(self.eloss,valdata,convwidth)

            subdata = self.signals[:,col] - valdata

            plt.plot(self.eloss,self.signals[:,col],self.eloss,self.C[:,col]*res,self.eloss,valdata,self.eloss,subdata)
            plt.legend(('data','scaled core compton','estimated valence','extracted data'))            
            plt.draw()

            self.valence[:,col] = valdata
            self.sqw[:,col]     = subdata/res # scale the extracted data back to fit the HF core profile

            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()    # close the figure to show the next one
        plt.ioff()

    def extractval(self,whichq,mirror=False,linrange1=None,linrange2=None):
        """
        extracts a valence profile from q-value(s) given in whichq by first fitting 
        the core HF profile to the data at places linrange1 and linrange2 (one, two, 
        or no ranges can be given), then subtracting the HF profile from the data. the
        resulting valence profile in the near-edge region can be replaced by a pearson 
        function (default) or by mirroring the negative side of the valence profile 
        (mirror=True). if mirror is set to False, also the asymmetry is fitted.
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        # set pz scale of first given q-val as 'master' grid
        absolutepz = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[columns[0]])[0]
        
        plt.cla()
        plt.ion()
        for col in columns: 
            # set the pz scale for each q
            pz = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[col])[0]

            if linrange1 and linrange2:
                range1   = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
                range2   = np.where(np.logical_and(self.eloss>=linrange2[0],self.eloss<=linrange2[1]))[0]
                linrange = np.append(range1,range2)
            elif linrange1:
                linrange = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
            else: 
                linrange = np.where(0.1*self.C[:,col] > self.V[:,col])[0]

            # simple minimization to subtract a linear from the data to fit ontop of the Compton profile:
            fitfct = lambda a: (self.signals[linrange,col] - np.polyval([a[0],a[1]],self.eloss[linrange]) ) - self.J[linrange,col]
            res    = optimize.leastsq(fitfct,[0.0,0.0])[0]

            # raw valence (when later extracted from the data, also a linear shoud be accounted for)
            val = self.signals[:,col] - np.polyval(res,self.eloss) - self.C[:,col]

            if mirror: # just replace the edgepart of the valence profile by the other half of the profile
                mirrorval = np.append(val[pz<=0.0],np.flipud(val[pz<=0]))
                mirrorpz  = np.append(pz[pz<=0],np.flipud(pz[pz<=0]*-1))
                order = np.argsort(mirrorpz)

                f = interpolate.interp1d(mirrorpz[order],mirrorval[order],bounds_error=False, fill_value=0.0)
                extractedval = f(absolutepz)
                
                plt.plot(absolutepz,val,absolutepz,extractedval)
                plt.legend(['exp. S(q,w) - HF core profile','mirrored extracted valence profile'])
                plt.xlabel('pz [a.u.]')
                plt.ylabel('S(q,w) [1/eV]')
                plt.draw()
                _ = input("Press [enter] to continue.") # wait for input from the user
                plt.close()

                self.valence[:,col] = extractedval
                self.valasymmetry = np.zeros_like(self.valence)

            else: # fit pearson to replace near edge part
                print ('select a point above which the valence profile should be replaced by a Pearson function')
                plt.plot(absolutepz,val)
                plt.ylim(np.amin(val)-np.amin(val)*1.1,val[absolutepz.flat[np.abs(absolutepz - 0.0).argmin()]]*1.5) 
                xyval       = pylab.ginput(1)[0]
                edgeregion  = np.where(pz < xyval[0])[0]
                start_param = [pz[val==np.amax(val[edgeregion])][0], 4.0, 1.0, np.amax(val[edgeregion]), 0.0 ]
                fitfct      = lambda a: val[edgeregion] - pearson7(pz[edgeregion],a)
                param = optimize.leastsq(fitfct,start_param)[0]
                # param   = optimize.curve_fit(pearson7_forcurvefit, pz[theregion], val[theregion],p0=start_param)[0]
                pearson = pearson7(pz,param)
                val[pz>xyval[0]] = pearson[pz>xyval[0]]
                extractedval = val;
                plt.close()

                # try fitting the valence asymmetry
                print( 'trying to extract valence asymmetry!')
                pzp = absolutepz[absolutepz >= 0.0]
                pzm = absolutepz[absolutepz < 0.0]
                jvalp  = extractedval[absolutepz >=0.0 ]
                f = interpolate.UnivariateSpline(-pzm,extractedval[absolutepz<0.0])
                jvalm  = f(pzp)
                fitfct = lambda a: jvalp-jvalm - a[0]*(np.tanh(pzp/a[1])*np.exp(-(pzp/np.absolute(a[2]))**4.0))
                res    = optimize.leastsq(fitfct,[0.0,1.0,1.0])[0]
                asym   = (res[0]*(np.tanh(pz/res[1])*np.exp(-(pz/np.absolute(res[2]))**4.0)))/2.0
                plt.plot(absolutepz,extractedval,absolutepz,asym,absolutepz,extractedval+asym)
                plt.legend(['extracted valence profile','fitted valence asymmetry','asymmetry corrected valence profile'])
                plt.draw()
                self.valencepz[:,col]      = extractedval-asym
                self.valasymmetrypz[:,col] = asym

        self.pzscale = absolutepz
        plt.ioff()

    def extractval_test(self,whichq,mirror=False,linrange1=None,linrange2=None):
        """
        extracts a valence profile from q-value(s) given in whichq by first fitting 
        the core HF profile to the data at places linrange1 and linrange2 (one, two, 
        or no ranges can be given), then subtracting the HF profile from the data. the
        resulting valence profile in the near-edge region can be replaced by a pearson 
        function (default) or by mirroring the negative side of the valence profile 
        (mirror=True). if mirror is set to False, also the asymmetry is fitted.
        """
        
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        plt.cla()
        plt.ion()

        for col in columns: 
            # set the pz scale for each q
            pz_dmy = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[col])[0]

            # find the regions in which to fit a linear function
            if linrange1 and linrange2:
                range1   = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
                range2   = np.where(np.logical_and(self.eloss>=linrange2[0],self.eloss<=linrange2[1]))[0]
                linrange = np.append(range1,range2)
            elif linrange1:
                linrange = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
            else: 
                linrange = np.where(0.1*self.C[:,col] > self.V[:,col])[0]

            # simple minimization to subtract a linear from the data to fit ontop of the HF Compton profile:
            fitfct = lambda a: (self.signals[linrange,col] - np.polyval([a[0],a[1]],self.eloss[linrange]) ) - a[2]*self.J[linrange,col]
            res    = optimize.leastsq(fitfct,[0.0,0.0,1.0])[0]
            self.background[:,col] = np.polyval(res[0:2],self.eloss) # save the linear 

            # raw valence (when later extracted from the data, also a linear shoud be accounted for)
            val = self.signals[:,col] - np.polyval(res[0:2],self.eloss) - res[2]*self.C[:,col]
            
            #plt.plot(self.eloss,self.signals[:,col],self.eloss,self.C[:,col]*res[2]+np.polyval(res[0:2],self.eloss))
            if mirror: # just replace the edgepart of the valence profile by the other half of the profile
                mirrorval = np.append(val[pz_dmy<=0.0],np.flipud(val[pz_dmy<=0]))
                mirrorpz  = np.append(pz_dmy[pz_dmy<=0],np.flipud(pz_dmy[pz_dmy<=0]*-1))
                order = np.argsort(mirrorpz)

                f = interpolate.interp1d(mirrorpz[order],mirrorval[order],bounds_error=False, fill_value=0.0)
                extractedval = f(pz_dmy)
                
                plt.plot(pz_dmy,val,pz_dmy,extractedval)
                plt.legend(['exp. S(q,w) - HF core profile','mirrored extracted valence profile'])
                plt.xlabel('pz [a.u.]')
                plt.ylabel('S(q,w) [1/eV]')
                plt.draw()
                _ = input("Press [enter] to continue.") # wait for input from the user
                plt.close()

                self.valence[:,col] = extractedval
                self.valasymmetry = np.zeros_like(self.valence)

            else: # fit pearson to replace near edge part
                print ('select a point above which the valence profile should be replaced by a Pearson function')
                plt.plot(pz_dmy,val)
                plt.ylim((np.amin(val[pz_dmy<2.0])*0.9,np.amax(val[pz_dmy<2.0])*1.1)) # np.amin(val)-np.amin(val)*1.1,val[pz_dmy.flat[np.abs(pz_dmy - 0.0).argmin()]]*1.5) 
                xyval       = pylab.ginput(1)[0]
                edgeregion  = np.where(pz_dmy < xyval[0])[0]
                start_param = [pz_dmy[val==np.amax(val[edgeregion])][0], 4.0, 1.0, np.amax(val[edgeregion]), 0.0 ]
                fitfct      = lambda a: val[edgeregion] - pearson7(pz_dmy[edgeregion],a)
                param = optimize.leastsq(fitfct,start_param)[0]
                # param   = optimize.curve_fit(pearson7_forcurvefit, pz_dmy[theregion], val[theregion],p0=start_param)[0]
                pearson = pearson7(pz_dmy,param)
                val[pz_dmy>xyval[0]] = pearson[pz_dmy>xyval[0]]
                extractedval = val;
                plt.close()

                # try fitting the valence asymmetry
                print( 'trying to extract valence asymmetry!')
                pzp = pz_dmy[pz_dmy >= 0.0]
                pzm = pz_dmy[pz_dmy < 0.0]
                jvalp  = extractedval[pz_dmy >=0.0 ]
                f = interpolate.interp1d(-pzm,extractedval[pz_dmy<0.0],bounds_error=False, fill_value=0.0)
                jvalm  = f(pzp)
                fitfct = lambda a: jvalp-jvalm - a[0]*(np.tanh(pzp/a[1])*np.exp(-(pzp/np.absolute(a[2]))**4.0))
                res    = optimize.leastsq(fitfct,[0.0,1.0,1.0])[0]
                print( res)
                asym   = -(res[0]*(np.tanh(pz_dmy/res[1])*np.exp(-(pz_dmy/np.absolute(res[2]))**4.0)))/2.0
                plt.plot(pz_dmy,extractedval,pz_dmy,asym,pz_dmy,extractedval+asym)
                #plt.plot(pz_dmy[pz_dmy<0],extractedval[pz_dmy<0]+asym[pz_dmy<0],-pz_dmy[pz_dmy>=0],extractedval[pz_dmy>=0]+asym[pz_dmy>=0])
                plt.legend(['extracted valence profile','fitted valence asymmetry','asymmetry corrected valence profile'])
                plt.draw()
                #self.valence[:,col]      = extractedval-asym
                #self.valasymmetry[:,col] = asym

                #print pz_dmy[0], pz_dmy[-1], absolutepz[0],absolutepz[-1]
                f = interpolate.interp1d(np.flipud(pz_dmy),np.flipud(extractedval-asym), kind='cubic',bounds_error=False, fill_value=0.0)
                #f = interpolate.UnivariateSpline(pz_dmy,extractedval-asym,k=3,s=10)
                absval = f(np.flipud(self.pzscale))
                self.valencepz[:,col] = np.flipud(absval)
                f = interpolate.interp1d(np.flipud(pz_dmy),np.flipud(asym), kind='cubic',bounds_error=False, fill_value=0.0)
                absasym = f(np.flipud(self.pzscale))
                self.valasymmetrypz[:,col] = np.flipud(absasym)
                #print absval
                #plt.cla()                
                #plt.plot(absolutepz,absval)
                #plt.cla()
                #plt.plot(absolutepz,np.interp(absolutepz,pz_dmy,extractedval-asym),pz_dmy,extractedval-asym)
                #break

            #self.pzscale[:,col] = absolutepz #pz_dmy
        plt.ioff()

    def getallvalprof(self,whichq,smoothgval=0.0,stoploop=True):
        """
        takes the extracted valence profile extracted from q-value whichq
        and transforms them onto the other q-values 
        whichq     = column from which the valence profile was extracted
        smoothgval = FWHM used for gaussian smoothing (default is 0.0, i.e. no smoothing)
        stoploop   = boolean, plots each result if set to True
        """
        newenergy  = np.zeros((len(self.pzscale),len(self.tth)))
        newvalence = np.zeros((len(self.eloss),len(self.tth)))
        newasym    = np.zeros((len(self.eloss),len(self.tth)))
        plt.ion()
        for n in range(len(self.tth)):
            newenergy[:,n]  = (pz2e1(self.E0,self.pzscale,self.tth[n]) - self.E0)*1e3 # each energy scale in [eV]
            f               = interpolate.interp1d(newenergy[:,n],self.valencepz[:,whichq],bounds_error=False, fill_value=0.0)
            newvalence[:,n] = f(self.eloss)/self.qvals[:,n]
            f               = interpolate.interp1d(newenergy[:,n],self.valasymmetrypz[:,whichq],bounds_error=False, fill_value=0.0)
            newasym[:,n]    = f(self.eloss)/self.qvals[:,n]

            plt.plot(self.eloss,newvalence[:,n],self.eloss,newasym[:,n])
            plt.draw()
            if stoploop:
                _ = input("Press [enter] to continue.") # wait for input from the user
            plt.close()

        if smoothgval > 0.0:
            for n in range(len(self.tth)):
                self.valence[:,n] = convg(self.eloss,newvalence[:,n],smoothgval) + newasym[:,n]
            self.valasymmetry = newasym
        else:
            self.valence = newvalence + newasym
            self.valasymmetry = newasym
        plt.ioff()

    def remvalenceprof(self,whichq,eoffset=0.0):
        """
        removes extracted valence profile from q-values given in list whichq by fitting the scaled
        valence profile plus a linear function plus the HF core compton profile to the data
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        inds = np.where(self.eloss>=self.prenormrange[0])[0]

        plt.ion()
        for col in columns:
            f = interpolate.interp1d(self.eloss+eoffset,self.valence[:,col],bounds_error=False, fill_value=0.0)
            #self.valence[:,col]  = f(self.eloss)
            thevalence = f(self.eloss)
            fitfct = lambda a: self.signals[inds,col]-self.C[inds,col]-a[0]*thevalence[inds]-np.polyval([a[1],a[2]],self.eloss[inds])
            res    = optimize.leastsq(fitfct,[6.0,0.0,0.0])[0]
            plt.plot(self.eloss,self.signals[:,col])
            plt.plot(self.eloss,self.C[:,col]+res[0]*thevalence)
            plt.plot(self.eloss,np.polyval(res[1:3],self.eloss))
            plt.plot(self.eloss,self.signals[:,col]-res[0]*thevalence-np.polyval([res[1],res[2]],self.eloss),self.eloss,self.C[:,col])
            plt.draw()
            self.sqw = self.signals[:,col]-res[0]*thevalence-np.polyval([res[1],res[2]],self.eloss)
            self.valence[:,col] = thevalence
        plt.ioff()

    def remvalenceprof_test(self,whichq,eoffset=0.0):
        """
        removes extracted valence profile from q-values given in list whichq by fitting the scaled
        valence profile plus a linear function plus the HF core compton profile to the data
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else:
            columns = whichq

        inds = np.where(self.eloss>=self.prenormrange[0])[0]

        plt.ion()
        for col in columns:
            f = interpolate.interp1d(self.eloss+eoffset,self.valence[:,col],bounds_error=False, fill_value=0.0)
            #self.valence[:,col]  = f(self.eloss)
            thevalence = f(self.eloss)

            fitfct = lambda a: self.signals[inds,col]-self.C[inds,col]-a[0]*np.interp(self.eloss,self.eloss+a[1],thevalence)[inds]-np.polyval([a[2],a[3]],self.eloss[inds])
            res    = optimize.leastsq(fitfct,[1.0,0.0,0.0,0.0])[0]
            plt.plot(self.eloss,self.signals[:,col])
            plt.plot(self.eloss,self.C[:,col]+res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)+np.polyval(res[2:4],self.eloss))
            plt.plot(self.eloss,np.polyval(res[2:4],self.eloss))
            plt.plot(self.eloss,self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss),self.eloss,self.C[:,col])
            plt.draw()
            self.sqw[:,col] = self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss)
            self.valence[:,col] = res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence) # reassign shifted and scaled valence profile
            #self.valasymmetry *= res[0] # reassign scaled valence asymmetry for plotting
            self.background[:,col] = np.polyval(res[2:4],self.eloss)
            # save some ascii files for paper
            thedata = np.zeros((len(self.eloss),6))
            thedata[:,0] = self.eloss
            thedata[:,1] = self.signals[:,col]
            thedata[:,2] = res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)+np.polyval(res[2:4],self.eloss)
            thedata[:,3] = np.polyval(res[2:4],self.eloss)
            thedata[:,4] = self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss)
            thedata[:,5] = self.C[:,col]
            filename = '/home/christoph/Dropbox/tool_paper/figures/analysis/val_extraction_det' + '%s' % col + '.dat'
            #np.savetxt(filename,thedata) 
        plt.ioff()

    def averageqs(self,whichq,errorweighing=True):
        """
        averages S(q,w) over the q-values given
        whichq        = list of q-values over which to average (index starts at zero)
        errorweighing = boolean, weighs sum by errors if set to True
        """
        if not isinstance(whichq,list):
            columns = []
            columns.append(whichq)
        else: 
            columns = whichq

        # build the matricies
        av    = np.zeros((len(self.eloss),len(whichq)))
        averr = np.zeros((len(self.eloss),len(whichq)))
        for n in range(len(columns)):
            # find data points with error = 0.0 and replace by 1.0
            inds = np.where(self.errors[:,columns[n]] == 0.0)[0]
            for ind in inds:
                self.errors[ind,columns[n]] = 1.0
            # arrange the desired columns into a matrix
            av[:,n]    = self.sqw[:,columns[n]]
            averr[:,n] = self.errors[:,columns[n]]
        # sum things up
        if errorweighing:
            self.sqwav    = np.sum( av/averr**2.0 ,axis=1)/( np.sum(1.0/averr**2.0,axis=1))
            self.sqwaverr = np.sqrt( 1.0/np.sum(1.0/(averr)**2.0,axis=1) )

        else: 
            self.sqwav    = np.sum(av,axis=1)
            self.sqwaverr = np.sqrt(np.sum(np.absolute(self.errors)**2.0,axis=1)) # check this again


        def save_state_hdf5(self,  filename, groupname, comment ="" ):
            import h5py
            h5 = h5py.File(filename,"a")
            if(  groupname  in list(h5.keys()) ):
                del h5[groupname]

            h5.require_group(groupname)
            h5group =  h5[groupname]
            for key in ["eloss" ,"sqwav" ,"sqwaverr" ,"avsignals","avC" ,"yres" ,"newspec" ]: 
                if hasattr( self, key ) :
                    data= getattr(self,key)
                    h5group[key]  =  data

            h5group["comment"]  = comment
            h5.flush()
            h5.close()

    def savetxtsqwav(self,filename, emin=None, emax=None, normrange=None):
        """
        save the S(q,w) into a filename (energy loss, sqw, error), save only part of the spectrum, if emin and emax are given, normalize to         area using np.trapezoid if normrange is given and a list of length 2 
        """
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
            assert type(normrange) is list and len(normrange) is 2, "normrange has to be a list of length two!"
            inds = np.where(np.logical_and(data[:,0]>=normrange[0],data[:,0]<=normrange[1]))[0]
            norm = trapz(data[inds,1],data[inds,0])
            data[:,1] /= norm
            data[:,2] /= norm
        np.savetxt(filename,data)

    def extract_highq_valence_profile(
        self,
        highq,
        linrange1=None,
        linrange2=None,
        mirror=False,
        pearson_cut_pz=None,
        fit_total_J=True,
        view=True
    ):
        """
        从指定高 q 通道提取价电子轮廓，并保存到统一的 pz 网格中。

        Parameters
        ----------
        highq : int
            用于提取价电子轮廓的高动量转移通道编号。

        linrange1 : list or tuple, optional
            第一个用于拟合线性背景和理论轮廓的能量区间 [emin, emax]。

        linrange2 : list or tuple, optional
            第二个用于拟合线性背景和理论轮廓的能量区间 [emin, emax]。

        mirror : bool
            True 时使用 pz <= 0 的价电子轮廓镜像构造另一半；
            False 时使用 Pearson 函数修正近边区域。

        pearson_cut_pz : float, optional
            Pearson 替换边界。若为 None 且 mirror=False，则通过图形交互选择。

        fit_total_J : bool
            True 时拟合理论总 Compton profile self.J；
            False 时只用 core profile self.C 参与扣除。
            原 extractval_test 中使用的是 self.J。

        view : bool
            是否绘图检查。

        Returns
        -------
        valence_pz : ndarray
            提取并修正后的价电子轮廓，定义在 self.pzscale 上。

        asymmetry_pz : ndarray
            价电子非对称项，定义在 self.pzscale 上。
        """

        import numpy as np
        import matplotlib.pyplot as plt
        from scipy import interpolate, optimize
        import pylab

        col = highq

        # 1. 当前高 q 通道的能量轴转换到 pz 轴
        pz = e2pz(
            self.eloss / 1.0e3 + self.E0,
            self.E0,
            self.tth[col]
        )[0]

        # 2. 确定线性背景拟合区间
        if linrange1 is not None and linrange2 is not None:
            range1 = np.where(
                np.logical_and(self.eloss >= linrange1[0], self.eloss <= linrange1[1])
            )[0]
            range2 = np.where(
                np.logical_and(self.eloss >= linrange2[0], self.eloss <= linrange2[1])
            )[0]
            fit_region = np.append(range1, range2)

        elif linrange1 is not None:
            fit_region = np.where(
                np.logical_and(self.eloss >= linrange1[0], self.eloss <= linrange1[1])
            )[0]

        else:
            # 默认选择 core 明显强于 valence 的区域
            fit_region = np.where(0.1 * self.C[:, col] > self.V[:, col])[0]

        if len(fit_region) == 0:
            raise ValueError("线性背景拟合区间为空，请检查 linrange1 / linrange2。")

        # 3. 拟合线性背景 + 理论 Compton profile
        if fit_total_J:
            theory_profile = self.J[:, col]
        else:
            theory_profile = self.C[:, col]

        def fit_function(a):
            linear_background = np.polyval([a[0], a[1]], self.eloss[fit_region])
            scaled_theory = a[2] * theory_profile[fit_region]
            return self.signals[fit_region, col] - linear_background - scaled_theory

        res = optimize.leastsq(fit_function, [0.0, 0.0, 1.0])[0]

        slope = res[0]
        intercept = res[1]
        theory_scale = res[2]

        linear_background = np.polyval([slope, intercept], self.eloss)

        self.background[:, col] = linear_background

        # 4. 扣除线性背景和 core profile，得到原始 valence
        raw_valence = (
            self.signals[:, col]
            - linear_background
            - theory_scale * self.C[:, col]
        )

        # 5. 修正 valence 轮廓
        if mirror:
            # 使用 pz <= 0 的部分进行镜像
            negative_mask = pz <= 0.0

            mirror_val = np.append(
                raw_valence[negative_mask],
                np.flipud(raw_valence[negative_mask])
            )

            mirror_pz = np.append(
                pz[negative_mask],
                np.flipud(-pz[negative_mask])
            )

            order = np.argsort(mirror_pz)

            f_val = interpolate.interp1d(
                mirror_pz[order],
                mirror_val[order],
                bounds_error=False,
                fill_value=0.0
            )

            corrected_valence = f_val(pz)
            asymmetry = np.zeros_like(corrected_valence)

        else:
            # Pearson 函数修正近边区域
            if pearson_cut_pz is None:
                if not view:
                    raise ValueError(
                        "mirror=False 且 view=False 时，必须显式提供 pearson_cut_pz。"
                    )

                plt.figure()
                plt.plot(pz, raw_valence)
                plt.xlabel("pz [a.u.]")
                plt.ylabel("raw valence profile")
                plt.title("请选择 Pearson 替换边界 pz")
                plt.grid(False)
                selected_point = pylab.ginput(1)[0]
                pearson_cut_pz = selected_point[0]
                plt.close()

            edge_region = np.where(pz < pearson_cut_pz)[0]

            if len(edge_region) == 0:
                raise ValueError("Pearson 拟合区间为空，请检查 pearson_cut_pz。")

            max_index_in_edge = edge_region[
                np.argmax(raw_valence[edge_region])
            ]

            start_param = [
                pz[max_index_in_edge],              # peak position
                4.0,                                # FWHM
                1.0,                                # Pearson shape
                raw_valence[max_index_in_edge],     # peak intensity
                0.0                                 # background
            ]

            def pearson_fit_function(a):
                return raw_valence[edge_region] - pearson7(pz[edge_region], a)

            pearson_param = optimize.leastsq(
                pearson_fit_function,
                start_param
            )[0]

            pearson_profile = pearson7(pz, pearson_param)

            corrected_valence = raw_valence.copy()
            corrected_valence[pz > pearson_cut_pz] = pearson_profile[pz > pearson_cut_pz]

            # 6. 拟合 valence 非对称项
            pzp = pz[pz >= 0.0]
            pzm = pz[pz < 0.0]

            val_plus = corrected_valence[pz >= 0.0]
            val_minus_raw = corrected_valence[pz < 0.0]

            if len(pzp) > 0 and len(pzm) > 0:
                f_minus = interpolate.interp1d(
                    -pzm,
                    val_minus_raw,
                    bounds_error=False,
                    fill_value=0.0
                )

                val_minus = f_minus(pzp)

                def asymmetry_fit_function(a):
                    model = a[0] * (
                        np.tanh(pzp / a[1])
                        * np.exp(-(pzp / np.absolute(a[2])) ** 4.0)
                    )
                    return val_plus - val_minus - model

                asym_res = optimize.leastsq(
                    asymmetry_fit_function,
                    [0.0, 1.0, 1.0]
                )[0]

                asymmetry = -(
                    asym_res[0]
                    * np.tanh(pz / asym_res[1])
                    * np.exp(-(pz / np.absolute(asym_res[2])) ** 4.0)
                ) / 2.0

            else:
                asymmetry = np.zeros_like(corrected_valence)

            corrected_valence = corrected_valence - asymmetry

        # 7. 插值到统一的 self.pzscale
        # 注意：原代码中 self.pzscale 是从正到负排列，因此这里保持 flip 逻辑
        order = np.argsort(pz)

        f_corrected_valence = interpolate.interp1d(
            pz[order],
            corrected_valence[order],
            kind="cubic",
            bounds_error=False,
            fill_value=0.0
        )

        f_asymmetry = interpolate.interp1d(
            pz[order],
            asymmetry[order],
            kind="cubic",
            bounds_error=False,
            fill_value=0.0
        )

        valence_pz = f_corrected_valence(self.pzscale)
        asymmetry_pz = f_asymmetry(self.pzscale)

        self.valencepz[:, col] = valence_pz
        self.valasymmetrypz[:, col] = asymmetry_pz

        if view:
            plt.figure()
            plt.plot(pz, raw_valence, label="raw valence")
            plt.plot(pz, corrected_valence, label="corrected valence")
            plt.plot(pz, asymmetry, label="asymmetry")
            plt.xlabel("pz [a.u.]")
            plt.ylabel("intensity")
            plt.legend()
            plt.title(f"High-q valence extraction, q column = {col}")
            plt.grid(False)
            plt.show()

        return valence_pz, asymmetry_pz


    def transform_valence_profile_to_all_q(
        self,
        highq,
        smoothgval=0.0,
        include_asymmetry=True,
        view=False
    ):
        """
        将高 q 提取到的 valence(pz) 转换到所有 q 通道的能量损失轴上。

        Parameters
        ----------
        highq : int
            已经完成 valence 提取的高 q 通道编号。

        smoothgval : float
            高斯平滑 FWHM，单位与 eloss 一致。0 表示不平滑。

        include_asymmetry : bool
            是否加入 valence asymmetry 项。

        view : bool
            是否逐通道绘图检查。

        Returns
        -------
        valence_energy : ndarray
            每个 q 通道对应的 valence(E)，形状与 self.signals 相同。
        """

        import numpy as np
        import matplotlib.pyplot as plt
        from scipy import interpolate

        n_energy = len(self.eloss)
        n_q = len(self.tth)

        new_energy = np.zeros((len(self.pzscale), n_q))
        new_valence = np.zeros((n_energy, n_q))
        new_asymmetry = np.zeros((n_energy, n_q))

        for n in range(n_q):
            # 1. 对当前 q 通道，将统一 pzscale 映射到能量损失轴
            new_energy[:, n] = (
                pz2e1(self.E0, self.pzscale, self.tth[n]) - self.E0
            ) * 1.0e3

            # 2. valence profile 插值到当前 eloss 轴
            f_valence = interpolate.interp1d(
                new_energy[:, n],
                self.valencepz[:, highq],
                bounds_error=False,
                fill_value=0.0
            )

            # 原代码中除以 qvals，保留该强度换算关系
            new_valence[:, n] = f_valence(self.eloss) / self.qvals[:, n]

            # 3. asymmetry profile 同样转换
            f_asymmetry = interpolate.interp1d(
                new_energy[:, n],
                self.valasymmetrypz[:, highq],
                bounds_error=False,
                fill_value=0.0
            )

            new_asymmetry[:, n] = f_asymmetry(self.eloss) / self.qvals[:, n]

            if view:
                plt.figure()
                plt.plot(self.eloss, new_valence[:, n], label="valence")
                plt.plot(self.eloss, new_asymmetry[:, n], label="asymmetry")
                plt.xlabel("energy loss [eV]")
                plt.ylabel("intensity")
                plt.title(f"Transformed valence profile, q column = {n}")
                plt.legend()
                plt.grid(False)
                plt.show()

        # 4. 可选平滑
        if smoothgval > 0.0:
            for n in range(n_q):
                smoothed = convg(self.eloss, new_valence[:, n], smoothgval)

                if include_asymmetry:
                    self.valence[:, n] = smoothed + new_asymmetry[:, n]
                else:
                    self.valence[:, n] = smoothed

                self.valasymmetry[:, n] = new_asymmetry[:, n]

        else:
            if include_asymmetry:
                self.valence = new_valence + new_asymmetry
            else:
                self.valence = new_valence

            self.valasymmetry = new_asymmetry

        return self.valence


    def subtract_valence_profile_from_q(
        self,
        whichq,
        eoffset=0.0,
        fit_start_energy=None,
        fit_core=True,
        allow_energy_shift=True,
        view=True
    ):
        """
        从指定 q 通道中扣除转换后的价电子轮廓。

        拟合模型：
            signal(E) =
                core(E)
                + A * valence(E + dE)
                + linear_background(E)

        最终保存：
            self.sqw[:, col] =
                signal(E)
                - A * valence(E + dE)
                - linear_background(E)

        Parameters
        ----------
        whichq : int or list
            需要扣除 valence 的 q 通道编号。

        eoffset : float
            初始能量偏移，单位 eV。

        fit_start_energy : float, optional
            拟合起始能量。若为 None，则使用 self.prenormrange[0]。

        fit_core : bool
            True 时在拟合残差中扣除 self.C；
            保持与原 remvalenceprof_test 逻辑一致。

        allow_energy_shift : bool
            True 时拟合 valence 能量平移；
            False 时只拟合强度缩放和线性背景。

        view : bool
            是否绘图检查。

        Returns
        -------
        results : dict
            每个 q 通道对应的拟合参数。
        """

        import numpy as np
        import matplotlib.pyplot as plt
        from scipy import interpolate, optimize

        if not isinstance(whichq, list):
            columns = [whichq]
        else:
            columns = whichq

        if fit_start_energy is None:
            fit_start_energy = self.prenormrange[0]

        fit_inds = np.where(self.eloss >= fit_start_energy)[0]

        if len(fit_inds) == 0:
            raise ValueError("扣除拟合区间为空，请检查 fit_start_energy 或 prenormrange。")

        results = {}

        for col in columns:
            # 1. 取出已经转换到该 q 通道的 valence profile
            f_valence = interpolate.interp1d(
                self.eloss + eoffset,
                self.valence[:, col],
                bounds_error=False,
                fill_value=0.0
            )

            base_valence = f_valence(self.eloss)

            # 2. 拟合 valence 缩放、能量平移、线性背景
            if allow_energy_shift:

                def fit_function(a):
                    valence_scale = a[0]
                    energy_shift = a[1]
                    linear_background = np.polyval([a[2], a[3]], self.eloss)

                    shifted_valence = np.interp(
                        self.eloss,
                        self.eloss + energy_shift,
                        base_valence
                    )

                    model_without_core = (
                        valence_scale * shifted_valence
                        + linear_background
                    )

                    if fit_core:
                        residual = (
                            self.signals[fit_inds, col]
                            - self.C[fit_inds, col]
                            - model_without_core[fit_inds]
                        )
                    else:
                        residual = (
                            self.signals[fit_inds, col]
                            - model_without_core[fit_inds]
                        )

                    return residual

                start = [1.0, 0.0, 0.0, 0.0]
                res = optimize.leastsq(fit_function, start)[0]

                valence_scale = res[0]
                energy_shift = res[1]
                slope = res[2]
                intercept = res[3]

            else:

                def fit_function(a):
                    valence_scale = a[0]
                    linear_background = np.polyval([a[1], a[2]], self.eloss)

                    model_without_core = (
                        valence_scale * base_valence
                        + linear_background
                    )

                    if fit_core:
                        residual = (
                            self.signals[fit_inds, col]
                            - self.C[fit_inds, col]
                            - model_without_core[fit_inds]
                        )
                    else:
                        residual = (
                            self.signals[fit_inds, col]
                            - model_without_core[fit_inds]
                        )

                    return residual

                start = [1.0, 0.0, 0.0]
                res = optimize.leastsq(fit_function, start)[0]

                valence_scale = res[0]
                energy_shift = 0.0
                slope = res[1]
                intercept = res[2]

            # 3. 构造最终用于扣除的 valence 和 background
            shifted_valence = np.interp(
                self.eloss,
                self.eloss + energy_shift,
                base_valence
            )

            fitted_valence = valence_scale * shifted_valence
            fitted_background = np.polyval([slope, intercept], self.eloss)

            # 4. 扣除 valence 和线性背景
            self.sqw[:, col] = (
                self.signals[:, col]
                - fitted_valence
                - fitted_background
            )

            self.valence[:, col] = fitted_valence
            self.background[:, col] = fitted_background

            results[col] = {
                "valence_scale": valence_scale,
                "energy_shift": energy_shift,
                "background_slope": slope,
                "background_intercept": intercept
            }

            if view:
                plt.figure()
                plt.plot(self.eloss, self.signals[:, col], label="raw signal")
                plt.plot(
                    self.eloss,
                    fitted_valence + fitted_background,
                    label="valence + linear background"
                )
                plt.plot(self.eloss, fitted_valence, label="fitted valence")
                plt.plot(self.eloss, fitted_background, label="linear background")
                plt.plot(self.eloss, self.sqw[:, col], label="subtracted S(q,w)")
                plt.plot(self.eloss, self.C[:, col], label="HF core profile")
                plt.xlabel("energy loss [eV]")
                plt.ylabel("intensity")
                plt.title(f"Valence subtraction, q column = {col}")
                plt.legend()
                plt.grid(False)
                plt.show()

        return results


    def run_valence_subtraction_workflow(
        self,
        highq,
        target_qs,
        linrange1=None,
        linrange2=None,
        mirror=False,
        pearson_cut_pz=None,
        smoothgval=0.0,
        eoffset=0.0,
        fit_start_energy=None,
        view=True
    ):
        """
        完整价电子轮廓提取与扣除流程。

        流程：
            1. 从 highq 通道提取 valence(pz)
            2. 将 valence(pz) 转换到所有 q 通道
            3. 对 target_qs 中的每个 q 通道拟合并扣除 valence
            4. 返回每个 q 通道的拟合参数

        Parameters
        ----------
        highq : int
            用于提取价电子轮廓的高 q 通道。

        target_qs : list
            需要扣除价电子贡献的 q 通道列表。

        linrange1, linrange2 : list or tuple
            高 q 提取 valence 时的拟合区间。

        mirror : bool
            是否用镜像法修正 valence。

        pearson_cut_pz : float
            Pearson 替换边界。

        smoothgval : float
            转换到所有 q 后的 valence 平滑宽度。

        eoffset : float
            扣除时的初始能量偏移。

        fit_start_energy : float
            扣除时的拟合起始能量。

        view : bool
            是否绘图检查。

        Returns
        -------
        subtraction_results : dict
            每个目标 q 通道的扣除拟合参数。
        """

        self.extract_highq_valence_profile(
            highq=highq,
            linrange1=linrange1,
            linrange2=linrange2,
            mirror=mirror,
            pearson_cut_pz=pearson_cut_pz,
            fit_total_J=True,
            view=view
        )

        self.transform_valence_profile_to_all_q(
            highq=highq,
            smoothgval=smoothgval,
            include_asymmetry=True,
            view=False
        )

        subtraction_results = self.subtract_valence_profile_from_q(
            whichq=target_qs,
            eoffset=eoffset,
            fit_start_energy=fit_start_energy,
            fit_core=True,
            allow_energy_shift=True,
            view=view
        )

        return subtraction_results






