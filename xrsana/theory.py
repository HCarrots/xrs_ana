#!/usr/bin/python
# Filename: theory.py

from .xrs_utilities import *
from .math_functions import *
import string
import numpy as np
from numpy import array
from scipy import interpolate, signal, integrate, constants, optimize


import matplotlib.pyplot as plt

__metaclass__ = type # new style classes


class HFspecpredict:
    def __init__(self,formulas,concentrations,rho_formu,correctasym=None,E0=9.68,eloss=np.arange(0,1,0.0001),alpha=None,beta=None,samthick=None):
        if not isinstance(formulas,list):
            theformulas = []
            theformulas.append(formulas)
        else:
            theformulas = formulas
        self.formulas       = theformulas
        self.concentrations = concentrations
        self.E0             = E0
        self.eloss          = eloss
        self.rho_formu      = rho_formu
        self.rho            = 0.0
        if len(rho_formu)>1:
            for n in range(len(rho_formu)):
                self.rho += rho_formu[n]*concentrations[n]
        else:
            self.rho = rho_formu
        if not correctasym:
            correctasym = []
            for formula in formulas:
                elements,stoichiometries = parseformula(formula)
                correctasym.append(np.zeros(len(elements)))
        self.correctasym    = correctasym
        self.alpha          = alpha
        self.beta           = beta
        if self.beta<0: # transmission geometry
            self.tth = alpha-beta
        else: # reflection geometry
            self.tth = 180.0 - (alpha+beta)
        self.thickness      = samthick # in [cm] now

        self.eloss,self.J,self.C,self.V,self.q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0,tth=self.tth,correctasym=self.correctasym)
        # sample self absorption
        self.alpha_r = alpha
        self.beta_r  = beta
        self.mu_in,self.mu_out = mpr_compds(self.eloss/1e3+self.E0,self.formulas,self.concentrations,self.E0,self.rho_formu)
        self.ac = abscorr2(self.mu_in,self.mu_out,self.alpha_r,self.beta_r,self.thickness)
        self.J  = self.J/self.ac*self.rho
        self.C  = self.C/self.ac*self.rho
        self.V  = self.V/self.ac*self.rho

    def plotHFspec(self):
        plt.plot(self.eloss,self.J,self.eloss,self.C,self.eloss,self.V)
        plt.legend(('sum','core contribution','valence contribution'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('S(q,w) [1/eV]')
        #plt.title('About as simple as it gets, folks')
        plt.grid(False)
        plt.show(block=False)

    def plotmurho(self):
        pass
    def plotresult(self):
        pass

class HFspecpredict_series:
    def __init__(self,formulas,concentrations,rho_formu,correctasym=None,E0=9.68,eloss=np.arange(0,1,0.0001),alpha=0.0,beta=-30.0,samthick=0.1):
        self.concentrations = concentrations
        self.eloss          = eloss
        self.rho_formu      = rho_formu
        self.rho            = 0
        if len(rho_formu)>1:
            for n in range(len(rho_formu)):
                self.rho += rho_formu[n]*concentrations[n]
        else:
            self.rho = rho_formu
        if not correctasym:
            correctasym = []
            for formula in formulas:
                elements,stoichiometries = parseformula(formula)
                correctasym.append(np.zeros(len(elements)))
        self.correctasym    = correctasym
        # make everything iterable
        # formulas
        if not isinstance(formulas,list):
            theformulas = []
            theformulas.append(formulas)
        else:
            theformulas = formulas
        self.formulas       = theformulas
        # E0
        if not isinstance(E0,list):
            theE0s = []
            theE0s.append(E0)
        else:
            theE0s = E0        
        self.E0             = theE0s

        # alpha
        if not isinstance(alpha,list):
            thealphas = []
            thealphas.append(alpha)
        else:
            thealphas = alpha
        self.alpha          = thealphas
        # beta
        if not isinstance(beta,list):
            thebetas = []
            thebetas.append(beta)
        else:
            thebetas = beta
        self.beta        = thebetas
        # tth
        self.tth = []
        for anglea in self.alpha:
            for angleb in self.beta:
                if angleb<0: # transmission geometry
                    tth = anglea-angleb
                    self.tth.append(tth)
                else: # reflection geometry
                    tth = 180.0 - (anglea+angleb)
                    self.tth.append(tth)
        # sample thickness
        if not isinstance(samthick,list):
            thethickness = []
            thethickness.append(samthick)
        else:
            thethickness = samthick
        self.thickness      = thethickness # in [cm] now

        # now calculate spectra for all possible configurations (several E0, several tth, several samthick)
        # one E0, one tth, one thickness
        if len(self.E0)==1 and len(self.alpha)==1 and len(self.beta)==1 and len(self.thickness)==1:
            self.eloss,self.J,self.C,self.V,self.q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
            # sample self absorption
            self.mu_in,self.mu_out = mpr_compds(self.eloss/1e3+self.E0,self.formulas,self.concentrations,self.E0,self.rho_formu)
            self.ac = abscorr2(self.mu_in,self.mu_out,self.alpha[0],self.beta[0],self.thickness[0])
            self.J  = self.J/self.ac*self.rho
            self.C  = self.C/self.ac*self.rho
            self.V  = self.V/self.ac*self.rho

        # several E0, one tth, one thickness
        if len(self.E0) > 1:
            eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
            self.eloss = eloss
            self.J     = np.zeros((len(eloss),len(self.E0)))
            self.C     = np.zeros((len(eloss),len(self.E0)))
            self.V     = np.zeros((len(eloss),len(self.E0)))
            self.q     = np.zeros((len(eloss),len(self.E0)))
            self.ac    = np.zeros((len(eloss),len(self.E0)))
            for n in range(len(self.E0)):
                eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[n],tth=self.tth[0],correctasym=self.correctasym)
                mu_in, mu_out = mpr_compds(eloss/1e3+self.E0[n],self.formulas,self.concentrations,self.E0[n],self.rho_formu)
                ac = abscorr2(mu_in,mu_out,self.alpha[0],self.beta[0],self.thickness[0])
                j  = j/ac*self.rho
                c  = c/ac*self.rho
                v  = v/ac*self.rho
                self.J[:,n]  = np.interp(self.eloss,eloss,j)
                self.C[:,n]  = np.interp(self.eloss,eloss,c)
                self.V[:,n]  = np.interp(self.eloss,eloss,v)
                self.q[:,n]  = np.interp(self.eloss,eloss,q)
                self.ac[:,n] = np.interp(self.eloss,eloss,ac)

        # several incidence angles:
        if len(self.alpha) > 1:
            eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
            self.eloss = eloss
            self.J     = np.zeros((len(eloss),len(self.alpha)))
            self.C     = np.zeros((len(eloss),len(self.alpha)))
            self.V     = np.zeros((len(eloss),len(self.alpha)))
            self.q     = np.zeros((len(eloss),len(self.alpha)))
            self.ac    = np.zeros((len(eloss),len(self.alpha)))
            for n in range(len(self.alpha)):
                eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[n],correctasym=self.correctasym)
                mu_in, mu_out = mpr_compds(eloss/1e3+self.E0[0],self.formulas,self.concentrations,self.E0[0],self.rho_formu)
                ac = abscorr2(mu_in,mu_out,self.alpha[n],self.beta[0],self.thickness[0])
                j  = j/ac*self.rho
                c  = c/ac*self.rho
                v  = v/ac*self.rho
                self.J[:,n]  = np.interp(self.eloss,eloss,j)
                self.C[:,n]  = np.interp(self.eloss,eloss,c)
                self.V[:,n]  = np.interp(self.eloss,eloss,v)
                self.q[:,n]  = np.interp(self.eloss,eloss,q)
                self.ac[:,n] = np.interp(self.eloss,eloss,ac)
    
        # several exit angles:
        if len(self.beta) > 1:
            eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
            self.eloss = eloss
            self.J     = np.zeros((len(eloss),len(self.beta)))
            self.C     = np.zeros((len(eloss),len(self.beta)))
            self.V     = np.zeros((len(eloss),len(self.beta)))
            self.q     = np.zeros((len(eloss),len(self.beta)))
            self.ac    = np.zeros((len(eloss),len(self.beta)))
            for n in range(len(self.beta)):
                eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[n],correctasym=self.correctasym)
                mu_in, mu_out = mpr_compds(eloss/1e3+self.E0[0],self.formulas,self.concentrations,self.E0[0],self.rho_formu)
                ac = abscorr2(mu_in,mu_out,self.alpha[0],self.beta[n],self.thickness[0])
                j  = j/ac*self.rho
                c  = c/ac*self.rho
                v  = v/ac*self.rho
                self.J[:,n]  = np.interp(self.eloss,eloss,j)
                self.C[:,n]  = np.interp(self.eloss,eloss,c)
                self.V[:,n]  = np.interp(self.eloss,eloss,v)
                self.q[:,n]  = np.interp(self.eloss,eloss,q)
                self.ac[:,n] = np.interp(self.eloss,eloss,ac)

        # several sample thicknesses:
        if len(self.thickness) > 1:
            eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
            self.eloss = eloss
            self.J     = np.zeros((len(eloss),len(self.thickness)))
            self.C     = np.zeros((len(eloss),len(self.thickness)))
            self.V     = np.zeros((len(eloss),len(self.thickness)))
            self.q     = np.zeros((len(eloss),len(self.thickness)))
            self.ac    = np.zeros((len(eloss),len(self.thickness)))
            for n in range(len(self.thickness)):
                eloss,j,c,v,q = makeprofile_compds(self.formulas,self.concentrations,E0=self.E0[0],tth=self.tth[0],correctasym=self.correctasym)
                mu_in, mu_out = mpr_compds(eloss/1e3+self.E0[0],self.formulas,self.concentrations,self.E0[0],self.rho_formu)
                ac = abscorr2(mu_in,mu_out,self.alpha[0],self.beta[0],self.thickness[n])
                j  = j/ac*self.rho
                c  = c/ac*self.rho
                v  = v/ac*self.rho
                self.J[:,n]  = np.interp(self.eloss,eloss,j)
                self.C[:,n]  = np.interp(self.eloss,eloss,c)
                self.V[:,n]  = np.interp(self.eloss,eloss,v)
                self.q[:,n]  = np.interp(self.eloss,eloss,q)
                self.ac[:,n] = np.interp(self.eloss,eloss,ac)

    def plotHFspec(self):
        plt.plot(self.eloss,self.J,self.eloss,self.C,self.eloss,self.V)
        plt.legend(('sum','core contribution','valence contribution'))
        plt.xlabel('energy loss [eV]')
        plt.ylabel('S(q,w) [1/eV]')
        #plt.title('About as simple as it gets, folks')
        plt.grid(False)
        plt.show(block=False)

class HFspectrum:



    def __init__(self,data,formulas,concentrations,correctasym=None,correctasym_pertth=None, initialise=True):
        """
        class for building S(q,w) from tabulated Hartree-Fock Compton profiles to use in
        the extraction algorithm.
        data = instance of one of the classes from the xrs_read module (like read_id20)
        formulas = single string or list of strings of chemical sum formulas of which the sample is made up
        concentrations = single value or list of concentrations of how the different chemical formulas are mixed
                 (sum should be 1)
        correctasym = single value or list of scaling values for the HR-correction to 
                  the 1s, 2s, and 2p shells. one value per element in the list of formulas
        correctasym_pertth = list of additional scaling values for each scattering angle, so that 
                a momentum transfer dependent asymmetry correction is possible (i.e. no 
                correction for low q, finite correction for high q data) 
        """
        if not initialise:
            return
        # from raw data
        self.eloss = data.eloss
        self.tth   = data.tth
        self.E0    = data.E0
        self.cenom = data.cenom
        # new info
        self.formulas       = formulas
        self.concentrations = concentrations
        if not correctasym:
            correctasym = []
            for formula in formulas:
                elements,stoichiometries = parseformula(formula)
                correctasym.append(np.zeros(len(elements)))
        self.correctasym = correctasym
        # make one profile per scattering angle and spline it onto the exp. energy loss scale
        self.J = np.zeros(np.shape(data.signals))
        self.C = np.zeros(np.shape(data.signals))
        self.V = np.zeros(np.shape(data.signals))
        self.q = np.zeros(np.shape(data.signals))
        if not correctasym_pertth:
            for n in [ nn for nn in range(len(data.signals[0,:])) if  ((data.signals[:,nn]).sum()>0) ]:
                

                el,j,c,v,q = makeprofile_compds(formulas,concentrations,E0=self.cenom[n],tth=data.tth[n],correctasym=self.correctasym)
                f = interpolate.interp1d(el,j, bounds_error=False, fill_value=0.0)
                self.J[:,n] = f(data.eloss)
                f = interpolate.interp1d(el,c, bounds_error=False, fill_value=0.0)
                self.C[:,n] = f(data.eloss)
                f = interpolate.interp1d(el,v, bounds_error=False, fill_value=0.0)
                self.V[:,n] = f(data.eloss)
                f = interpolate.interp1d(el,q, bounds_error=False, fill_value=0.0)
                self.q[:,n] = f(data.eloss)
        else:
            if len(correctasym_pertth) != len(self.tth):
                print( 'Please provide a Python list of one scaling factor [0,1] per scattering angle!')
                print( 'Currently %d scattering angles defined, but %d scaling factors provided!' % (len(self.tth), len(correctasym_pertth)))
                return
            else:
                for n in [ nn for nn in  range(len(data.signals[0,:])) if  ((data.signals[:,nn]).sum()>0)     ]:
                    print( self.correctasym, correctasym_pertth[n])

                    el,j,c,v,q = makeprofile_compds(formulas,concentrations,E0=self.cenom[n],tth=data.tth[n],correctasym=np.array(self.correctasym)*correctasym_pertth[n])
                    f = interpolate.interp1d(el,j, bounds_error=False, fill_value=0.0)
                    self.J[:,n] = f(data.eloss)
                    f = interpolate.interp1d(el,c, bounds_error=False, fill_value=0.0)
                    self.C[:,n] = f(data.eloss)
                    f = interpolate.interp1d(el,v, bounds_error=False, fill_value=0.0)
                    self.V[:,n] = f(data.eloss)
                    f = interpolate.interp1d(el,q, bounds_error=False, fill_value=0.0)
                    self.q[:,n] = f(data.eloss)

        # correct interpolation errors in q (first couple of values are 0.0, replace by smallest value) THIS NEEDS TO BE FIXED
        for n in [ nn for nn in range(len(data.signals[0,:])) if  ((data.signals[:,nn]).sum()>0) ]:
            inds = np.where(self.q[:,n] == 0)
            self.q[inds,n] = self.q[np.where(self.q[:,n]>0)[0][0],n]

        def save_state_hdf5(self, filename, groupname, comment=""):
            import h5py
            h5 = h5py.File(filename,"a")

            h5.require_group(groupname)
            h5group =  h5[groupname]
            if(   "J" in list(h5group.keys()) ):
                raise Exception(" Read data already present in  " + filename+":"+groupname)

            for key in ["eloss" ,"tth" ,"E0" ,"cenom","formulas" ,"concentrations" ,"correctasym" ,"J" ,"C" ,"V", "q"]: 
                data = getattr(self,key)
                if key == "concentrations":
                    s=data[0]
                    for t in data[1:]:
                        s=s+" "+t
                    data=s
                h5group[key]  =  data

            h5group["comment"]  = comment
            h5.flush()
            h5.close()

        def load_state_hdf5(self, filename, groupname):
            import h5py
            h5 = h5py.File(filename,"r")

            h5group =  h5[groupname]
            chiavi = {
            "eloss":array ,
            "tth":array ,
            "E0":float ,
            "cenom":array,
            "formulas":array ,
            "concentrations":array ,
            "correctasym":array ,
            "J":array ,
            "C":array ,
            "V":array, 
            "q":array
            }
            
            for key in  chiavi: 
                if key== "concentrations":
                    data=str(h5group[key])
                    data=string.split(data," ")
                    setattr(self,key, data)
                else:
                    setattr(self,key, chiavi[key](array(h5group[key]))  )
            h5.flush()
            h5.close()


    def plotHFC(self):
        plt.clf()
        for n in range(len(self.C[0,:])):
            plt.plot(self.eloss,self.C[:,n])
        plt.title('core HF profile')
        plt.xlabel('energy loss [eV]')
        plt.ylabel('S(q,w) [1/eV]')
        plt.grid(False)
        plt.show()

    def plotHFJ(self):
        plt.clf()
        for n in range(len(self.J[0,:])):
            plt.plot(self.eloss,self.J[:,n])
        plt.title('total HF profile')
        plt.xlabel('energy loss [eV]')
        plt.ylabel('S(q,w) [1/eV]')
        plt.grid(False)
        plt.show()

    def plotHFV(self):
        plt.clf()
        for n in range(len(self.V[0,:])):
            plt.plot(self.eloss,self.V[:,n])
        plt.title('valence HF profile')
        plt.xlabel('energy loss [eV]')
        plt.ylabel('S(q,w) [1/eV]')
        plt.grid(False)
        plt.show()








