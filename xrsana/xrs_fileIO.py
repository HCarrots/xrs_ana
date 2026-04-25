#!/usr/bin/python
# Filename: xrs_utilities.py

__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import numpy as np


def readbiggsdata(filename,element):
    """
    Reads Hartree-Fock Profile of element 'element' from values tabulated 
    by Biggs et al. (Atomic Data and Nuclear Data Tables 16, 201-309 (1975))
    as provided by the DABAX library (http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/ComptonProfiles.dat).
    input:
    filename = path to the ComptonProfiles.dat file (the file should be distributed with this package)
    element  = string of element name
    returns:
    data     = the data for the according element as in the file:
        #UD  Columns: 
        #UD  col1: pz in atomic units 
        #UD  col2: Total compton profile (sum over the atomic electrons
        #UD  col3,...coln: Compton profile for the individual sub-shells
    occupation = occupation number of the according shells
    bindingen  = binding energies of the accorting shells
    colnames   = strings of column names as used in the file
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








