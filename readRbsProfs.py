#! /usr/bin/python

from read_EFIT import *
from read_iterdb_file import *
#from calc_fields_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys
from write_iterdb import *

def readRbsProfs(rbsProfs_fileName):

    rbsProfs = np.genfromtxt(rbsProfs_fileName, skip_header = 3)
    RBSdict = {}
    RBSdict['rhotn'] = rbsProfs[:, 0]
    RBSdict['psipn'] = rbsProfs[:, 1]
    RBSdict['R'] = rbsProfs[:, 24]

    RBSdict['Pres'] = rbsProfs[:, 2]    # total pressure N / m^2
    RBSdict['ne'] = rbsProfs[:, 3]    # unit: 10^20 m^-3
    RBSdict['ti'] = rbsProfs[:, 4] * 1.E03    # unit: eV
    RBSdict['te'] = rbsProfs[:, 5] * 1.E03    # unit: eV

    RBSdict['Er'] = rbsProfs[:, 16]    # unit: V / m
    RBSdict['qpsi'] = rbsProfs[:, 23]    # safety factor q
    RBSdict['Bpol'] = rbsProfs[:, 25]    # B poloidal
    RBSdict['Btor'] = abs(rbsProfs[:, 26])    # B toroidal
    
    return RBSdict

