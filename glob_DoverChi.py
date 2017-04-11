#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from finite_differences import *
from interp import *
import optparse as op

def readProfiles(speciesName, suffix):

    profilesFilename = 'profiles_'+speciesName+'_'+ suffix
    fluxprofFilename = 'fluxprof'+speciesName+'_'+ suffix+'.dat'

    profiles = np.genfromtxt(profilesFilename)
    rho = profiles[:,0]
    tempProf = profiles[:,2]
    densProf = profiles[:,3]

    fluxprof = np.genfromtxt(fluxprofFilename,skip_footer=3)
    fluxrho = fluxprof[:,0]
    Gamma_es = fluxprof[:,1]
    Qheat_es = fluxprof[:,2]
    Gamma_em = fluxprof[:,3]
    Qheat_em = fluxprof[:,4]
    Gamma_tot = Gamma_es + Gamma_em
    Qheat_tot = Qheat_es + Qheat_em
    Qheat_tot = Qheat_tot - 3./2.*Gamma_tot

    return rho, tempProf, densProf, fluxrho, Gamma_tot, Qheat_tot

def gradProf(rho, t, n, fluxrho):

    dtdrho = fd_d1_o4(t,rho)
    dndrho = fd_d1_o4(n,rho)

    dndrho_fluxrho = interp(rho, dndrho, fluxrho)
    ndtdrho_fluxrho = interp(rho, n*dtdrho, fluxrho)

    return dndrho_fluxrho, ndtdrho_fluxrho
    
def weightedAvg(fluxrho, Gamma, Qheat, dndrho, ndtdrho):

    weight = np.abs(Qheat)
    norm = np.sum(weight)

    avgGamma = np.sum(Gamma * weight) / norm
    avgQheat = np.sum(Qheat * weight) / norm

    D = avgGamma / (np.sum(dndrho * weight) / norm)
    chi = avgQheat / (np.sum(ndtdrho * weight) / norm)

    return avgGamma, avgQheat, D, chi

parser = op.OptionParser()
options, args = parser.parse_args()
suffix = args[0]

species = ['i','e','z']
Gamma = np.empty(len(species),dtype='float')
Qheat = np.empty(len(species),dtype='float')
D = np.empty(len(species),dtype='float')
chi = np.empty(len(species),dtype='float')

for i in range(len(species)):
    s = species[i]
    rho, t, n, fluxrho, Gamma, Qheat = readProfiles(s, suffix)
    dndrho, ndtdrho = gradProf(rho, t, n, fluxrho)
    Gamma[i], Qheat[i], D[i], chi[i] = weightedAvg(fluxrho, Gamma, Qheat, dndrho, ndtdrho)
    outstr = 'D_'+s+'/ chi_'+s+' ='
    print outstr, D[i] / chi[i]

chi_tot = np.sum(chi)
print 'D_i / chi_tot =', D[0] / chi_tot
print 'D_e / chi_tot =', D[1] / chi_tot
print 'D_z / chi_tot =', D[2] / chi_tot


