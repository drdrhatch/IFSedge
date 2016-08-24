#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from interp import *

def first_derivative(f_in,x_in):
    x = x_in.flatten()
    f = f_in.flatten()
    dx = x[1]-x[0]
    dfdx = np.empty(len(f))
    for i in range(len(f)):
	if i == 0:
	    dfdx[i] = (f[i+1] - f[i])/dx
	elif i == 1 or i == len(f)-2:
	    dfdx[i] = (f[i+1]-f[i-1])/2./dx
	elif i == len(f)-1:
	    dfdx[i] = (f[i] - f[i-1])/dx
	else:
	    dfdx[i] = (-f[i+2]+8.*f[i+1]-8.*f[i-1]+f[i-2])/12./dx
	
    return dfdx

def rhot2psip(rbsProfs_file_name,flux_surface):
    rhot_fs = float(flux_surface)
    rbsProfs_data = np.genfromtxt(rbsProfs_file_name)
    psip_data_full = rbsProfs_data[:,1]
    rhot_data_full = rbsProfs_data[:,0]
    ind = np.argmin(abs(rhot_data_full-rhot_fs))
    psip_data = rbsProfs_data[ind-10:ind+10,1]
    rhot_data = rbsProfs_data[ind-10:ind+10,0]
    rhot_psip_spl = interpolate.UnivariateSpline(rhot_data,psip_data)
    psip_fs = rhot_psip_spl(rhot_fs)
    print 'rhot = ',rhot_fs
    print 'psip = ',psip_fs
    return psip_fs

def read_rbsProfs(rbsProfs_file_name,flux_surface):

    n_spl = 10

    rbsProfs_data = np.genfromtxt(rbsProfs_file_name)

    rhot_fs = float(flux_surface)
    rhot_data_full = rbsProfs_data[:,0]
    ind = np.argmin(abs(rhot_data_full-rhot_fs))

    psip_data = rbsProfs_data[ind-n_spl:ind+n_spl,1]
    rhot_data = rbsProfs_data[ind-n_spl:ind+n_spl,0]

    #rhot_data_unif = np.linspace(rhot_data[0],rhot_data[-1],len(rhot_data)*1000)
    #psip_data_unif = interp(rhot_data,psip_data,rhot_data_unif)

    R_obmp = rbsProfs_data[ind,24]
    print 'R_obmp = ', R_obmp

    n_data = rbsProfs_data[ind-n_spl:ind+n_spl,3]
    n_rhot_spl = interpolate.UnivariateSpline(rhot_data,n_data)
    fprime_obmp = abs(n_rhot_spl(rhot_fs,nu=1)/n_rhot_spl(rhot_fs))
    print 'n = ', n_rhot_spl(rhot_fs)
    print 'fprime=',fprime_obmp 
    
    Ti_data = rbsProfs_data[ind-n_spl:ind+n_spl,4]
    Ti_rhot_spl = interpolate.UnivariateSpline(rhot_data,Ti_data)
    tprime_obmp = abs(Ti_rhot_spl(rhot_fs,nu=1)/Ti_rhot_spl(rhot_fs))
    print 'Ti = ', Ti_rhot_spl(rhot_fs)
    print 'tprime=',tprime_obmp

    q_data = rbsProfs_data[ind-n_spl:ind+n_spl,23]
    q_rhot_spl = interpolate.UnivariateSpline(rhot_data,q_data)
    shat_obmp = q_rhot_spl(rhot_fs,nu=1)/q_rhot_spl(rhot_fs)*rhot_fs
    print 'q = ',q_rhot_spl(rhot_fs)
    print 'shat=',shat_obmp

    Bp_obmp = rbsProfs_data[ind,25]
    Bt_obmp = rbsProfs_data[ind,26]
    gamE_obmp = rbsProfs_data[ind,9]
    print 'gamE_obmp = ', gamE_obmp

    B_tot_obmp = np.sqrt(Bp_obmp**2+Bt_obmp**2)

    return R_obmp,Bp_obmp,abs(Bt_obmp),B_tot_obmp,gamE_obmp,abs(tprime_obmp),abs(fprime_obmp),shat_obmp
