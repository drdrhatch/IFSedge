#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *
from read_EFIT_file_new import *
#from fields_along_fluxsurface import *
from read_rbsProfs import *
from calc_fields_from_EFIT import *

parser = op.OptionParser()
options,args = parser.parse_args()
efit_file_name = args[0]
rbsProfs_file_name = args[1]
rhot_fs = args[2]
#psip_fs = rhot_fs
psip_fs = rhot2psip(rbsProfs_file_name,rhot_fs)

R_fs, B_pol, B_tor, B_tot = flux_surface_B_fields(efit_file_name,psip_fs)

plt.plot(R_fs,B_pol,label='B_pol_x')
plt.plot(R_fs,abs(B_tor),label='B_tor_x')
##plt.title(efit_file_name+' at flux surface of '+psip_fs[0]+' +- 0.01')
plt.xlabel('R')
plt.ylabel('B')
plt.legend()
plt.show()

R_obmp,Bp_obmp,Bt_obmp,B_tot_obmp,gamE_obmp,tprime_obmp,fprime_obmp,shat_obmp = read_rbsProfs(rbsProfs_file_name,rhot_fs)
print R_obmp,Bp_obmp,Bt_obmp,B_tot_obmp,gamE_obmp,tprime_obmp,fprime_obmp,shat_obmp
tprime_fs = R_fs*B_pol*tprime_obmp/R_obmp/Bp_obmp
fprime_fs = R_fs*B_pol*fprime_obmp/R_obmp/Bp_obmp
gammaE_fs = (R_fs*B_pol)**2/B_tot*gamE_obmp*B_tot_obmp/(R_obmp*Bp_obmp)**2

print 'min abs(gamE)=',min(abs(gammaE_fs))
plt.plot(R_fs,tprime_fs,'.',label='tprime')
plt.plot(R_fs,fprime_fs,'.',label='fprime')
f = open('gamE.dat','w')
f.write('#1.R_fs 2.gamE \n#\n')
np.savetxt(f,np.column_stack((R_fs,abs(gammaE_fs))))
f.close()
plt.xlabel('R')
plt.legend(loc=2)
##plt.title(rbsProfs_file_name+' fs= '+str(fs))
plt.show()
plt.plot(R_fs,abs(gammaE_fs),'.',label='abs(gammaE)')
plt.xlabel('R')
plt.legend(loc=4)
#plt.title(rbsProfs_file_name+' fs= '+str(fs))
plt.show()

