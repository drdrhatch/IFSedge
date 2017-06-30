#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np
import optparse as op
import matplotlib.pyplot as plt
from finite_differences import *
from interp import *
from read_EFIT_file import *

parser=op.OptionParser(description='Calculates Er/(R B_pol) from a gene profiles file and an efit file.  Arguments: gene_profiles_file_name_i, efit_file_name.')
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include gene_profiles_file_name_i and efit_file_name."
    \n""")

pfile = args[0]
efit_file_name = args[1]

data = np.genfromtxt(pfile) 

dummy = raw_input("WARNING: profile file must have rho_poloidal as second column! Press any key to continue.\n")
dummy = raw_input("Assuming ion charge is 1 and reference mass is deuterium (press any key to continue).\n")

rhot = data[:,0]
psi = data[:,1]**2
Ti = data[:,2]
ni = data[:,3]

mi = 1.673e-27
ee = 1.602e-19
mref = 2.0*mi
Z = 1.0

Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,0.9)

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)

print "R_major",R_major
print "Bref",Bref
print "psisep",psisep
print "psiax",psiax
psisep0 = psisep-psiax
print "psisep0",psisep0


#Get everything on the same even psi grid
psi0 = np.linspace(0.0,1.0,num=3000)
#interopolate rho_tor, n, T, qpsi
ni0 = interp(psi,ni,psi0)
Ti0 = interp(psi,Ti,psi0)
Ti0J = Ti0*1000.0*ee
ni00 = ni0*1.0e19
pi0 = Ti0J*ni00
rhot0 = interp(psi,rhot,psi0)
qpsi0 = interp(psip_n,qpsi,psi0)
R0 = interp(psip_n,Rgrid,psi0)


trpeps = rhot0*Lref/R_major
vti = (0.5*Ti0*1000.0*ee/mref)**0.5
#nuii = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-log(sqrt(ne*1.0E13)/Te*0.001))
coll = 2.3031E-5*Lref*(ni0)/(Ti0)**2*(24.0-log(sqrt(ni0*1.0E13)/Ti0*0.001))
nustar_i=8.0/3.0/pi**0.5*qpsi0/trpeps**1.5*(R_major/Lref)*Z**4*coll
#plt.plot(psi0,nustar_i)
#plt.title('nustar_i')
#plt.show()

ft = trpeps**0.5*(1.0+0.46*(1.0-trpeps))
fc = 1.0 - ft

a = 1.0/(1.0+0.5*nustar_i**0.5)
b = -1.17*fc/(1.0-0.22*ft-0.19*ft**2) + 0.25*(1-ft**2)*nustar_i**0.5
c =    0.315*nustar_i**2*ft**6
d = 1.0/(1.0+0.15*nustar_i**2*ft**6)
kpar = -(a*b+c)*d
#plt.plot(psi0,kpar)
#plt.title('kpar')
#plt.show()

dTdpsi = fd_d1_o4(Ti0J,psi0)
dndpsi = fd_d1_o4(ni00,psi0)
dpdpsi = fd_d1_o4(pi0,psi0)
dRdpsi = fd_d1_o4(R0,psi0)

omegator = (1.0/psisep0/ee)*((1-kpar)*dTdpsi + Ti0J/ni00*dndpsi)
omegator0 = (1.0/psisep0/ee)*1/ni00*dpdpsi

f=open('omegator'+pfile,'w')
f.write('#1.rho_tor 2.psiN 3.Er/(R Bpol) 4.Er0(gradP)/(R Bpol) \n')
np.savetxt(f,np.column_stack((rhot0,psi0,omegator,omegator0)))
f.close()



