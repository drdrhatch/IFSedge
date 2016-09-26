#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *


parser=op.OptionParser(description='Calculates gradient scale lengths and eta from gene profile files (outputs full info to data files and outputs to screen info at selected rhotor).  Arguments: profiles_e, profiles_i, rho_tor')
options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include profiles_e and profiles_i file names."
    \n""")

f_ge = args[0]
f_gi = args[1]
rhot0 = float(args[2])

gene_e = np.genfromtxt(f_ge)
gene_i = np.genfromtxt(f_gi)

rhote = gene_e[:,0]
te = gene_e[:,2]
ne = gene_e[:,3]

rhoti = gene_i[:,0]
ti = gene_i[:,2]
ni = gene_i[:,3]

rhot1 = np.arange(1000)/999.0*(rhote[-1]-rhote[0])+rhote[0]
n1 = interp(rhote,ne,rhot1)
t1 = interp(rhote,te,rhot1)
rhot2 = np.arange(1000)/999.0*(rhoti[-1]-rhoti[0])+rhoti[0]
n2 = interp(rhoti,ni,rhot2)
t2 = interp(rhoti,ti,rhot2)

omt1 = abs(1.0/t1*fd_d1_o4(t1,rhot1))
omn1 = abs(1.0/n1*fd_d1_o4(n1,rhot1))
omt2 = abs(1.0/t2*fd_d1_o4(t2,rhot2))
omn2 = abs(1.0/n2*fd_d1_o4(n2,rhot2))


#plt.figure(figsize=[8.0,3.5])
#plt.subplot(121)
#plt.plot(rhot1,n1,linewidth=2,color='black',label='n e')
#plt.plot(rhot2,n2,linewidth=2,color='red',label='n i')
#plt.xlabel(r'$\rho_{tor}$',labelpad=10,size=14)
#plt.ylabel(r'$n(10^{19}m^{-3})$',size=14)
#ax=plt.axis()
##plt.axis((0.96,1.0,0.0,0.7))
#plt.xticks((0.96,0.97,0.98,0.99,1.0))
#plt.legend(loc='lower left')
#plt.subplot(122)
#plt.plot(rhot1,omt1/omn1,linewidth=2,color='black',label='eta e')
#plt.plot(rhot1,omt2/omn2,linewidth=2,color='red',label='eta i')
#ax=plt.axis()
##plt.axis((0.96,1.0,0.0,8))
#plt.xticks((0.96,0.97,0.98,0.99,1.0))
#plt.xlabel(r'$\rho_{tor}$',labelpad=10,size=14)
#plt.ylabel(r'$\eta_e$',size=14)
#plt.legend(loc='upper right')
#plt.tight_layout()
#plt.show()

print len(rhot1)
print len(n1)
print len(omt1)
f = open('profile_info_'+f_ge,'w')
f.write('#1.rhot 2.dummy 3.Te 4.ne 5.omte 6.omne 7.etae \n')
np.savetxt(f,np.column_stack((rhot1,rhot1,t1,n1,omt1,omn1,omt1/omn1)))
f.close()


f = open('profile_info_'+f_gi,'w')
f.write('#1.rhot 2.dummy 3.Ti 4.ni 5.omti 6.omni 7.etai \n')
np.savetxt(f,np.column_stack((rhot2,rhot2,t2,n2,omt2,omn2,omt2/omn2)))
f.close()

irhot1= np.argmin(abs(rhot1[:]  - rhot0))
irhot2= np.argmin(abs(rhot2[:]  - rhot0))

print 'omte at rho_tor = '+str(rhot0)+': ',omt1[irhot1]
print 'omti at rho_tor = '+str(rhot0)+': ',omt2[irhot2]
print 'omne at rho_tor = '+str(rhot0)+': ',omn1[irhot1]
print 'omni at rho_tor = '+str(rhot0)+': ',omn2[irhot2]

