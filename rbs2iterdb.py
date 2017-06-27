#! /usr/bin/python

from read_EFIT import *
from read_iterdb_file import *
#from calc_fields_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys
from write_iterdb import *
from readRbsProfs import *

EFIT_file_name = sys.argv[1]
rbsProfs_fileName = sys.argv[2]
er_file_name = sys.argv[3]

zave = 10.
zeff = 2.2

EFITdict = read_EFIT(EFIT_file_name)
RBSdict = readRbsProfs(rbsProfs_fileName)

for quantity in ['R', 'qpsi', 'Pres', 'Bpol', 'Btor']:
    plt.plot(EFITdict['rhotn'], EFITdict[quantity], label = 'EFIT')
    plt.plot(RBSdict['rhotn'], RBSdict[quantity], label = 'rbsProfs')
    plt.title('compare EFIT and rbsProfs')
    plt.xlabel('rhot')
    plt.ylabel(quantity)
    plt.legend()
    plt.show()

psipn = EFITdict['psipn']
rhotn = EFITdict['rhotn']

uni_rhot = np.linspace(rhotn[0], rhotn[-1], len(rhotn)*10)
#psip_fs = 0.98
#R_fs, Z_fs, B_pol, B_tor, B_tot, psipn_obmp, R_obmp, Bp_obmp, Bt_obmp = \
#                                          Bfields(EFIT_file_name,psip_fs)
psipn_unirhot = interp(rhotn, psipn, uni_rhot)
R_unirhot = interp(psipn, EFITdict['R'], psipn_unirhot)

rhot1 = RBSdict['rhotn']
ti1 = RBSdict['ti']
ne1 = RBSdict['ne']
ni1 = ne1 * (zave - zeff) / (zave - 1.)
nz1 = ne1 * (zeff - 1.) / (zave - 1.) / zave

ti_unirhot = interp(rhot1, ti1, uni_rhot)
ni_unirhot = interp(rhot1, ni1, uni_rhot)
pi_unirhot = ti_unirhot * ni_unirhot

uni_R = np.linspace(R_unirhot[0], R_unirhot[-1], len(R_unirhot))
pi_uR = interp(R_unirhot, pi_unirhot, uni_R)
ni_uR = interp(R_unirhot, ni_unirhot, uni_R)
rhot_uR = interp(R_unirhot, uni_rhot, uni_R)
Er_gradPi = - first_derivative(pi_uR, uni_R) / ni_uR
er = interp(rhot_uR, Er_gradPi, rhot1)

er1 = RBSdict['Er']

if 1 == 1:
    plt.plot(rhot1, -er1, label = 'rbsProfs')
    plt.plot(rhot_uR, Er_gradPi, label = 'grad Pi / ni / e')
    plt.plot(rhot1, er, label = 'grad Pi / ni / e out')
    plt.xlabel('rhot')
    plt.ylabel('V/m')
    plt.axis([0.9, 1.0, -100000., 300000.])
    plt.legend(loc = 2)
    plt.show()
if 1 == 0:
    f = open('rbsProfs_er_gradPi', 'w')
    f.write('#1.rhot 2.Er(V/m)\n')
    np.savetxt(f,np.column_stack((rhot1,er)))
    f.close()

erData0 = np.genfromtxt(er_file_name, dtype = 'float')
rhotEr0 = erData0[:, 0]
er6 = - erData0[:, 1]
er7 = interp(rhotEr0 - 0.011, er6, rhot1)
er8 = interp(rhotEr0 - 0.005, er6, rhot1)
#Er_exp = interp(rEr, er2, R_unirhot)
#er3 = interp(uni_rhot, Er_exp, rhot1)
#er4 = interp(uni_rhot - 0.011, Er_exp, rhot1)
#er5 = interp(uni_rhot - 0.005, Er_exp, rhot1)

if 1 == 1:
    R1 = RBSdict['R']
    Bpol1 = RBSdict['Bpol']
    plt.plot(rhot1, er8/ R1/ Bpol1, label = 'exp er8 out')
    plt.plot(rhot1, er/ R1/ Bpol1, label = 'grad Pi / ni / e out')
    plt.plot(rhot1, - er1/ R1/ Bpol1, label = 'rbsProfs')
    plt.xlabel('rhot')
    plt.ylabel('rad/s')
    plt.legend(loc = 2)
    plt.axis([0.9, 1.0, -100000., 400000.])
    plt.show()
if 1 == 0:
    f = open('rbsProfs_er_exp_left011', 'w')
    f.write('#1.rhot 2.Er(V/m)\n')
    np.savetxt(f,np.column_stack((rhot1,er8)))
    f.close()
omega_tor = er / R1 / Bpol1
omega_tor1 = er7 / R1 / Bpol1
if 1 == 1:
    file_out_base = 'CmodH2_n48_allTi' 
    base_number = 'gradPi' + 'test'
    time_str = '01075'
    #input T in KeV, n in 10^19 m^-3
    output_iterdb(rhot1,RBSdict['psipn'],ne1*10.,RBSdict['te']*1.E-3,ni1*10.,ti1*1.E-3,file_out_base+base_number,base_number,time_str,vrot=omega_tor,nimp=nz1*10.)
    file_out_base = 'CmodH2_n48_allTi' 
    base_number = 'ExpEr' + 'test'
    time_str = '01075'
    #input T in KeV, n in 10^19 m^-3
    output_iterdb(rhot1,RBSdict['psipn'],ne1*10.,RBSdict['te']*1.E-3,ni1*10.,ti1*1.E-3,file_out_base+base_number,base_number,time_str,vrot=omega_tor1,nimp=nz1*10.)
