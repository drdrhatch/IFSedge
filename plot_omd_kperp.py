#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
import optparse as op
from read_write_geometry import *
from subprocess import call
import sys
import math
from finite_differences import *
import scipy.special as ss

parser=op.OptionParser(description='Plots eigenfunction, kperp and omega_d.')
parser.add_option('--plot_phi','-p',action='store_const',const=1,help='plot',default=0)
parser.add_option('--plot_phi_kz','-k',action='store_const',const=1,help='plot',default=0)
options,args=parser.parse_args()
suffix = args[0]
plot_phi = options.plot_phi
plot_phi_kz = options.plot_phi_kz

suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

if 'Lref' in pars:
   Lref = pars['Lref']
else:
   Lref = 1.
if 'Bref' in pars:
   Bref = pars['Bref']
else:
   Bref = 1.
if 'nref' in pars:
   nref = pars['nref']
else:
   nref = 1.
if 'Tref' in pars:
   Tref = pars['Tref']
else:
   Tref = 1.

omn = float(raw_input('Enter omn: \n'))
omt = float(raw_input('Enter omt: \n'))
beta = 4.03E-3*nref*Tref/Bref**2

geom = pars['magn_geometry'][1:-1]

field = fieldfile('field'+suffix,pars)
field.set_time(field.tfld[-1])
print 't = ', field.tfld[-1]

if field.nx%2:
   zgrid = np.linspace(-field.nx,field.nx,field.nx*field.nz,endpoint=False)
else:
   zgrid = np.linspace(-(field.nx-1),(field.nx+1),field.nx*field.nz,endpoint=False)
dz = zgrid[1]-zgrid[0]

ikx_grid = np.arange(-field.nx/2+1,field.nx/2+1)
theta_grid = np.zeros(field.nx*field.nz,dtype='float128')
if 'edge_opt' in pars:
	theta_even = np.linspace(-np.pi,np.pi,field.nz,endpoint=False)
	N = np.arcsinh(pars['edge_opt']*theta_even[0])/theta_even[0]
	theta = 1./pars['edge_opt']*np.sinh(N*theta_even)
	for i in ikx_grid:
		this_theta_grid = i*2+theta/np.pi
		theta_grid[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_theta_grid

kperp = np.zeros(field.nx*field.nz,dtype='float128')
omega_d = np.zeros(field.nx*field.nz,dtype='float128')
phi = np.zeros(field.nx*field.nz,dtype='complex128')
phi2 = np.zeros(field.nx*field.nz,dtype='complex128')
Jac = np.zeros(field.nx*field.nz,dtype='complex128')


if 'x_local' in pars:
    if pars['x_local']:
        x_local = True
    else:
        x_local = False 
else:
    x_local = True

if 'kx_center' in pars:
    kx_center = pars['kx_center']
else:
    kx_center = 0.

if x_local:
	zgrid_c = np.linspace(-1,1,field.nz,endpoint=False)
	kymin = pars['kymin']
	lx = pars['lx']
	#geom_file = 'tracer_efit'+suffix
	geom_file = geom+suffix
	#dkx = 2*np.pi/lx
	dkx = 2*np.pi*pars['shat']*kymin
	gpars,geometry = read_geometry_local(geom_file)
	ggxx = geometry['ggxx']
	ggxy = geometry['ggxy']
	ggyy = geometry['ggyy']
	ggxz = geometry['ggxz']
	ggyz = geometry['ggyz']
	gdBdx = geometry['gdBdx']
	gdBdy = geometry['gdBdy']
	gdBdz = geometry['gdBdz']
	gBfield = geometry['gBfield']
	jacob = geometry['gjacobian']
	zprime_grid_c = zgrid_c*np.pi/gBfield/jacob

for i in ikx_grid:
	kx = i*dkx+kx_center
	this_kperp = np.sqrt(ggxx*kx**2+2.*ggxy*kx*kymin+ggyy*kymin**2)

	gamma1 = ggxx*ggyy-ggxy**2
	gamma2 = ggxx*ggyz-ggxy*ggxz
	gamma3 = ggxy*ggyz-ggyy*ggxz
	Kx = -(gdBdy+gamma2/gamma1*gdBdz)*Lref/Bref
	Ky = gdBdx-gamma3/gamma1*gdBdz*Lref/Bref
	omega_d_y = Ky*kymin
	omega_d_x = Kx*kx
	this_omegad =  -(omega_d_y + omega_d_x)/gBfield+beta*(omn+omt)*\
	kymin/2./gBfield**2

	kperp[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_kperp
	omega_d[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_omegad
	Jac[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=1./np.pi/jacob/gBfield

if 'n0_global' in pars:
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars['n0_global']*pars['q0'])
else:
        phase_fac = -1.0
    
if pars['shat'] > 0.:
    for i in ikx_grid:
	this_phi = field.phi()[:,0,i]*phase_fac**i
	phi[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_phi
else:
    for i in ikx_grid:
	this_phi = field.phi()[:,0,-i]*phase_fac**i
	phi[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_phi
phi = phi/field.phi()[field.nz/2,0,0]
phi = phi/max(np.abs(phi))

plt.plot(zgrid,np.real(phi),label=r'$Re[\phi]$',color='red')
plt.plot(zgrid,np.imag(phi),label=r'$Im[\phi]$',color='blue')
plt.plot(zgrid,np.abs(phi),label=r'$|\phi|$',color='black')
if pars['n_fields']==1:
	plt.plot(zgrid,kperp/60,label='kperp',color='lawngreen')
	plt.plot(zgrid,omega_d/10,label='omega_d',color='aqua')
else:
	plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
	plt.plot(zgrid,omega_d,label='omega_d',color='aqua')
plt.title('ky = '+str(pars['kymin'])+', kx_center = '+\
str(kx_center))
#plt.axis([-1.,1.,-6.,6.])
plt.grid()
plt.show()
