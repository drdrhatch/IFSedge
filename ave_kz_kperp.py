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
from finite_differences_x import *
from finite_differences import *

parser=op.OptionParser(description='Plots eigenfunction, kperp and omega_d.')
options,args=parser.parse_args()
suffix = args[0]
#geom = args[1]
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
geom = pars['magn_geometry'][1:-1]

if geom == 's_alpha' or 'miller':
	z0 = 0.
	kz0 = 0.
elif geom == 'tracer_efit':
	z0 = float(raw_input('Enter offset z0: \n'))
	kz0 = float(raw_input('Enter offset kz0: \n'))

field = fieldfile('field'+suffix,pars)
field.set_time(field.tfld[-1])
print 't = ', field.tfld[-1]

if field.nx%2:
   zgrid = np.linspace(-field.nx,field.nx,field.nx*field.nz,endpoint=False)

else:
   zgrid = np.linspace(-(field.nx-1),(field.nx+1),field.nx*field.nz,endpoint=False)

ikx_grid = np.arange(-field.nx/2+1,field.nx/2+1)
theta_grid = np.zeros(field.nx*field.nz,dtype='float128')
if 'edge_opt' in pars:
	theta_even = np.linspace(-np.pi,np.pi,field.nz,endpoint=False)
	N = np.arcsinh(pars['edge_opt']*theta_even[0])/theta_even[0]
	theta = 1./pars['edge_opt']*np.sinh(N*theta_even)
	#plt.plot(theta_even,'.',label='theta_even')
	#plt.plot(theta,'.',label='theta')
	#plt.legend()
	#plt.show()
	for i in ikx_grid:
		this_theta_grid = i*2+theta/np.pi
		theta_grid[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_theta_grid
kperp = np.zeros(field.nx*field.nz,dtype='float128')
omega_d = np.zeros(field.nx*field.nz,dtype='float128')
phi = np.zeros(field.nx*field.nz,dtype='complex128')
phi2 = np.zeros(field.nx*field.nz,dtype='complex128')
B0 = np.zeros(field.nx*field.nz,dtype='float128')
Jac = np.zeros(field.nx*field.nz,dtype='float128')


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
	this_omegad =  -(omega_d_y + omega_d_x)/gBfield

	kperp[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_kperp
	omega_d[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=this_omegad
	B0[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=gBfield
	Jac[(i-ikx_grid[0])*field.nz:(i-ikx_grid[0]+1)*field.nz]=jacob


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
phi = phi/max(np.abs(phi))*10.




if 'edge_opt' in pars:
	plt.plot(theta_grid,np.real(phi),label=r'$Re[\phi]$',color='red')
	plt.plot(theta_grid,np.imag(phi),label=r'$Im[\phi]$',color='blue')
	plt.plot(theta_grid,np.abs(phi),label=r'$|\phi|$',color='black')
	plt.plot(theta_grid,kperp,label='kperp',color='lawngreen')
	plt.plot(theta_grid,omega_d,label='omega_d',color='aqua')
	plt.xlabel('z/pi')
	#plt.legend()
	plt.title('ky = '+str(kymin))
	plt.grid()
	plt.show()

	#phi[:len(phi)/2] = 0.
	phi[len(phi)/2:] = 0.
	i_phi = 0.
	denom = 0.
	for i in np.arange(len(phi)-1):
    	    denom = denom + (abs(phi[i])+abs(phi[i+1]))/2.*(theta_grid[i+1]-theta_grid[i])
    	    i_phi = i_phi + ((theta_grid[i]-z0)**2*abs(phi[i])+\
		(theta_grid[i+1]-z0)**2*abs(phi[i+1]))/2.*(theta_grid[i+1]-theta_grid[i])
	i_phi = np.sqrt(i_phi/denom*2.)
	print 'eigenfunction width (pi)', i_phi

	ave_kperp = 0.
	ave_omegad = 0.
	denom = 0.
	for i in np.arange(len(phi)-1):
		ave_kperp = ave_kperp + (kperp[i]**2*abs(phi[i])**2 +\
			kperp[i+1]**2*abs(phi[i+1])**2)/2.*\
			(theta_grid[i+1]-theta_grid[i])
		ave_omegad = ave_omegad + (omega_d[i]*abs(phi[i])**2 +\
			omega_d[i+1]*abs(phi[i])**2)/2.*\
			(theta_grid[i+1]-theta_grid[i])
		denom = denom + (abs(phi[i])**2 +abs(phi[i+1])**2)/2.*\
			(theta_grid[i+1]-theta_grid[i])
	ave_kperp = np.sqrt(ave_kperp/denom)
	print 'weighted k_perp', ave_kperp,\
		'(gyy', (ave_kperp/pars['kymin'])**2,')'

	ave_omegad = ave_omegad/denom
	print 'weighted omega_d', ave_omegad,\
		'(omd', (ave_omegad/pars['kymin']),')'
		#'(omd', (ave_omegad/pars['kymin']*pars['major_R']),')'

	zi=complex(0,1)
	phi_kz = np.empty(0,dtype='complex128')
	nkz = 500
	lkz = field.nz/2
	kz_grid = np.linspace(-lkz,lkz,500,endpoint=False)
	#phi[:field.nz*field.nx/2] = 0.
	#phi[field.nz*field.nx/2:] = 0.
	for k in np.arange(len(kz_grid)):
    		this_phi_kz = 0.
    		for i in np.arange(len(zgrid)-1):
        	    this_phi_kz = this_phi_kz + 0.5*(phi[i]*\
 		    np.exp(zi*kz_grid[k]*theta_grid[i])\
        	    +phi[i+1]*np.exp(zi*kz_grid[k]*theta_grid[i+1]))*\
        	    (theta_grid[i+1]-theta_grid[i])
    		phi_kz = np.append(phi_kz,this_phi_kz)
	plt.plot(kz_grid,np.abs(phi_kz),label='abs(phi_kz)')
	plt.plot(kz_grid,np.real(phi_kz),label='real(phi_kz)')
	plt.plot(kz_grid,np.imag(phi_kz),label='imag(phi_kz)')
	plt.xlabel('kz')
	plt.title('ky = '+str(kymin))
	plt.legend()
	plt.show()
	ave_kz = np.sqrt(sum((kz_grid-kz0)**2*abs(phi_kz)**2)/sum(abs(phi_kz)**2))
	print 'weighted kz', ave_kz
	print 'disp rel kz', ave_kz/np.pi/pars['q0']/pars['major_R']

else:
	ave_kperp = np.sqrt(sum(kperp**2*abs(phi)**2)/sum(abs(phi)**2))
	print 'weighted k_perp', ave_kperp,\
		'(gyy', (ave_kperp/pars['kymin'])**2,')'

	ave_omegad = sum(omega_d*abs(phi)**2)/sum(abs(phi)**2)
	print 'weighted omega_d', ave_omegad,\
		'(omd', (ave_omegad/pars['kymin']),')'
		#'(omd', (ave_omegad/pars['kymin'])*pars['major_R'],')'
	phi_gaussian = np.exp(-zgrid**2)*zgrid
	phi_gaussian = phi_gaussian/max(phi_gaussian)
	#plt.plot(zgrid,np.abs(phi_gaussian),label=r'$Re[\phi]$',color='green')
	plt.plot(zgrid,np.real(phi),label=r'$Re[\phi]$',color='red')
	plt.plot(zgrid,np.imag(phi),label=r'$Im[\phi]$',color='blue')
	plt.plot(zgrid,np.abs(phi),label=r'$|\phi|$',color='black')
	plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
	plt.plot(zgrid,omega_d,label='omega_d',color='aqua')
	plt.title('shat = '+str(pars['shat']))
	#plt.legend()
	#plt.axis([-3.,3.,-.5,1.5])
	plt.grid()
	plt.show()

	#For eigenfunction with two bumps
	#phi[:len(phi)/2] = 0.
	#phi[len(phi)/2:] = 0.
	i_phi = 0.
	denom = 0.
	for i in np.arange(len(phi)-1):
    	    denom = denom + (abs(phi[i])+abs(phi[i+1]))/2.*(zgrid[i+1]-zgrid[i])
    	    i_phi = i_phi + ((zgrid[i]-z0)**2*abs(phi[i])+(zgrid[i+1]-z0)**2*abs(phi[i+1]))/2.*(zgrid[i+1]-zgrid[i])
	i_phi = np.sqrt(i_phi/denom/2.)
	print 'eigenfunction width (pi)', i_phi

	zi=complex(0,1)
	phi_kz = np.empty(0,dtype='complex128')
	nkz = 500
	lkz = field.nz/2
	kz_grid = np.linspace(-lkz,lkz,500,endpoint=False)
	for k in np.arange(len(kz_grid)):
    		this_phi_kz = 0.
    		for i in np.arange(len(zgrid)-1):
        	    this_phi_kz = this_phi_kz + 0.5*(phi[i]*\
 		    np.exp(zi*kz_grid[k]*zgrid[i])\
        	    +phi[i+1]*np.exp(zi*kz_grid[k]*zgrid[i+1]))*\
        	    (zgrid[i+1]-zgrid[i])
    		phi_kz = np.append(phi_kz,this_phi_kz)
	plt.plot(kz_grid,np.abs(phi_kz),label='abs(phi_kz)')
	plt.plot(kz_grid,np.real(phi_kz),label='real(phi_kz)')
	plt.plot(kz_grid,np.imag(phi_kz),label='imag(phi_kz)')
	plt.xlabel('kz')
	plt.title('ky = '+str(kymin))
	plt.legend()
	plt.show()
	ave_kz = np.sqrt(sum((kz_grid-kz0)**2*abs(phi_kz)**2)/sum(abs(phi_kz)**2))
	print 'ave_kz = ', ave_kz
	print 'd.r. kz = ', ave_kz/np.pi/pars['q0']/pars['major_R']

	#dphidz = first_derivative(phi,zgrid)
	dphidz = fd_d1_o4(phi,zgrid)/np.pi/B0/Jac
	plt.plot(zgrid,np.abs(dphidz),label='abs(dphi/dz)')
	plt.plot(zgrid,np.real(dphidz),label='real(dphi/dz)')
	plt.plot(zgrid,np.imag(dphidz),label='imag(dphi/dz)')
	plt.legend()
	plt.xlabel('z')
	plt.show()
	ave_kz2 = np.sqrt(sum(abs(dphidz)**2)/sum(abs(phi)**2))
	print 'ave_kz2 = ', ave_kz2
#	print 'd.r. kz2 = ', ave_kz2/np.pi/pars['q0']/pars['major_R']
