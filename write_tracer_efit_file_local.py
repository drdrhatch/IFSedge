#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import sys

def write_tracer_efit_file_local(parameters,geometry,file_name):
    f = open(file_name,'w')
    f.write('&parameters\n')
    f.write('gridpoints =    '+str(int(float(parameters['gridpoints'])))+'\n')
    f.write('q0 =    '+str(parameters['q0'])+'\n')
    f.write('shat =    '+str(parameters['shat'])+'\n')
    f.write('s0 =    '+str(parameters['s0'])+'\n')
    f.write('minor_r =    '+str(parameters['minor_r'])+'\n')
    f.write('major_R =    '+str(parameters['major_R'])+'\n')
    f.write('trpeps =    '+str(parameters['trpeps'])+'\n')
    f.write('beta =    '+str(parameters['beta'])+'\n')
    f.write('Lref =    '+str(parameters['Lref'])+'\n')
    f.write('Bref =    '+str(parameters['Bref'])+'\n')
    f.write('magn_geometry =    '+'\''+str(parameters['magn_geometry'])+'\'\n/\n')

    nz0 = int(float(parameters['gridpoints']))

    geom_list = ['ggxx','ggxy','ggxz','ggyy','ggyz','ggzz','gBfield','gdBdx','gdBdy','gdBdz','gjacobian','gl_R','gl_phi','gl_z','gl_dxdR','gl_dxdZ']
    for i in range(nz0):
        for j in geom_list:
            f.write("%20.12E"% geometry[j][i])
        f.write('\n')
   
    f.close()


