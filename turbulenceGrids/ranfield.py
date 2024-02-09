import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.size'] = 18.0
import matplotlib.pyplot as plt
from grf import vector_grf

nx = 80
nd = (nx,nx,nx)
beta = 2.

vf = vector_grf(nd,beta,field_type='solenoidal')

velfield = vf.signal.real

gridsize = velfield.shape[0]

vfile = open('solenoidalVel3D.bin','w')

vfile.write('{:d}'.format(gridsize))

vfile.write('\n')

for i in range(gridsize):
    for j in range(gridsize):
        for k in range(gridsize):
            line = ''
            for n in range(3):
                line += '{:.4e}'.format(velfield[i,j,k,n])
                line += ' '
            line += '\n'
            vfile.write(line)

vfile.close()

