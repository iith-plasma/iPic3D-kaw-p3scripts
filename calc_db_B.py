import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

###################User defined parameters#####################################
B0x = 0.0
B0y = 0.0
B0z = 0.01
Fields = pv.read('../3D_runs/120_runs/341_Lz20_mb/data/fields00000000-0.vtr')
###############################################################################

meanB0 = np.sqrt(B0x**2+B0y**2+B0z**2)
nd = Fields.GetDataDimension()
xx = Fields.GetDimensions()

xlen = xx[0]
ylen = xx[1]
zlen = xx[2]

E = Fields['E']
Efield = np.reshape(E,(zlen,ylen,xlen,3))

B = Fields['B']
Bfield = np.reshape(B,(zlen,ylen,xlen,3))

dBx = Bfield[:,:,:,0]-B0x
dBy = Bfield[:,:,:,1]-B0y
dBz = Bfield[:,:,:,2]-B0z

dbrms = np.sqrt(np.sum(dBx**2+dBy**2+dBz**2)/(xlen*ylen*zlen))

print("db_B= ",dbrms/meanB0)



