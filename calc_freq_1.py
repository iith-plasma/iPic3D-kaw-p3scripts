import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

###########################User defined parameters##########################
exptdpks = 4
tini=0
tfin=1500
directory = "../3D_runs/120_runs/341_Lz20_mb/data/"
dt = 0.5
nsmooth=30
############################################################################

def find_peaks(arr):
    nelts = np.size(arr)
    buff = 5
    peaks = np.array([])
    #pack the array
    arr = np.append(arr,arr[0:buff])
    arr = np.append(arr[nelts-buff:nelts],arr)
    for uu in range (buff,nelts+buff):
        if ((arr[uu]>=arr[uu-1]) and (arr[uu]>=arr[uu+1])):
            if ((arr[uu]>arr[uu-2]) and (arr[uu]>arr[uu+2])):
                if ((arr[uu]>arr[uu-3]) and (arr[uu]>arr[uu+3])):
                    if ((arr[uu]>arr[uu-4]) and (arr[uu]>arr[uu+4])):
                        peaks=np.append(peaks,uu-buff)
    return(peaks)

def find_valleys(arr):
    nelts = np.size(arr)
    buff = 5
    valleys = np.array([])
    #pack the array
    arr = np.append(arr,arr[0:buff])
    arr = np.append(arr[nelts-buff:nelts],arr)
    for uu in range (buff,nelts+buff):
        if ((arr[uu]<=arr[uu-1]) and (arr[uu]<=arr[uu+1])):
            if ((arr[uu]<arr[uu-2]) and (arr[uu]<arr[uu+2])):
                if ((arr[uu]<arr[uu-3]) and (arr[uu]<arr[uu+3])):
                    if ((arr[uu]<arr[uu-4]) and (arr[uu]<arr[uu+4])):
                        valleys=np.append(valleys,uu-buff)
    return(valleys)


def smooth(arr):
    newarr = 0.25*np.roll(arr,-1)+0.5*arr+0.25*np.roll(arr,1)
    return(newarr)

filename  = directory+"fields"+str(tini).zfill(8)+"-0.vtr"
print("filename= ",filename)
Fields = pv.read(filename)
nd = Fields.GetDataDimension()
xx = Fields.GetDimensions()
xlen = xx[0]
ylen = xx[1]
zlen = xx[2]
E = Fields['E']
Ei = np.reshape(E,(zlen,ylen,xlen,3))
B = Fields['B']
Bi = np.reshape(B,(zlen,ylen,xlen,3))


filename  = directory+"fields"+str(tfin).zfill(8)+"-0.vtr"
print("filename= ",filename)
Fields = pv.read(filename)
nd = Fields.GetDataDimension()
xx = Fields.GetDimensions()
xlen = xx[0]
ylen = xx[1]
zlen = xx[2]
E = Fields['E']
Ef = np.reshape(E,(zlen,ylen,xlen,3))
B = Fields['B']
Bf = np.reshape(B,(zlen,ylen,xlen,3))

filename = directory+"currents"+str(tini).zfill(8)+"-0.vtr"
Currents = pv.read(filename)
Jfield = Currents['J0']
Je = np.reshape(Jfield,(zlen,ylen,xlen,3)) 
Jex = Je[:,:,:,0]
Jey = Je[:,:,:,1]
Jez = Je[:,:,:,2]
Jfield = Currents['J1']
Ji = np.reshape(Jfield,(zlen,ylen,xlen,3))
Jix = Ji[:,:,:,0]
Jiy = Ji[:,:,:,1]
Jiz = Ji[:,:,:,2]

Ji_ini = Je

filename = directory+"currents"+str(tfin).zfill(8)+"-0.vtr"
Currents = pv.read(filename)
Jfield = Currents['J0']
Je = np.reshape(Jfield,(zlen,ylen,xlen,3)) 
Jex = Je[:,:,:,0]
Jey = Je[:,:,:,1]
Jez = Je[:,:,:,2]
Jfield = Currents['J1']
Ji = np.reshape(Jfield,(zlen,ylen,xlen,3))
Jix = Ji[:,:,:,0]
Jiy = Ji[:,:,:,1]
Jiz = Ji[:,:,:,2]

Ji_fin = Je

#signal = Bi[:,3,0,0]    
#for i in range (0,100):
#    signal = smooth(signal)
#plt.plot(signal)
#plt.show()

freq_gather = 0.0
count = 0
for slc in range(0,zlen):

    signal0 = Bi[slc,:,0,0]
    for i in range (0,nsmooth):
        signal0 = smooth(signal0)

    signal1 = Bf[slc,:,0,0]
    for i in range (0,nsmooth):
        signal1 = smooth(signal1)

    """
    plt.plot(signal0, linestyle = 'solid')
    plt.plot(signal1, linestyle = "--")
    plt.savefig("test.eps")
    plt.show()
    """

    pp0=find_peaks(signal0)
    if (np.size(pp0) != exptdpks) :
        print("not getting correct number of peaks, skipping, np.size(pp0) ",slc)
        continue

    vv0=find_valleys(signal0)
    if (np.size(vv0) != exptdpks) :
        print("not getting correct number of valleys, skipping, no.size(vv0)",slc)
        continue

    pp1=find_peaks(signal1)
    if (np.size(pp1) != exptdpks) : 
        print("not getting correct number of peaks, skipping, np.size(pp1)", slc)
        continue

    vv1=find_valleys(signal1)
    if (np.size(vv1) != exptdpks) :
        print("not getting correct number of valleys, skipping, np.size(vv1)", slc)
        continue

    if (pp1[0]<pp0[0]) :
        pp1[0]=pp1[0]+ylen
        pp1 = np.roll(pp1,-1)

    if (vv1[0]<vv0[0]) :
        vv1[0]=vv1[0]+ylen
        vv1 = np.roll(vv1,-1)

    for p in range (0,exptdpks):
        distp = pp1[p]-pp0[p]
        distv = vv1[p]-vv0[p] 
        dist = 0.5*(distp+distv)
        freq = 2*np.pi*dist*exptdpks/(ylen*(tfin-tini)*dt)
        freq_gather = freq_gather + freq
        count = count+1

print("count, final freq= ", count, freq_gather/count)

"""
Fields = pv.read('../kaw3/data/fields00000500-0.vtr')

nd = Fields.GetDataDimension()
xx = Fields.GetDimensions()

xlen = xx[0]
ylen = xx[1]
zlen = xx[2]

E = Fields['E']
Efield = np.reshape(E,(zlen,ylen,xlen,3))

B = Fields['B']
B500 = np.reshape(B,(zlen,ylen,xlen,3))

"""

"""
Rho = pv.read('../simulations/test8/data/rho00000100-0.vtr')
rho0 = Rho['rho0']
ne = np.reshape(rho0,(zlen,ylen,xlen))
rho1 = Rho['rho1']
ni = np.reshape(rho1,(zlen,ylen,xlen))

Curr = pv.read('../simulations/test8/data/currents00000100-0.vtr')
J0 = Curr['J0']
je = np.reshape(J0,(zlen,ylen,xlen,3))
J1 = Curr['J1']
ji = np.reshape(J1,(zlen,ylen,xlen,3))
"""

"""
fig = plt.figure(figsize=(5, 5))
gs = gridspec.GridSpec(1, 1, hspace=0.1, wspace=0.2)
ax0 = plt.subplot(gs[0])
temp = Bfield[:,:,0,0]
z_min, z_max = temp.min(), temp.max()
plt.imshow(temp, cmap='rainbow', vmin=z_min, vmax=z_max,
           interpolation='nearest', origin='lower',extent=[-5.0, 5.0, -5.0, 5.0])
#ax0.set_xlabel('x(di)')
ax0.set_ylabel('z(di)')
ax0.set_title('vx')
cbar = plt.colorbar(pad=0.1)
plt.savefig('test.eps', transparent=False)
plt.show() 
"""

"""
def smooth(arr):
    newarr = 0.25*np.roll(arr,-1)+0.5*arr+0.25*np.roll(arr,1)
    return(newarr)

signal0 = B0[:,0,0,0]
for i in range (0,10):
    signal0 = smooth(signal0)

signal1 = B500[:,0,0,0]
for i in range (0,10):
    signal1 = smooth(signal1)

peak1 = np.argmax(signal0)
peak2 = np.argmax(signal1)
if (peak2 < peak1) :
    peak2 = peak2+(zlen-1)
print("peaks = ", peak1, peak2)
dist = peak2-peak1
freq = 2*np.pi*dist/((zlen-1)*500.0*0.5)
print("freq try= ", freq)

freqavg = 0.0
for k in range (ylen):

    signal0 = B0[:,k,0,0]
    for i in range (0,10):
        signal0 = smooth(signal0)

    signal1 = B500[:,k,0,0]
    for i in range (0,10):
        signal1 = smooth(signal1)

    peak1 = np.argmax(signal0)
    peak2 = np.argmax(signal1)
    if (peak2 < peak1) :
        peak2 = peak2+zlen-1

    dist = peak2-peak1
    freq = 2*np.pi*dist/((zlen-1)*500.0*0.5)
    print("freq= ", freq)
    freqavg = freqavg + freq

print("average fereq = ", freqavg/ylen)

    

plt.plot(signal0, linestyle = 'solid')
plt.plot(signal1, linestyle = 'dashed')
plt.show()

"""
