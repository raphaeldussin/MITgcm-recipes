import numpy as np
import sys
import matplotlib.pylab as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap 

# from SIZE.h, we read that :
sNx =  46
sNy =  34
nSx =   1
nSy =   1
nPx =   4
nPy =   6

# which leads to :
nx = sNx * nSx * nPx
ny = sNy * nSy * nPy

# read bathy file and reshape
dirin = '/Users/raphael/WORK/preprocessing_mitgcm/'
filein = dirin + 'bathy_from_Depth.bin'
raw = np.fromfile(filein, dtype='>f')
bathy = np.reshape(raw, (ny,nx))
bathy = np.ma.masked_values(bathy,0)

# plot in i,j coordinates
x = np.arange(nx)
y = np.arange(ny)
x2,y2 = np.meshgrid(x,y)
plt.figure() ; plt.pcolormesh(x2,y2,bathy,cmap=cm.jet)
plt.plot(x2[::sNy,:],y2[::sNy,:],'k.')
plt.plot(x2[:,::sNx],y2[:,::sNx],'k.')
plt.colorbar()
plt.title('north atlantic model and processor layout')

# from MITgcm data (= namelist) parameters, we read:
phiMin=-20.25
thetaMin=-72.25
dX=0.5 # from delX, take just the step in x
dY=0.5

theta = thetaMin + dX * np.arange(nx)
phi = phiMin + dY * np.arange(ny)

# plot in geographical coordinates :

def setup_map(plt_topo=False,hide_grid=False):
        ''' set the map, with etopo optionally '''
        # background
        bmap = Basemap(projection='cyl',llcrnrlat=-22,urcrnrlat=85,\
                                         llcrnrlon=-105,urcrnrlon=35,resolution='l')
        parallels = np.arange(-20.,80.,20.)
        meridians = np.arange(-120.,60.,30.)
        if hide_grid:
                bmap.drawparallels(parallels,labels=[True,False,False,True],linewidth=1,color=[1.,1.,1.])
                bmap.drawmeridians(meridians,labels=[True,False,False,True],linewidth=1,color=[1.,1.,1.])
        else:
                bmap.drawparallels(parallels,labels=[True,False,False,True],linewidth=1,color=[0.6,0.6,0.6],fontsize=16)
                bmap.drawmeridians(meridians,labels=[True,False,False,True],linewidth=1,color=[0.6,0.6,0.6],fontsize=16)
        bmap.drawcoastlines()
        if plt_topo:
		bmap.bluemarble()
        return bmap

plt.figure()
m = setup_map()
x2,y2 = np.meshgrid(theta,phi)
m.contourf(x2,y2,bathy,cmap=cm.jet)
plt.colorbar()
plt.title('north atlantic model as computed')

# bring back the Gulf of Mexico where it belongs:
x2[55:125,132:] = x2[55:125,132:] - ((nx-1) * 0.5)

plt.figure()
m = setup_map(plt_topo=True)
m.contourf(x2,y2,bathy,cmap=cm.jet)
plt.colorbar()
plt.title('north atlantic model with Gulf at its geographical location')
plt.show()
