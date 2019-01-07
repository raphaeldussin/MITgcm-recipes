#!/usr/bin/env python

import sys
import numpy as np
import xarray as xr

filein=sys.argv[-1]
fname=filein.replace('/',' ').split()[-1]
var=fname.replace('_',' ').split()[1]

fileout= fname + '.nc'

print('loading file', filein)

data = np.fromfile(filein,dtype='>f')
nx = 640 ; ny = 320
data = np.reshape(data, [-1,ny,nx])

dsgrid = xr.open_dataset('jra_grid.nc')

lon = dsgrid['longitude'].values
lat = dsgrid['latitude'].values
nt = data.shape[0]
time = np.arange(nt)

print('creating file', fileout)
print('min/max', data.min() , data.max() )

da = xr.DataArray(data, coords=[time, lat, lon], dims=['time', 'lat', 'lon'])
ds = xr.Dataset({var: (['time', 'lat', 'lon'],  da)},
                 coords={'lon': (['x', 'y'], lon),
                         'lat': (['x', 'y'], lat)})
ds.to_netcdf(fileout)
