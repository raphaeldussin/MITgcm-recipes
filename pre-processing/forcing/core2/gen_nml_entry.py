#!/usr/bin/env python
import xarray as xr
import numpy as np

def write_grid_info(ds, varname, fileout):
    ''' write formatted data for data.exf '''
    lon0 = ds.lon.values.min()
    lon_inc = ds.lon.values[1] - ds.lon.values[0]
    lat0 = ds.lat.values.min()
    lat_inc = ds.lat.values[1:] - ds.lat.values[:-1]
    nlon = len(ds.lon.values)
    nlat = len(ds.lat.values)

    line1 = ' ' + varname + '_lon0       = %.3fD0,' % (lon0) 
    line2 = ' ' + varname + '_lon_inc    = %.3fDO,' % (lon_inc)
    line3 = ' ' + varname + '_lat0       = %.6fDO,' % (lat0)
    line4 = ' ' + varname + '_lat_inc    = '
    ct=1
    for dlat in lat_inc:
        line4 = line4 + '%.6fD0, ' % (dlat) 
        if np.mod(ct,5) == 0:
            line4 = line4 + '\n                     '
        ct=ct+1
    line5 = ' ' + varname + '_nlon       = %3i' % (nlon)
    line6 = ' ' + varname + '_nlat       = %2i' % (nlat)

    fid = open(fileout,'a')
    fid.write(line1 + '\n')
    fid.write(line2 + '\n')
    fid.write(line3 + '\n')
    fid.write(line4 + '\n')
    fid.write(line5 + '\n')
    fid.write(line6 + '\n')
    fid.write('#\n')
    fid.close()
    return None

dir_core='/local/data/artemis/workspace/rdussin/CORE2/'

fileout='data.exf_CORE2'

ds = xr.open_mfdataset(dir_core + 'rain.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'precip', fileout)

ds = xr.open_mfdataset(dir_core + 'Tair.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'atemp', fileout)

ds = xr.open_mfdataset(dir_core + 'Qair.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'aqh', fileout)

ds = xr.open_mfdataset(dir_core + 'swrad.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'swdown', fileout)

ds = xr.open_mfdataset(dir_core + 'lwrad.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'lwdown', fileout)

ds = xr.open_mfdataset(dir_core + 'Uwind.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'uwind', fileout)

ds = xr.open_mfdataset(dir_core + 'Vwind.1948-2007.nc', decode_times=False)
write_grid_info(ds, 'vwind', fileout)

