#!/usr/bin/env python
import xarray as xr
import numpy as np

def write_grid_info(ds, varname, fileout):
    ''' write formatted data for data.exf '''
    lon0 = ds.LON.values.min()
    lon_inc = ds.LON.values[1] - ds.LON.values[0]
    lat0 = ds.LAT.values.min()
    lat_inc = ds.LAT.values[1:] - ds.LAT.values[:-1]
    nlon = len(ds.LON.values)
    nlat = len(ds.LAT.values)

    line1 = ' ' + varname + '_lon0       = %.3fD0,' % (lon0) 
    line2 = ' ' + varname + '_lon_inc    = %.3fD0,' % (lon_inc)
    line3 = ' ' + varname + '_lat0       = %.6fD0,' % (lat0)
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

dir_core='/local/data/artemis/workspace/rdussin/CORE2-GFDL/'

fileout='data.exf_CORE2'

ds = xr.open_mfdataset(dir_core + 'ncar_precip.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'precip', fileout)

ds = xr.open_mfdataset(dir_core + 't_10.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'atemp', fileout)

ds = xr.open_mfdataset(dir_core + 'q_10.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'aqh', fileout)

ds = xr.open_mfdataset(dir_core + 'ncar_rad.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'swdown', fileout)

ds = xr.open_mfdataset(dir_core + 'ncar_rad.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'lwdown', fileout)

ds = xr.open_mfdataset(dir_core + 'u_10.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'uwind', fileout)

ds = xr.open_mfdataset(dir_core + 'v_10.1948-2009.23OCT2012.nc', decode_times=False)
write_grid_info(ds, 'vwind', fileout)

