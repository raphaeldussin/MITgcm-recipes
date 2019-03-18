import xarray as xr
import numpy as np
import sys
import calendar

def write_to_binary(ds, variable, fileout, precision='single'):
    ''' write variable from dataset ds to fileout with precision '''
    # write data to binary files
    fid   = open(fileout, "wb")
    flatdata = ds[variable].transpose(*('TIME','LAT','LON')).values.flatten()
    if precision == 'single':
        if sys.byteorder == 'little':
            tmp = flatdata.astype(np.dtype('f')).byteswap(True).tobytes()
        else:
            tmp = flatdata.astype(np.dtype('f')).tobytes()
        fid.write(tmp)
    elif precision == 'double':
        if sys.byteorder == 'little':
            tmp = flatdata.astype(np.dtype('d')).byteswap(True).tobytes()
        else:
            tmp = flatdata.astype(np.dtype('d')).tobytes()
    fid.close()
    return None

core_dir='/local/data/artemis/workspace/rdussin/CORE2-GFDL/yearly/'

# 6-hourly
listvar = ['Q_10_MOD','T_10_MOD','U_10_MOD','V_10_MOD']

for year in np.arange(1948,2009+1):
    print('6-hourly, working on year', year)
    filelist=[]
    for var in listvar:
        filelist.append(core_dir + var + '_CORE2_y' + str(year) + '.nc')

    dsyear = xr.open_mfdataset(filelist, decode_times=False)

    # add a second january 1st if leap year
    if calendar.isleap(year):
        ds = xr.concat([dsyear.isel(TIME=range(4)), dsyear], dim='TIME')
    else:
        ds = dsyear

    for var in listvar:
        print('working on variable', var)
        fileout = core_dir + 'binaries/' + var + '_CORE2_' + str(year)
        write_to_binary(ds, var, fileout, precision='single')

# daily
listvar = ['LWDN_MOD','SWDN_MOD']

for year in np.arange(1948,2009+1):
    print('daily, working on year', year)
    filelist=[]
    for var in listvar:
        filelist.append(core_dir + var + '_CORE2_y' + str(year) + '.nc')

    dsyear = xr.open_mfdataset(filelist, decode_times=False)

    # add a second january 1st if leap year
    if calendar.isleap(year):
        ds = xr.concat([dsyear.isel(TIME=0), dsyear], dim='TIME')
    else:
        ds = dsyear

    for var in listvar:
        print('working on variable', var)
        fileout = core_dir + 'binaries/' + var + '_CORE2_' + str(year)
        write_to_binary(ds, var, fileout, precision='single')

# monthly
listvar = ['RAIN','SNOW']

for year in np.arange(1948,2009+1):
    print('monthly, working on year', year)
    filelist=[]
    for var in listvar:
        filelist.append(core_dir + var + '_CORE2_y' + str(year) + '.nc')

    dsyear = xr.open_mfdataset(filelist, decode_times=False)
    # total precip + unit conversion
    dsyear['PRECIP'] = (dsyear['RAIN'] + dsyear['SNOW']) /1000. # [kg.m-2.s-1] / [kg/m3] = [m/s]

    fileout = core_dir + 'binaries/' + 'PRECIP' + '_CORE2_' + str(year)
    write_to_binary(dsyear, 'PRECIP', fileout, precision='single')

