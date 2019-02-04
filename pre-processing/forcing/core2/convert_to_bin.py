import xarray as xr
import numpy as np
import sys

def write_to_binary(ds, variable, fileout, precision='single'):
    ''' write variable from dataset ds to fileout with precision '''
    # write data to binary files
    fid   = open(fileout, "wb")
    flatdata = ds[variable].values.flatten()
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

core_dir='/local/data/artemis/workspace/rdussin/CORE2/by_year/'

listvar = ['Pair','Qair','Tair','Uwind','Vwind','lwrad','rain','swrad']
for year in np.arange(1948,2007+1):
    print('working on year', year)
    for var in listvar:
        print('working on variable', var)
        filein  = core_dir + var + '_CORE2_y' + str(year) + '.nc'
        fileout = core_dir + 'binaries/' + var + '_CORE2_y' + str(year)

        if var == 'lwrad':
            var = 'lwrad_down'

        ds = xr.open_dataset(filein, decode_times=False)

        # unit conversion
        if var == 'Tair':
            ds['Tair'] = 273.15 + ds['Tair'] # from degC to K

        if var == 'rain':
            ds['rain'] = ds['rain'] / 1000. # [kg.m-2.s-1] / [kg/m3] = [m/s]

        write_to_binary(ds, var, fileout, precision='single')
