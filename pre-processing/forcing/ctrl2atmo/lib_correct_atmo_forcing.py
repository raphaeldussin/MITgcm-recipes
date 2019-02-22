import xarray as xr
import xmitgcm 
import pyresample as _pyresample
import numpy as np
import pandas as pd
import struct, sys
import MITgcm_recipes

def read_forcing_grid(filein, lonvar='lon', latvar='lat'):
    ''' read coordinates of forcing grid '''
    ds = xr.open_dataset(filein)
    lon = ds[lonvar].values
    lat = ds[latvar].values

    return lon, lat

def read_forcing_from_binary(filein, varname, lon, lat, 
                             dateref='1900-1-1', startdate='1980-1-1',
                             freq='3H', precision='single'):
    ''' read forcing file from binary and return dataset '''

    if precision == 'single':
       data = np.fromfile(filein,'>f')
    else:
       data = np.fromfile(filein,'>d')
    
    ny = len(lat)
    nx = len(lon)

    data = np.reshape(data, (-1, ny, nx)) # nt variable (leap years)
    nt = data.shape[0]
    ds = xr.Dataset({varname: (['time', 'lat', 'lon'], data)},
                 coords={'lon': (['lon'], lon),
                         'lat': (['lat'], lat),
                         'time': pd.date_range(startdate, freq=freq, periods=nt),
                         'reference_time': pd.Timestamp(dateref)})
    return ds
       

def read_ctrl_file(varname, optim_cycle, ntimes, grid_dir, refdate='2002-01-01', ncname=None): 
    ''' read the content of the XX file '''
    astemd = xmitgcm.utils.get_extra_metadata(domain='aste',nx=270)

    tmp = xmitgcm.utils.read_mds(varname, iternum=optim_cycle, use_mmap=False,
                                  endian='>', dtype=np.dtype('f'), use_dask=True, 
                                  extra_metadata=astemd, chunks='2D', llc=True, legacy=True,
                                  shape=(ntimes,astemd['ny'],astemd['nx']))
    
    grid = xmitgcm.open_mdsdataset(grid_dir, prefix=['XC','YC'], geometry='llc', 
                                   extra_metadata=astemd, nx=astemd['nx'])

    #print('reading', tmp[list(tmp.keys())[0]])
    data = xr.DataArray(tmp[list(tmp.keys())[0]], dims=['time', 'face', 'j', 'i'],
                        coords={
                             'time': pd.date_range(refdate, freq='14D', periods=419),
                             'reference_time': pd.Timestamp(refdate)})

    data.name = ncname

    return data, grid

def regrid_to_forcing(lon_atm, lat_atm, grid_ctrl, data_ctrl,
                      varname, refdate='2002-01-01', freq='14D'):
    ''' regrid from model grid to atmospheric grid '''

    nt = data_ctrl.shape[0]
    nyout = lat_atm.shape[0]
    nxout = lon_atm.shape[0]

    # input grid
    lon_ctrl = grid_ctrl['XC'].values.ravel()
    lat_ctrl = grid_ctrl['YC'].values.ravel()
    # lon_ctrl = np.mod(lon_ctrl+360,360) does not work, see below

    grid_src = _pyresample.geometry.SwathDefinition(lons=lon_ctrl, 
                                                    lats=lat_ctrl)

    # output grid (but shifted to -180/180 to fit ASTE)
    # Note: for some unexplained reason, moving ASTE to 0-360 does not do the job
    # and has all values between 180/360 not computed
    lons_atm, lats_atm = np.meshgrid(lon_atm, lat_atm)
    grid_target = _pyresample.geometry.GridDefinition(lons=lons_atm-180.,
                                                      lats=lats_atm)

    # we need to be able to reset the data on the 0-360 forcing grid
    ind_greenwich= (np.abs(lon_atm-180.)).argmin()

    # let's do that regridding!
    data = np.empty((nt, nyout, nxout))
    for kt in np.arange(nt): 
        print('working on step' + str(kt))
        # get current frame
        data_step = data_ctrl.isel(time=kt).values.ravel()

        # interpolate
        data_interp = _pyresample.kd_tree.resample_nearest(grid_src, data_step,
                                                           grid_target,
                                                           radius_of_influence=50000,
                                                           fill_value=0)
        # roll back from -180/180 to 0/360
        data[kt,:,:] = np.roll(data_interp, -ind_greenwich)

    ds = xr.Dataset({varname: (['time', 'lat', 'lon'], data)},
                     coords={'lon': (['lon'], lon_atm),
                             'lat': (['lat'], lat_atm),
                             'time': pd.date_range(refdate, freq='14D', periods=419),
                             'reference_time': pd.Timestamp(refdate)})

    return ds

def regrid_to_model(ds, filegrid, varname):

    grid = xr.open_dataset(filegrid)

    mask_face2 = np.ones((270, 270))
    mask_face2[130:140,130:140] = 0
    mask_north_pole = {2: mask_face2}

    out = MITgcm_recipes.regridding.regrid_2_mitgcm_llc(ds, grid, [varname], 
                                                        method='bilinear',
                                                        faces2blend=[2], blend_mask=mask_north_pole,
                                                        lonname='lon', latname='lat', point='T', 
                                                        periodic=True, reuse_weights=True)
    out = out.transpose(*('time','face','j','i'))
    return out

def resample_ctrl_to_forcing_dates(ds_ctrl, year, freq='3H'):
    ''' resample the ctrl to the actual forcing date/time '''
    #cfirstyear = str(ds_ctrl['time'].dt.year.min().values)
    
    # we want to make sure we don't miss the first and last ctrl
    # for the year. They can be up to 13 days before Jan 01 and
    # 13 days after Dec 12, rounding up to middle of month

    startdate=str(year-1) + '-12-15'
    enddate=str(year+1) + '-01-15'

    ds_out=ds_ctrl.sel(time=slice(startdate, enddate)).resample(time=freq).interpolate('linear')

    return ds_out

def add_correction(ds_forcing, ds_ctrl_resampled, mult_fact):
    ''' add the resampled correction onto forcing field '''
    # xarray is going to broadcast ds_ctrl_resampled onto the forcing time
    ds_out = ds_forcing + mult_fact * ds_ctrl_resampled
    ds_out = ds_out.transpose(*('time','face','j','i'))
    return ds_out

def write_to_binary(ds, variable, fileout, precision='single'):
    ''' write variable from dataset ds to fileout with precision '''
    # write data to binary files
    fid   = open(fileout, "wb")
    print(ds[variable].values.shape)
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
