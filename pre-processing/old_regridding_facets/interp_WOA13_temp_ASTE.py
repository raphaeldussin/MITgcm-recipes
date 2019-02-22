import xarray as xr
from lib_regrid2MITgcm import *
import seawater
import numpy as np

# artemis
dirwoa   = '/local/data/artemis/workspace/rdussin/Observations/WOA13/'
dirgrid  = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dirout   = '/local/data/artemis/workspace/rdussin/Observations/WOA13/interp_ASTE/'

woa_global = xr.open_mfdataset(dirwoa + 'woa13_decav_[t,s]??_01.nc', decode_times=False)
# compute pot temperature
pres = np.zeros(woa_global['s_an'].values.shape)
for k in np.arange(woa_global['depth'].shape[0]):
    pres[:,k,:,:] = woa_global['depth'].values[k]

woa_global['theta'] = xr.DataArray(seawater.ptmp(woa_global['s_an'], woa_global['t_an'], pres),
coords=woa_global['s_an'].coords, dims=woa_global['s_an'].dims)
variables2regrid = ['theta']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirgrid + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated1.to_netcdf(dirout + 'WOA_THETA_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirgrid + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=False)
interpolated3.to_netcdf(dirout + 'WOA_THETA_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirgrid + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated4.to_netcdf(dirout + 'WOA_THETA_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirgrid + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated5.to_netcdf(dirout + 'WOA_THETA_ASTE_FACET5.nc')

