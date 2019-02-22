import xarray as xr
from lib_regrid2MITgcm import *

# artemis
dirin    = '/local/data/artemis/workspace/rdussin/ASTE/Iron_dust/'
dirgrid  = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dirout   = '/local/data/artemis/workspace/rdussin/ASTE/Iron_dust/'

bling_global = xr.open_dataset(dirin + 'dust_deposition_global.nc')
variables2regrid = ['iron_dust']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirgrid + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(bling_global, variables2regrid, 'lon', 'lat', depthvarname='',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated1.to_netcdf(dirout + 'IronDust_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirgrid + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(bling_global, variables2regrid, 'lon', 'lat', depthvarname='',
                                    method='bilinear', blend_missing=True, periodicity=0, periodic=False)
interpolated3.to_netcdf(dirout + 'IronDust_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirgrid + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(bling_global, variables2regrid, 'lon', 'lat', depthvarname='',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated4.to_netcdf(dirout + 'IronDust_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirgrid + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(bling_global, variables2regrid, 'lon', 'lat', depthvarname='',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated5.to_netcdf(dirout + 'IronDust_ASTE_FACET5.nc')

