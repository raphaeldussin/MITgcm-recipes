import xarray as xr
from lib_regrid2MITgcm import *

# laptop
dirbling = '/Volumes/L2/ASTE/global_BLINGv2/'
dirgrid  = '/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/'
dirout   = '/Volumes/L2/ASTE/global_BLINGv2_regridded_ASTE/'

# artemis
dirbling = '/local/data/artemis/workspace/rdussin/ASTE/global_BLINGv2/'
dirgrid  = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dirout   = '/local/data/artemis/workspace/rdussin/ASTE/global_BLINGv2_regridded_ASTE/'

bling_global = xr.open_dataset(dirbling + 'C1p5.control.3201-3300.ocean_bling_avg.nc')
variables2regrid = ['alk','dic','o2','no3','po4','fed','n_org','p_org']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirgrid + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(bling_global, variables2regrid, 'geolon_t', 'geolat_t',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated1.to_netcdf(dirout + 'bling_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirgrid + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(bling_global, variables2regrid, 'geolon_t', 'geolat_t',
                                    method='bilinear', blend_missing=True, periodicity=0, periodic=False)
interpolated3.to_netcdf(dirout + 'bling_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirgrid + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(bling_global, variables2regrid, 'geolon_t', 'geolat_t',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated4.to_netcdf(dirout + 'bling_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirgrid + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(bling_global, variables2regrid, 'geolon_t', 'geolat_t',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated5.to_netcdf(dirout + 'bling_ASTE_FACET5.nc')

