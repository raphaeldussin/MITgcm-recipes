import xarray as xr
from lib_regrid2MITgcm import *

# artemis
dirglodap = '/local/data/artemis/workspace/rdussin/Observations/GLODAPv2/merged_for_regridding/'
dirgrid  = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dirout   = '/local/data/artemis/workspace/rdussin/ASTE/GLODAPv2_regridded_ASTE/'

# read dataset
glodap = xr.open_dataset(dirglodap + 'GLODAPv2_4_regridding.nc')
# add vertical coordinate
glodap['depth_surface'] = glodap['Depth']

variables2regrid = ['TAlk','TCO2','oxygen','NO3','PO4']

def convert_units_clip(dataset, variables):
    for var in variables:
        # convert from GLOBAL to BLING units (umol/kg -> mol/m3)
        dataset[var] = dataset[var] * 1035. / 1.0e+6
        # remove negative values
        dataset[var] = dataset[var].clip(min=0)
    return dataset

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirgrid + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(glodap, variables2regrid, 'lon', 'lat', depthvarname='depth_surface', reuse_weights=False,
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)

interpolated1 = convert_units_clip(interpolated1, variables2regrid)
interpolated1.to_netcdf(dirout + 'GLODAPv2_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirgrid + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(glodap, variables2regrid, 'lon', 'lat', depthvarname='depth_surface', reuse_weights=False,
                                    method='bilinear', blend_missing=True, periodicity=0, periodic=False)

interpolated3 = convert_units_clip(interpolated3, variables2regrid)
interpolated3.to_netcdf(dirout + 'GLODAPv2_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirgrid + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(glodap, variables2regrid, 'lon', 'lat', depthvarname='depth_surface', reuse_weights=False,
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)

interpolated4 = convert_units_clip(interpolated4, variables2regrid)
interpolated4.to_netcdf(dirout + 'GLODAPv2_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirgrid + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(glodap, variables2regrid, 'lon', 'lat', depthvarname='depth_surface', reuse_weights=False,
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)

interpolated5 = convert_units_clip(interpolated5, variables2regrid)
interpolated5.to_netcdf(dirout + 'GLODAPv2_ASTE_FACET5.nc')

