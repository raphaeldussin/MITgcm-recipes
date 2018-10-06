import xarray as xr
from lib_regrid2MITgcm import *

# artemis
dirwoa   = '/local/data/artemis/workspace/rdussin/ASTE/Sillicate_input/'
dirgrid  = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dirout   = '/local/data/artemis/workspace/rdussin/ASTE/Sillicate_input/'

woa_global = xr.open_mfdataset(dirwoa + 'woa13_all_i??_01.nc', decode_times=False)
variables2regrid = ['i_an']

def convert_units_clip(dataset, variables):
    for var in variables:
        # convert from WOA to MITgcm BLING units (umol/L -> mol/m3)
        dataset[var] = dataset[var] * 1000 / 1.0e+6 
        # remove negative values
        dataset[var] = dataset[var].clip(min=0)
    return dataset

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirgrid + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated1 = convert_units_clip(interpolated1, variables2regrid)
interpolated1.to_netcdf(dirout + 'WOA_SiO4_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirgrid + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=True, periodicity=0, periodic=False)
interpolated3 = convert_units_clip(interpolated3, variables2regrid)
interpolated3.to_netcdf(dirout + 'WOA_SiO4_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirgrid + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated4 = convert_units_clip(interpolated4, variables2regrid)
interpolated4.to_netcdf(dirout + 'WOA_SiO4_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirgrid + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(woa_global, variables2regrid, 'lon', 'lat', depthvarname='depth',
                                    method='bilinear', blend_missing=False, periodicity=0, periodic=True)
interpolated5 = convert_units_clip(interpolated5, variables2regrid)
interpolated5.to_netcdf(dirout + 'WOA_SiO4_ASTE_FACET5.nc')

