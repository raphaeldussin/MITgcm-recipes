#!/usr/bin/env python

import xarray as xr
from lib_regrid2MITgcm import *
import sys

# artemis
dirclim='/local/data/artemis/workspace/rdussin/Observations/RG_ARGO_clim/originals/'
dirout='/local/data/artemis/workspace/rdussin/Observations/RG_ARGO_clim/interp_ASTE/'
dirastegrd='/local/data/artemis/workspace/rdussin/ASTE/GRID/'

year=str(sys.argv[-1])

ds_TS_in = xr.open_mfdataset([dirclim + 'RG_ArgoClim_Temperature_2017.nc', 
                             dirclim + 'RG_ArgoClim_Salinity_2017.nc'], decode_times=False)
variables2regrid_mean = ['ARGO_TEMPERATURE_MEAN', 'ARGO_SALINITY_MEAN']
variables2regrid_anom = ['ARGO_TEMPERATURE_ANOMALY', 'ARGO_SALINITY_ANOMALY']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirastegrd + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(ds_TS_in, variables2regrid_mean, 'LONGITUDE', 'LATITUDE',
                                    'PRESSURE', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='TIME')
interpolated1.to_netcdf(dirout + 'RG_Argo_Clim_ASTE_FACET1.nc')

interpolated1a = Regrid2ASTE1.regrid(ds_TS_in, variables2regrid_anom, 'LONGITUDE', 'LATITUDE',
                                    'PRESSURE', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='TIME')
interpolated1a.to_netcdf(dirout + 'RG_Argo_Anom_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
# not enough data in Arctic, that'd be a waste a time

#--------------- FACET 4 ---------------------
# not enough data here too

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirastegrd + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(ds_TS_in, variables2regrid_mean, 'LONGITUDE', 'LATITUDE',
                                    'PRESSURE', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='TIME')
interpolated5.to_netcdf(dirout + 'RG_Argo_Clim_ASTE_FACET5.nc')

interpolated5a = Regrid2ASTE5.regrid(ds_TS_in, variables2regrid_anom, 'LONGITUDE', 'LATITUDE',
                                    'PRESSURE', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='TIME')
interpolated5a.to_netcdf(dirout + 'RG_Argo_Anom_ASTE_FACET5.nc')
