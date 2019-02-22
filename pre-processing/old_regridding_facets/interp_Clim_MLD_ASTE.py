#!/usr/bin/env python

import xarray as xr
from lib_regrid2MITgcm import *

# artemis
dirclim='/local/data/artemis/workspace/rdussin/Observations/MLD_deBoyerMontegut/originals/'
dirout='/local/data/artemis/workspace/rdussin/Observations/MLD_deBoyerMontegut/interp_ASTE/'
dirastegrd='/local/data/artemis/workspace/rdussin/ASTE/GRID/'

ds_MLD_in = xr.open_dataset(dirclim + 'mld_DR003_c1m_reg2.0.nc', 
                            decode_times=False)
ds_MLD_in['mld'] = ds_MLD_in['mld'].where(ds_MLD_in['mld'] < 1.e+8)

variables2regrid = ['mld']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirastegrd + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(ds_MLD_in, variables2regrid, 'lon', 'lat',
                                    '', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='time')
interpolated1.to_netcdf(dirout + 'MLD_Clim_DR003_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
# not enough data in Arctic, that'd be a waste a time

#--------------- FACET 4 ---------------------
# not enough data here too

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirastegrd + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(ds_MLD_in, variables2regrid, 'lon', 'lat',
                                    '', method='bilinear', blend_missing=False, drown=False,
                                    periodicity=0, periodic=True, timevarname='time')
interpolated5.to_netcdf(dirout + 'MLD_Clim_DR003_ASTE_FACET5.nc')

