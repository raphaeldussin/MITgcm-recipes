#!/usr/bin/env python

import xarray as xr
from lib_regrid2MITgcm import *
import sys

# artemis
dirsoda='/local/data/artemis/simulations/SODA/3.3.1/monthly/'
dirout='/local/data/artemis/workspace/rdussin/ASTE/SODA3.3.1_regridded_ASTE/'
dirastegrd='/local/data/artemis/workspace/rdussin/ASTE/GRID/'

year=str(sys.argv[-1])

soda_global = xr.open_dataset(dirsoda  + 'soda3.3.1_mn_ocean_reg_' + year +'.nc')
variables2regrid = ['temp','salt','u','v']
#variables2regrid = ['ssh','taux','tauy','mlt','mls','mlp','anompb','net_heating','salt_flux_total']
tag = '3d' # '3d'

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirastegrd + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(soda_global,variables2regrid,'longitude','latitude','depth',
                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
interpolated1.to_netcdf(dirout + 'SODA3.3.1_' + tag + '_ASTE_FACET1_' + year + '.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset(dirastegrd + 'ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(soda_global,variables2regrid,'longitude','latitude','depth',
                                    method='bilinear',blend_missing=True,periodicity=0,periodic=False)
interpolated3.to_netcdf(dirout + 'SODA3.3.1_' + tag + '_ASTE_FACET3_' + year + '.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset(dirastegrd + 'ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(soda_global,variables2regrid,'longitude','latitude','depth',
                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
interpolated4.to_netcdf(dirout + 'SODA3.3.1_' + tag + '_ASTE_FACET4_' + year + '.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset(dirastegrd + 'ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(soda_global,variables2regrid,'longitude','latitude','depth',
                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
interpolated5.to_netcdf(dirout + 'SODA3.3.1_' + tag + '_ASTE_FACET5_' + year + '.nc')

