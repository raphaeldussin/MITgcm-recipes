#!/usr/bin/env python

import xarray as xr
from lib_regrid2MITgcm import *

# artemis
dirgecco2='/local/data/artemis/simulations/GECCO2/'
dirout='/local/data/artemis/workspace/rdussin/ASTE/GECCO2_regridded_ASTE/'
dirastegrd='/local/data/artemis/workspace/rdussin/ASTE/GRID/'

gecco2_global = xr.open_mfdataset(dirgecco2  + '*.nc', decode_times=False)
variables2regrid = ['zeta','temp','salt','u','v']
variables2regrid = ['zeta','temp']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset(dirastegrd + 'ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(gecco2_global,variables2regrid,'lon','lat','Depth',
                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
interpolated1.to_netcdf(dirout + 'GECCO2_' + tag + '_ASTE_FACET1_' + year + '.nc')

##--------------- FACET 3 ---------------------
#aste_grid3 = xr.open_dataset(dirastegrd + 'ASTE_FACET3.nc')
#Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
#interpolated3 = Regrid2ASTE3.regrid(gecco2_global,variables2regrid,'lon','lat','Depth',
#                                    method='bilinear',blend_missing=True,periodicity=0,periodic=False)
#interpolated3.to_netcdf(dirout + 'GECCO2_' + tag + '_ASTE_FACET3_' + year + '.nc')
#
##--------------- FACET 4 ---------------------
#aste_grid4 = xr.open_dataset(dirastegrd + 'ASTE_FACET4.nc')
#Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
#interpolated4 = Regrid2ASTE4.regrid(gecco2_global,variables2regrid,'lon','lat','Depth',
#                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
#interpolated4.to_netcdf(dirout + 'GECCO2_' + tag + '_ASTE_FACET4_' + year + '.nc')
#
##--------------- FACET 5 ---------------------
#aste_grid5 = xr.open_dataset(dirastegrd + 'ASTE_FACET5.nc')
#Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
#interpolated5 = Regrid2ASTE5.regrid(gecco2_global,variables2regrid,'lon','lat','Depth',
#                                    method='bilinear',blend_missing=False,periodicity=0,periodic=True)
#interpolated5.to_netcdf(dirout + 'GECCO2_' + tag + '_ASTE_FACET5_' + year + '.nc')
#
