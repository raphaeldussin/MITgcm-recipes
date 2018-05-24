import xarray as xr
from lib_regrid2MITgcm import *

#bling_global = xr.open_dataset('/Volumes/L2/ASTE/global_BLINGv2/C1p5.control.3201-3300.ocean_bling_avg.nc')
bling_global = xr.open_dataset('C1p5.control.3201-3300.ocean_bling_avg.nc')
variables2regrid = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
#variables2regrid = ['alk']

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated1.to_netcdf('bling_ASTE_FACET1.nc')

#--------------- FACET 2 ---------------------
#aste_grid2 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET2.nc')
#Regrid2ASTE2 = Regrid2MITgcm(aste_grid2)
#interpolated2 = Regrid2ASTE2.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
#interpolated2.to_netcdf('bling_ASTE_FACET2.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0,periodic=False)
interpolated3.to_netcdf('bling_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated4.to_netcdf('bling_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated5.to_netcdf('bling_ASTE_FACET5.nc')

