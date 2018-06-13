import xarray as xr
from lib_regrid2MITgcm import *

dir_glodap = '/Volumes/L2/Observations/GLODAP/processed/'
variables2regrid = ['CFC11','CFC12']
dirout='/Volumes/L2/ASTE/glodap1_regridded_ASTE/'

glodap_cfc = xr.open_mfdataset([dir_glodap+'CFC11_CF.nc',dir_glodap+'CFC12_CF.nc'])
unit_conversion = 1035 * 1.e-12 # conversion from pmol/kg to mol/m3
glodap_cfc['CFC11'] = glodap_cfc['CFC11'] * unit_conversion
glodap_cfc['CFC12'] = glodap_cfc['CFC12'] * unit_conversion

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated1.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET1.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='nearest_s2d',blend_missing=False,periodicity=0,periodic=True)
interpolated3.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET3.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated4.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',blend_missing=True,periodicity=0,periodic=True)
interpolated5.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET5.nc')

