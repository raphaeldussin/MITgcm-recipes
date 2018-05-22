import xarray as xr
from lib_regrid2MITgcm import *

aste_grid = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET1.nc')
Regrid2ASTE = Regrid2MITgcm(aste_grid)

bling_global = xr.open_dataset('/Volumes/L2/ASTE/global_BLINGv2/C1p5.control.3201-3300.ocean_bling_avg.nc')
variables2regrid = ['alk','dic']

interpolated = Regrid2ASTE.regrid(bling_global,variables2regrid,'geolon_t','geolat_t',method='bilinear',blend_missing=True,periodicity=0)

interpolated.to_netcdf('bling_ASTE_FACET1.nc')
