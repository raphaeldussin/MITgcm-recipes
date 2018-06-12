from lib_create_obc import *
import xarray as xr

datasetdir='/Volumes/L2/ASTE/global_BLINGv2_regridded_ASTE/'
nameroot='bling_ASTE_FACET<FACET>.nc'
obc_dirout='/Volumes/L2/ASTE/INPUTS/inputs_bgc_from_global_BLINGv2/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/Users/raphael/WORK/MITgcm_runs/ASTE-TEST01/input_BE2_dthetadr/data.obcs')

for var in list_bgc_variables:
	aste_obc.create_obc_regional_llc(var,'south',list_dataset,obc_dirout + var + '_obc_south_monthly_global_bling_v2.bin')
	aste_obc.create_obc_regional_llc(var,'east',list_dataset,obc_dirout + var + '_obc_east_monthly_global_bling_v2.bin')
