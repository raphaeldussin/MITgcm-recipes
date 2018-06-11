from lib_create_obc import *
import xarray as xr

datasetdir='../regridding/'
nameroot='bling_ASTE_FACET<FACET>.nc'

list_facets = ['1','3','4','5']
list_bgc_variables = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/Users/raphael/WORK/MITgcm_runs/ASTE-TEST01/input_BE2_dthetadr/data.obcs')

for var in list_bgc_variables:
	aste_obc.create_obc_regional_llc(var,'south',list_dataset,var + '_south_bling.bin')
	aste_obc.create_obc_regional_llc(var,'east',list_dataset,var + '_east_bling.bin')
