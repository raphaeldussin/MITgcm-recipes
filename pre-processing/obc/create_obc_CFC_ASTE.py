from lib_create_obc import *
import xarray as xr

datasetdir='/Volumes/L2/ASTE/glodap1_regridded_ASTE/'
nameroot='GLODAP_CFC_ASTE_FACET<FACET>.nc'
obc_dirout='/Volumes/L2/ASTE/INPUTS/inputs_cfc_from_glodap+carmack/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['CFC11','CFC12']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/Users/raphael/WORK/MITgcm_runs/ASTE-TEST01/input_BE2_dthetadr/data.obcs')

for var in list_bgc_variables:
	aste_obc.create_obc_regional_llc(var,'south',list_dataset,obc_dirout + var + '_obc_south_monthly_glodap+carmack.bin',precision='single',repeats=289)
	aste_obc.create_obc_regional_llc(var,'east',list_dataset,obc_dirout + var + '_obc_east_monthly_glodap+carmack.bin',precision='single',repeats=289)
