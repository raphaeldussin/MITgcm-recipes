from lib_create_obc import *
import xarray as xr

datasetdir='../regridding/'
nameroot='bling_ASTE_FACET<FACET>.nc'

list_dataset = []
for facet in ['1','3','4','5']:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/Users/raphael/WORK/MITgcm_runs/ASTE-TEST01/input_BE2_dthetadr/data.obcs')
aste_obc.create_obc_regional_llc('alk','south',list_dataset,'alk_south_bling.bin')
aste_obc.create_obc_regional_llc('alk','east',list_dataset,'alk_east_bling.bin')
