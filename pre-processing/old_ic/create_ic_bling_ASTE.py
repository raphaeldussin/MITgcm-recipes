from lib_create_ic import *
import xarray as xr

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/global_BLINGv2_regridded_ASTE/'
nameroot='bling_ASTE_FACET<FACET>.nc'

dir_ic_out='/local/data/artemis/workspace/rdussin/ASTE/INPUTS/inputs_bgc_from_global_BLINGv2/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_ic = create_ic()

for var in list_bgc_variables:
	aste_ic.create_ic_regional_llc(var,list_dataset,dir_ic_out + var + '_IC_January_global_blingv2_f64.bin', \
                                       timestep=0, precision='double')
