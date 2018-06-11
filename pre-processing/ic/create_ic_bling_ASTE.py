from lib_create_ic import *
import xarray as xr

datasetdir='../regridding/'
nameroot='bling_ASTE_FACET<FACET>.nc'

list_facets = ['1','3','4','5']
#list_bgc_variables = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
list_bgc_variables = ['alk']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_ic = create_ic()

for var in list_bgc_variables:
	aste_ic.create_ic_regional_llc(var,list_dataset,var + '_IC__bling.bin',timestep=0)
