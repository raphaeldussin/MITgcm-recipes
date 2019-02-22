from lib_create_ic import *
import xarray as xr

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/GLODAPv2_regridded_ASTE/'
nameroot='GLODAPv2_ASTE_FACET<FACET>.nc'

dir_ic_out='/local/data/artemis/workspace/rdussin/ASTE/INPUTS/inputs_bgc_from_GLODAPv2/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['TAlk','TCO2','oxygen','NO3','PO4']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_ic = create_ic()

for var in list_bgc_variables:
	aste_ic.create_ic_regional_llc(var,list_dataset,dir_ic_out + var + '_IC_annual_mean_GLODAPv2_f64.bin',
                                       precision='double')
