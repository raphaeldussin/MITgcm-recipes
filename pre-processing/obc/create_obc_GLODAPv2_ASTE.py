from lib_create_obc import *
import xarray as xr

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/GLODAPv2_regridded_ASTE/'
nameroot='GLODAPv2_ASTE_FACET<FACET>.nc'
obc_dirout='/local/data/artemis/workspace/rdussin/ASTE/INPUTS/inputs_bgc_from_GLODAPv2/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['TAlk','TCO2','oxygen','NO3','PO4']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/home/rdussin/MITgcm_runs/ASTE-Phy-REF01/input_BE2_dthetadr/data.obcs')

for var in list_bgc_variables:
	aste_obc.create_obc_regional_llc(var,'south',list_dataset,obc_dirout + \
                                         var + '_obc_south_annual_GLODAPv2.bin')
	aste_obc.create_obc_regional_llc(var,'east',list_dataset,obc_dirout + \
                                         var + '_obc_east_annual_GLODAPv2.bin')
