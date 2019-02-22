from lib_create_ic import *
import xarray as xr

datasetdir='/Volumes/L2/ASTE/glodap1_regridded_ASTE/'
nameroot='GLODAP_CFC_ASTE_FACET<FACET>.nc'

dir_ic_out='/Volumes/L2/ASTE/INPUTS/inputs_cfc_from_glodap+carmack/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['CFC11','CFC12']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_ic = create_ic()

for var in list_bgc_variables:
	aste_ic.create_ic_regional_llc(var,list_dataset,dir_ic_out + var + '_IC_glodap+carmack+bckgrd.bin',precision='double',timestep=0,background=1e-15)
