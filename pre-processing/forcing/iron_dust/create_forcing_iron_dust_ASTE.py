from lib_create_2d_forcing import *
import xarray as xr
import os
import numpy as np

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/Iron_dust/'
nameroot='IronDust_ASTE_FACET<FACET>.nc'

dir_ic_out='/local/data/artemis/workspace/rdussin/ASTE/Iron_dust/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['iron_dust']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_iron = create_ic()

for var in list_bgc_variables:
	for kt in np.arange(12):
		aste_iron.create_ic_regional_llc(var,list_dataset,'iron_month' + str(kt+1).zfill(2),precision='double',timestep=kt)

# concat all the files
fileout = dir_ic_out + 'iron_deposition.bin'
os.system('cat iron_month?? > ' + fileout)
os.system('rm iron_month??')
