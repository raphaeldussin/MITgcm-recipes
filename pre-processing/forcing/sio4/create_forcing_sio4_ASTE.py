from lib_create_2d_forcing2 import *
import xarray as xr
import os
import numpy as np

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/Sillicate_input/'
nameroot='WOA_SiO4_ASTE_FACET<FACET>.nc'

dir_ic_out='/local/data/artemis/workspace/rdussin/ASTE/Sillicate_input/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['i_an']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_iron = create_ic()

for var in list_bgc_variables:
	for kt in np.arange(12):
		aste_iron.create_ic_regional_llc(var,list_dataset,'SiO4_month' + str(kt+1).zfill(2),precision='single',timestep=kt,level=0)

# concat all the files
fileout = dir_ic_out + 'sillev_raf.bin'
os.system('cat SiO4_month?? > ' + fileout)
os.system('rm SiO4_month??')
