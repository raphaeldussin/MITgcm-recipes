from lib_create_obc import *
import xarray as xr
import os

datasetdir='/local/data/artemis/workspace/rdussin/ASTE/global_BLINGv2_regridded_ASTE/'
nameroot='bling_ASTE_FACET<FACET>.nc'
obc_dirout='/local/data/artemis/workspace/rdussin/ASTE/INPUTS/inputs_bgc_from_global_BLINGv2/'

list_facets = ['1','3','4','5']
list_bgc_variables = ['alk','dic','o2','no3','po4','fed','n_org','p_org']
list_dataset = []

for facet in list_facets:
	list_dataset.append(xr.open_dataset(datasetdir + nameroot.replace('<FACET>',facet)))

aste_obc = create_obc('/home/rdussin/MITgcm_runs/ASTE-TEST01/input_BE2_dthetadr/data.obcs')

for var in list_bgc_variables:
	for bdry in ['south','east']:
		# december at the beginning
		fileout1 = obc_dirout + var + '_obc_' + bdry + '_monthly_global_bling_v2_rec0.bin'
		aste_obc.create_obc_regional_llc(var,bdry,list_dataset,fileout1,record=11)
	
		# jan to december for 24 years
		fileout2 = obc_dirout + var + '_obc_' + bdry + '_monthly_global_bling_v2_rec1-289.bin'
		aste_obc.create_obc_regional_llc(var,bdry,list_dataset,fileout2,repeats=23)
	
		# january at the end
		fileout3 = obc_dirout + var + '_obc_' + bdry + '_monthly_global_bling_v2_rec290.bin'
		aste_obc.create_obc_regional_llc(var,bdry,list_dataset,fileout3, record=0)

		# concat all the files
		fileout = obc_dirout + var + '_obc_' + bdry + '_monthly_global_bling_v2.bin'
		os.system('cat ' + fileout1 + ' ' + fileout2 + ' ' + fileout3 + ' > ' + fileout )
		os.system('rm ' + fileout1 + ' ' + fileout2 + ' ' + fileout3 )
