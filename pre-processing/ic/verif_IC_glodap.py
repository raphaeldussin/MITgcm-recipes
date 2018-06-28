from lib_regional_llc import *
from matplotlib import cm

dict_plt_cfc = {'mult_fact':1e+9,'vmin':0,'vmax':6,'colorbar':cm.jet}

dir_xy='/Users/raphael/WORK/ASTE/test_results'
dir_inputs='/Volumes/L2/ASTE/INPUTS/inputs_cfc_from_glodap+carmack/'

aste = regional_llc()

aste.plot_map_field(dir_inputs,'',dict_plt_cfc,level=0,filealt='CFC11_IC_glodap+carmack.bin',precision='single',\
dirxy=dir_xy)

aste.plot_map_field(dir_inputs,'',dict_plt_cfc,level=0,filealt='CFC12_IC_glodap+carmack.bin',precision='single',\
dirxy=dir_xy)
