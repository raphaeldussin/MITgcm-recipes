from lib_regional_llc import *
from matplotlib import cm

dict_plt_alk = {'vmin':0.001,'vmax':0.003,'colorbar':cm.jet}

dir_xy='/Users/raphael/WORK/ASTE/test_results'
dir_inputs='/Volumes/L2/ASTE/INPUTS/inputs_bgc_from_global_BLINGv2/'

aste = regional_llc()

aste.plot_map_field(dir_inputs,'',dict_plt_alk,level=20,filealt='alk_IC_January_global_blingv2.bin',precision='single',\
dirxy=dir_xy)

