from lib_regional_llc import *
from matplotlib import cm

dict_plt_sst = {'vmin':-2,'vmax':32,'colorbar':cm.gist_ncar}
dict_plt_sss = {'vmin':23,'vmax':40,'colorbar':cm.jet}
dict_plt_ssh = {'vmin':-2,'vmax':1,'colorbar':cm.jet}
dict_plt_depth = {'vmin':0,'vmax':6000,'colorbar':cm.jet}
dict_plt_bathy = {'vmin':-6000,'vmax':0,'colorbar':cm.jet}

dir_results='/Users/raphael/WORK/ASTE/test_results'
dir_inputs='/Users/raphael/WORK/ASTE/input_files'

aste = regional_llc()
aste.plot_map_field(dir_results,'T',dict_plt_sst,level=0,timestep=2307,precision='single')
aste.plot_map_field(dir_results,'S',dict_plt_sss,level=0,timestep=2307,precision='single')
aste.plot_map_field(dir_results,'Eta',dict_plt_ssh,timestep=2307,precision='single')
aste.plot_map_field(dir_results,'Depth',dict_plt_depth,precision='single')

aste.plot_map_field(dir_inputs,'',dict_plt_sst,level=0,filealt='WOA09v2_T_llc270_JAN.bin',precision='double',\
dirxy=dir_results)

aste.plot_map_field(dir_results,'',dict_plt_bathy,filealt='bathy_fill9iU42Ef_noStLA_v1.bin',precision='double',\
dirxy=dir_results)

# for High Resolution ASTE
aste_hr = regional_llc(domain='ASTE_hr')
#...
