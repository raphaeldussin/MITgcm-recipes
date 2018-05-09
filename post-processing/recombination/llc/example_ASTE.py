from lib_regional_llc import *
from matplotlib import cm

dict_plt_sst = {'vmin':-2,'vmax':32,'colorbar':cm.gist_ncar}
dict_plt_sss = {'vmin':23,'vmax':40,'colorbar':cm.jet}
dict_plt_ssh = {'vmin':-2,'vmax':1,'colorbar':cm.jet}
dict_plt_bathy = {'vmin':0,'vmax':5000,'colorbar':cm.jet}

aste = regional_llc()
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results','T',dict_plt_sst,level=0,timestep=2307,precision='single')
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results','S',dict_plt_sss,level=0,timestep=2307,precision='single')
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results','Eta',dict_plt_ssh,timestep=2307,precision='single')

# bug plot
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results/','',dict_plt_sst,level=0,filealt='../input_files/WOA09v2_T_llc270_JAN.bin',precision='double')

# on llc global grid
#aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results/','',dict_plt_bathy,filealt='../input_files/bathy_fill9iU42Ef_noStLA_v1.bin',precision='single')

