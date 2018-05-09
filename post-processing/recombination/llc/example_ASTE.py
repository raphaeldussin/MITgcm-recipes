from lib_regional_llc import *
from matplotlib import cm

dict_plt_sst = {'vmin':-2,'vmax':32,'colorbar':cm.gist_ncar}
dict_plt_sss = {'vmin':23,'vmax':40,'colorbar':cm.jet}

aste = regional_llc()
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results','T',dict_plt_sst,level=0,timestep=2307,precision='single')
aste.plot_map_field('/Users/raphael/WORK/ASTE/test_results','S',dict_plt_sss,level=0,timestep=2307,precision='single')
