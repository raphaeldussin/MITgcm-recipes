from matplotlib import cm
import xmitgcm
import MITgcm_recipes
import numpy as np

# copy available_diagnostics.log from run CFC into data_dir
# ln -s icfile TRAC01.0000000000.data

astemd = xmitgcm.utils.get_extra_metadata(domain='aste', nx=270)

cfc11 = xmitgcm.open_mdsdataset(data_dir='/local/data/artemis/workspace/rdussin/ASTE/INPUTS/inputs_cfc_from_glodap+carmack/', 
                                grid_dir='/local/data/artemis/workspace/rdussin/ASTE/RUNS/ASTE-BLING-Run01/outputs/',
                                prefix=['TRAC01'], iters=[0],geometry='llc', extra_metadata=astemd,
                                nx=270, default_dtype=np.dtype('d'))


dict_plt = {'vmin': 0,'vmax': 6e-9, 'cmap': cm.jet, 'figsize': [7,7], 'title':'Initial Condition CFC11 (mol/m3)'}
fig = MITgcm_recipes.plot_ASTE(cfc11['TRAC01'].sel(k=0).isel(time=0), dict_plt)

fig.savefig('IC_CFC11.png')

