import xmitgcm 

dirroot = '/local/data/artemis/workspace/rdussin/ASTE/GRID/'
dir_mitgrid = dirroot + 'grid_aste270/'
dir_nc = dirroot + 'nc/'

# load aste metadata
astemd = xmitgcm.utils.get_extra_metadata(domain='aste', nx=270)

# load grid files into xarray data structure
grid = xmitgcm.utils.get_grid_from_input(dir_mitgrid + 'tile<NFACET>.mitgrid',
                                         extra_metadata=astemd)

# write dataset to netcdf file
grid.to_netcdf(dir_nc + 'aste_grid.nc')
