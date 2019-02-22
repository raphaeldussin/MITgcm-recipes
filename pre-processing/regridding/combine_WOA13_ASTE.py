import MITgcm_recipes
import xarray as xr

#first open one model dataset and pick a dataarray with the dims we want
model_dir='/local/data/artemis/workspace/rdussin/ASTE/POST/climato/'
model_clim=model_dir + 'monthly_clim_ASTE_Run1.nc'
ds = xr.open_dataset(model_clim)
da = ds['THETA']

#combine the facets into one xarray dataset
woa_dir='/local/data/artemis/workspace/rdussin/Observations/WOA13/interp_ASTE/'
woa_no3 = MITgcm_recipes.facets2faces_aste270(da, 'n_an', 
                                              facet1=woa_dir + 'WOA_NO3_ASTE_FACET1.nc',
                                              facet3=woa_dir + 'WOA_NO3_ASTE_FACET3.nc',
                                              facet4=woa_dir + 'WOA_NO3_ASTE_FACET4.nc',
                                              facet5=woa_dir + 'WOA_NO3_ASTE_FACET5.nc')


woa13 = xr.Dataset({'no3': (['time', 'k', 'face', 'j', 'i'],  woa_no3) 
                   },
                   coords = da.coords)

woa13.to_netcdf('test.nc')
