import xarray as xr
from lib_regrid2MITgcm import *
import scipy.ndimage.filters as smoothing
import numpy as np

dir_glodap = '/Volumes/L2/Observations/GLODAP/processed/'
variables2regrid = ['CFC11','CFC12']
dirout='/Volumes/L2/ASTE/glodap1_regridded_ASTE/'
dir_cfc_carmack = '/Volumes/L2/Observations/CFCs/results_cruise_carmack/'
dir_aste_grid = '/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/'

glodap_cfc = xr.open_mfdataset([dir_glodap+'CFC11_CF.nc',dir_glodap+'CFC12_CF.nc'])
unit_conversion = 1035 * 1.e-12 # conversion from pmol/kg to mol/m3
glodap_cfc['CFC11'] = glodap_cfc['CFC11'] * unit_conversion
glodap_cfc['CFC12'] = glodap_cfc['CFC12'] * unit_conversion

#--------------- FACET 1 ---------------------
aste_grid1 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET1.nc')
Regrid2ASTE1 = Regrid2MITgcm(aste_grid1)
interpolated1 = Regrid2ASTE1.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',\
blend_missing=True,periodicity=0,periodic=True)
# nordic seas need a bit of cleaning up spurious values from extrapolation
interpolated1['CFC11'][:,:,430:,160:] = interpolated1['CFC11'][:,:,430:,160] * interpolated1['lsm'][:,430:,160:]
interpolated1['CFC12'][:,:,430:,160:] = interpolated1['CFC12'][:,:,430:,160] * interpolated1['lsm'][:,430:,160:]
interpolated1['CFC11'] = xr.where(interpolated1['CFC11'] < 0,0,interpolated1['CFC11'])
interpolated1['CFC12'] = xr.where(interpolated1['CFC12'] < 0,0,interpolated1['CFC12'])
interpolated1.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET1.nc')

#--------------- FACET 4 ---------------------
aste_grid4 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET4.nc')
Regrid2ASTE4 = Regrid2MITgcm(aste_grid4)
interpolated4 = Regrid2ASTE4.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',\
blend_missing=True,periodicity=0,periodic=True)
# some small spurious gradient that I want to smooth out
interpolated4['CFC11'][:,:,:110,:20] = (interpolated4['CFC11'][:,:,110,:20] * \
                                        interpolated4['lsm'][:,:110,:20] ).transpose('time','z','y','x')
interpolated4['CFC12'][:,:,:110,:20] = (interpolated4['CFC12'][:,:,110,:20] * \
                                        interpolated4['lsm'][:,:110,:20] ).transpose('time','z','y','x')
interpolated4['CFC11'] = xr.where(interpolated4['CFC11'] < 0,0,interpolated4['CFC11'])
interpolated4['CFC12'] = xr.where(interpolated4['CFC12'] < 0,0,interpolated4['CFC12'])
interpolated4.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET4.nc')

#--------------- FACET 5 ---------------------
aste_grid5 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET5.nc')
Regrid2ASTE5 = Regrid2MITgcm(aste_grid5)
interpolated5 = Regrid2ASTE5.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='bilinear',\
blend_missing=True,periodicity=0,periodic=True)
# more clean up
interpolated5['CFC11'][:,:,:105,:15] = interpolated5['CFC11'][:,:,111:116,1:5].mean(dim=('x','y')) * \
interpolated5['lsm'][:,:105,:15]
interpolated5['CFC12'][:,:,:105,:15] = interpolated5['CFC12'][:,:,111:116,1:5].mean(dim=('x','y')) * \
interpolated5['lsm'][:,:105,:15]
interpolated5['CFC11'] = xr.where(interpolated5['CFC11'] < 0,0,interpolated5['CFC11'])
interpolated5['CFC12'] = xr.where(interpolated5['CFC12'] < 0,0,interpolated5['CFC12'])
interpolated5.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET5.nc')

#--------------- FACET 3 ---------------------
aste_grid3 = xr.open_dataset('/Users/raphael/TOOLS/MITgcm-recipes/pre-processing/grid/ASTE_FACET3.nc')
Regrid2ASTE3 = Regrid2MITgcm(aste_grid3)
interpolated3 = Regrid2ASTE3.regrid(glodap_cfc,variables2regrid,'longitude','latitude',method='nearest_s2d',\
blend_missing=False,periodicity=0,periodic=True)
# For the Arctic, there's basically everything to fix
# We use CFC average profiles from a cruise for the Arctic
cfc_profiles_Carmack = xr.open_dataset(dir_cfc_carmack + 'CFC_profiles_Arctic_from_Carmack_depth_ASTE.nc')
interpolated3['CFC11'][:] = cfc_profiles_Carmack['CFC11_arctic_profile'][:]
interpolated3['CFC12'][:] = cfc_profiles_Carmack['CFC12_arctic_profile'][:]
# We pad the sides with results from the other facets
# -> FACET 1 (NorthEast Atlantic) becomes western boundary condition for FACET 3
interpolated1 = xr.open_dataset(dirout + 'GLODAP_CFC_ASTE_FACET1.nc')
bdry_west = interpolated1.isel(y=-1)#.rename({'x':'y','y':'x0'},inplace=True)
interpolated3['CFC11'][:,:,:,0] = bdry_west['CFC11'].values[:,:,::-1]
interpolated3['CFC12'][:,:,:,0] = bdry_west['CFC12'].values[:,:,::-1]
# -> FACET 4 (Pacific) becomes eastern boundary condition for FACET 3
interpolated4 = xr.open_dataset(dirout + 'GLODAP_CFC_ASTE_FACET4.nc')
bdry_east = interpolated4.isel(x=0)
interpolated3['CFC11'][:,:,:,-1] = bdry_east['CFC11'].values[:,:,:]
interpolated3['CFC12'][:,:,:,-1] = bdry_east['CFC12'].values[:,:,:]
# -> FACET 5 (Northwest Atlantic/Canadian basin) is the northern boundary condition for FACET 3
interpolated5 = xr.open_dataset(dirout + 'GLODAP_CFC_ASTE_FACET5.nc')
bdry_north = interpolated5.isel(x=0)
interpolated3['CFC11'][:,:,-1,:] = bdry_north['CFC11'].values[:,:,::-1]
interpolated3['CFC12'][:,:,-1,:] = bdry_north['CFC12'].values[:,:,::-1]

interpolated3['CFC11'][:] = interpolated3['CFC11'][:] * interpolated3['lsm'][:]
interpolated3['CFC12'][:] = interpolated3['CFC12'][:] * interpolated3['lsm'][:]

# -> FACET 1 transition
padsize=30
for ji in np.arange(1,padsize):
    #alpha = interpolated3['lsm'].values[:,:,0] * float(padsize-ji) / padsize # no need for linear transition (smoothing)
    alpha = interpolated3['lsm'].values[:,:,0] ; beta = 1- alpha
    interpolated3['CFC11'][0,:,:,ji] = alpha * bdry_west['CFC11'].values[0,:,::-1] + \
                                       beta * interpolated3['CFC11'].values[0,:,:,ji]
    interpolated3['CFC12'][0,:,:,ji] = alpha * bdry_west['CFC12'].values[0,:,::-1] + \
                                       beta * interpolated3['CFC12'].values[0,:,:,ji]

# -> FACET 4 transition
padsize=30 ; nx=270
for ji in np.arange(nx-padsize,nx):
    #alpha = interpolated3['lsm'].values[:,:,-1] * float(ji-(nx-padsize)) / padsize
    alpha = interpolated3['lsm'].values[:,:,-1] ; beta = 1 - alpha
    interpolated3['CFC11'][0,:,:,ji] = alpha * bdry_east['CFC11'].values[0,:,:] + \
                                       beta * interpolated3['CFC11'].values[0,:,:,ji]
    interpolated3['CFC12'][0,:,:,ji] = alpha * bdry_east['CFC12'].values[0,:,:] + \
                                       beta * interpolated3['CFC12'].values[0,:,:,ji]

# -> FACET 5 transition
padsize=30 ; ny=270
for jj in np.arange(ny-padsize,ny):
    #alpha = interpolated3['lsm'].values[:,-1,:] * float(jj-(ny-padsize)) / padsize
    alpha = interpolated3['lsm'].values[:,-1,:] ; beta = 1 -alpha
    interpolated3['CFC11'][0,:,jj,:] = alpha * bdry_north['CFC11'].values[0,:,::-1] + \
                                       beta * interpolated3['CFC11'].values[0,:,jj,:]
    interpolated3['CFC12'][0,:,jj,:] = alpha * bdry_north['CFC12'].values[0,:,::-1] + \
                                       beta * interpolated3['CFC12'].values[0,:,jj,:]

# mask the partial result
interpolated3['CFC11'][:] = interpolated3['CFC11'][:] * interpolated3['lsm'][:]
interpolated3['CFC12'][:] = interpolated3['CFC12'][:] * interpolated3['lsm'][:]

# Smooth the fields
npasses=1
for n in np.arange(npasses):
	for k in np.arange(50):
		# CFC-11
		tmp = interpolated3['CFC11'].isel(z=k).values.squeeze()
		tmp[np.where(tmp ==0)] = tmp[130,130] # set land value to interior value
		smoothed = smoothing.gaussian_filter(tmp,4)
		interpolated3['CFC11'][:,k,:,:] = smoothed * interpolated3['lsm'][k,:,:]
		# CFC-12
		tmp = interpolated3['CFC12'].isel(z=k).values.squeeze()
		tmp[np.where(tmp ==0)] = tmp[130,130]
		smoothed = smoothing.gaussian_filter(tmp,4)
		interpolated3['CFC12'][:,k,:,:] = smoothed * interpolated3['lsm'][k,:,:]

interpolated3['CFC11'] = xr.where(interpolated3['CFC11'] < 0,0,interpolated3['CFC11'])
interpolated3['CFC12'] = xr.where(interpolated3['CFC12'] < 0,0,interpolated3['CFC12'])
# Always mask the end result
interpolated3['CFC11'][:] = interpolated3['CFC11'][:] * interpolated3['lsm'][:]
interpolated3['CFC12'][:] = interpolated3['CFC12'][:] * interpolated3['lsm'][:]

interpolated3.to_netcdf(dirout + 'GLODAP_CFC_ASTE_FACET3.nc')
