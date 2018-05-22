import xesmf as xe
import xarray as xr
import mod_drown_sosie
import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as spint
import mod_akima_1d

class Regrid2MITgcm():

	def __init__(self,model_grid):
		''' init the regridding object with model grid
		type(model_grid) = xarray '''
		self.model_grid = model_grid
		self.spval = -9999.
		return None

	def regrid(self,dataset,list_variables,lonvarname,latvarname,method='bilinear',blend_missing=True,periodicity=0,target_point='center',mask_output=True):
		''' regrid the input xarray dataset
		list_vars = list of strings
		'''
		# xESMF wants lon/lat to be the names of coords for interpolation
		inputds = dataset.rename({lonvarname: 'lon', latvarname: 'lat'}, inplace=True)
		# define an output dataset from grid with correct destination points
		if target_point == 'center':
			target_grid = self.model_grid.rename({'XC': 'lon', 'YC': 'lat'}, inplace=True)
		else:
			pass # TO DO

		for variable in list_variables:
			# build land sea mask
			maskvar = self._mask_from_missing_value(inputds[variable])
			# extrapolate onto land for coastal values
			inputds[variable] = self._drown_field(inputds[variable],maskvar,periodicity=periodicity)

			#plt.figure()
			#inputds[variable][0,0,:,:].plot() ; plt.show()

		hremapped = self._perform_interpolation(inputds,target_grid,list_variables,method=method,blend_missing=blend_missing)

		# vertical interpolation
		outputds = self. _vertical_interpolation(inputds,hremapped,target_grid,list_variables)

		# TO DO

		if mask_output:
			outputds = self._mask_output(outputds,target_grid,list_variables)

		return outputds


	def _mask_from_missing_value(self,dataarray):
		''' create binary 1/0 mask for masked values in dataarray '''
		mask = np.ones(dataarray.shape)
		mask[np.isnan(dataarray.values)] = 0
		return mask

	def _drown_field(self,dataarray,mask,periodicity=0):
		''' periodicity : 0 = periodic, 2 = 2 points overlap,...
		periodicity : -1 =non-periodic

		'''
		# we assume data array contains 2 horizontal dimensions
		# previous needs to be iterated on (t,z,) or (z,) or (t,)
		ndims_loop = len(dataarray.dims) - 2

		if ndims_loop == 2:
			for k0 in np.arange(len(dataarray.coords[dataarray.dims[0]])):
				for k1 in np.arange(len(dataarray.coords[dataarray.dims[1]])):

					# extract the 2d fields for drown
					tmpin = dataarray[k0,k1,:,:].values.squeeze()
					masklevel = mask[k0,k1,::].squeeze()
					# the NaNs in xarray don't make fortran happy
					tmpin[np.isnan(tmpin)] = self.spval
					tmpout = mod_drown_sosie.mod_drown.drown(periodicity,tmpin.transpose(),masklevel.transpose(),\
					                                         nb_inc=100,nb_smooth=20)
					dataarray[k0,k1,:,:] = tmpout.transpose()
		else:
			pass # TO DO


		return dataarray

	def _blend(self,da1,da2,missing):
		tmp1 = da1.values
		tmp2 = da2.values
		mask = np.ones(tmp1.shape)
		mask[np.where(tmp1 == missing)] = 0
		tmp3 = (mask * tmp1) + ((1-mask) * tmp2)
		da3 = xr.DataArray(tmp3,dims=da1.dims)
		return da3

	def _perform_interpolation(self,input_dataset,target_grid,list_variables,method='bilinear',blend_missing=True):

		# create interpolators
		# not sure if I want to re-use the weights (source of potential mistakes)
		regridder = xe.Regridder(input_dataset, target_grid, method) #,reuse_weights=True)
		if blend_missing:
			backup_regridder = xe.Regridder(input_dataset, target_grid, 'nearest_s2d') #,reuse_weights=True)

		#hremapped = target_grid.copy()
		hremapped = xr.Dataset()
		hremapped.update({'XC':target_grid['lon']})
		hremapped.update({'YC':target_grid['lat']})

		for variable in list_variables:
			out_da = regridder(input_dataset[variable])
			if blend_missing:
				backup_da = backup_regridder(input_dataset[variable])
				out_da = self._blend(out_da,backup_da,0)

			hremapped.update({variable:out_da})

		print(hremapped)
		return hremapped


	def _vertical_interpolation(self,inputds,hremapped,target_grid,list_variables):

		outputds = target_grid.copy()
		for variable in list_variables:
			# needs failsafe for 2d data
			ndims_loop = len(hremapped[variable].dims) - 3
			# needs failsafe for 2d data
			if len(hremapped[variable].dims) == 4:
				vert_axis = 1
				nt = len(hremapped[variable].coords[hremapped[variable].dims[0]])
				#nz = len(hremapped[variable].coords[hremapped[variable].dims[1]])
				nz = len(target_grid['Z'].values) # hardcoded
				ny = len(hremapped[variable].coords[hremapped[variable].dims[2]])
				nx = len(hremapped[variable].coords[hremapped[variable].dims[3]])

				dest_dims = (hremapped[variable].dims[0],'z',hremapped[variable].dims[2],hremapped[variable].dims[3])
			#elif len(hremapped[variable].dims) == 3:
			#	vert_axis = 0
			#	nt = 1
			#	nz = len(hremapped[variable].coords[hremapped[variable].dims[0]])
			#	ny = len(hremapped[variable].coords[hremapped[variable].dims[1]])
			#	nx = len(hremapped[variable].coords[hremapped[variable].dims[2]])

			# get data for src and dest depth
			name_zold = hremapped[variable].dims[vert_axis]

			#print(name_zold)
			zold = inputds[name_zold].values
			#print(zold)
			znew = target_grid['Z'].values # hardcoded

			#print(outputds[variable].dims)
			#print(ndims_loop,vert_axis)
			#print(znew)
			#data_regridded = xr.DataArray(

			remapped = xr.DataArray(np.zeros((nt,nz,ny,nx)),dims=dest_dims)

			if ndims_loop == 1:
				for k0 in np.arange(len(hremapped[variable].coords[hremapped[variable].dims[0]])):
					remapped[k0,:] = mod_akima_1d.mod_akima_1d.vertical_interpolation(zold,hremapped[variable][k0,:].values,znew)
					#outputds[variable][k0,:] = mod_akima_1d.mod_akima_1d.vertical_interpolation(zold,hremapped[variable][k0,:].values,znew)
					outputds.update({variable:remapped})
			#elif ndims_loop == 0:
			#	outputds[variable][:] = mod_akima_1d.mod_akima_1d.vertical_interpolation(zold,hremapped[variable][:].values,znew)


		return outputds


	def _mask_output(self,outputds,target_grid,list_variables):
		''' mask the interpolated field with land sea mask from target grid '''
		for variable in list_variables:
			ndims_mask = len(target_grid.lsm.dims)
			ndims_array = len(outputds[variable].dims)
			ndims_loop = ndims_array - ndims_mask

			if ndims_loop == 2:
				pass # TOdO
			elif ndims_loop == 1:
				for k0 in np.arange(len(outputds[variable].coords[outputds[variable].dims[0]])):
					outputds[variable][k0,:] = outputds[variable][k0,:] * target_grid.lsm[:]

		return outputds
