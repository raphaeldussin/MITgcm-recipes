import xarray as xr
import numpy as np
import sys
import struct

class create_ic():

	def __init__(self):
		'''  '''
		return None

	def create_ic_regional_llc(self,variable,datasets_facets,fileout,timestep=None,nx=270,clim=False):
		''' create the IC file from list of datasets '''
		# in LLC configuration, we need to set nx
		self.nx = nx

		# Parse data from data.obcs
		#namelist_info = self._parse_namelist(boundary)

		# Compute length of OBC vector from N datasets combined
		#ntotal = self._verif_size_combined_facets(boundary,datasets_facets)

		# Sanity Check
		#if ntotal != namelist_info['ntotal']:
		#	message = 'datasets combined length do not match data.obcs \ntotal found = ' + \
		#	str(ntotal) + ' instead of ' + str(namelist_info['ntotal'])
		#	exit(message)
		#else:
		#	print('working on boundary', boundary, 'of length', ntotal)
		#	pass

		#datasets_facets = self._create_global_indices(boundary,datasets_facets)

		#obc = self._extract_obc_from_datasets(variable,boundary,datasets_facets,namelist_info)

		#self._write_obc_to_binary(obc,fileout)

		# compute global dimensions
		self.nt,self.nz,self.ny = self._compute_dimensions(datasets_facets)

		print(self.nx)
		print(self.ny)

		self.output = self._create_output_array(clim=clim)
		self._fill_array(datasets_facets,variable,timestep=timestep,clim=clim)

		self._write_ic_to_binary(self.output,fileout,precision='single')

		return None


	def _compute_dimensions(self,datasets_facets):
		''' compute the total size of domain in y '''
		# we can't assume facets are given in the right order so we need
		# to figure out which dimension correspond to the MITgcm y coord
		# the y of the facet is the good one for facets 1,2,3 but not for 4,5
		ny_total=0
		self.nylist = []
		for ds in datasets_facets:
			nx_facet = len(ds['x'])
			ny_facet = len(ds['y'])
			dims = [nx_facet,ny_facet]
			dims.remove(self.nx)
			ny_total = ny_total + dims[0]
			self.nylist.append(dims[0])

		print(self.nylist)
		self.nylist = np.array(self.nylist)
		for ds in datasets_facets:
			nz = len(ds['z'])

		for ds in datasets_facets:
			nt = len(ds['time'])

		return nt,nz,ny_total

	def _create_output_array(self,clim=False):
		if clim:
			placeholder = np.empty((self.nt,self.nz,self.ny,self.nx))
		else:
			placeholder = np.empty((1,self.nz,self.ny,self.nx))
		return placeholder


	def _fill_array(self,datasets_facets,variable,timestep=None,clim=False):
		# quick and dirty
		self.nyend = self.nylist.cumsum()
		if not clim:
			if self.nt > 1:
				self.output[0,:,:self.nyend[0],:]               = datasets_facets[0][variable][timestep,:,:,:].values
				self.output[0,:,self.nyend[0]:self.nyend[1],:] = datasets_facets[1][variable][timestep,:,:,:].values
				self.output[0,:,self.nyend[1]:self.nyend[2],:] = self._reorder_array(datasets_facets[2][variable][timestep,:,:,:].values)
				self.output[0,:,self.nyend[2]:self.nyend[3],:] = self._reorder_array(datasets_facets[3][variable][timestep,:,:,:].values)
		return None


	def _reorder_array(self,datain):
		''' reorder the major axis '''
		ndims = len(datain.shape)
		if ndims == 4:
			nt,nz,ny,nx = datain.shape
			dataout = np.empty((nt,nz,nx,ny))
			for kt in np.arange(nt):
				for kz in np.arange(nz):
					dataout[kt,kz,:,:] = np.reshape(datain[kt,kz,:,:],[nx,ny])
		if ndims == 3:
			nz,ny,nx = datain.shape
			dataout = np.empty((nz,nx,ny))
			for kz in np.arange(nz):
				dataout[kz,:,:] = np.reshape(datain[kz,:,:],[nx,ny])

		return dataout

#	def _parse_namelist(self,boundary):
#		if boundary == 'north':
#			keyword = 'OB_Jnorth'
#			#idx_length=0 ; idx_value=1
#		elif boundary == 'south':
#			keyword = 'OB_Jsouth'
#		elif boundary == 'east':
#			keyword = 'OB_Ieast'
#		elif boundary == 'west':
#			keyword = 'OB_Iwest'
#		else:
#			print('Unkown boundary condition')
#			pass
#
#		fid = open(self.namelist,'r')
#		lines = fid.readlines()
#		fid.close()
#		for line in lines:
#			if line.find(keyword) != -1:
#				useful_line = line
#
#		#print(useful_line)
#		segment_length = []
#		segment_value  = []
#		segments = useful_line.replace(keyword,'').replace('=','').replace(',','').split()
#		#print(segments)
#		for segment in segments:
#			segment_length.append(int(segment.replace('*',' ').split()[0]))
#			segment_value.append(int(segment.replace('*',' ').split()[1]))
#
#		segment_length = np.array(segment_length) ; segment_value = np.array(segment_value)
#		ntotal = segment_length.sum() ; nsegments = len(segment_length)
#		boundary_vector_global = np.arange(ntotal)
#		boundary_on_off_global = 9999 * np.ones(ntotal)
#
#		seg_start=0
#		for k in np.arange(nsegments):
#			seg_end = seg_start + segment_length[k]
#			boundary_on_off_global[seg_start:seg_end] = segment_value[k]
#			seg_start = seg_end
#
#		namelist_info = {'boundary_vector_global':boundary_vector_global,'boundary_on_off_global':boundary_on_off_global,\
#		'ntotal':ntotal}
#		return namelist_info
#
#	def _verif_size_combined_facets(self,boundary,datasets_facets):
#		len_total=0
#		if boundary in ['north','south']:
#			for ds in datasets_facets:
#				len_total = len_total + len(ds['x'])
#		elif boundary in ['east','west']:
#			for ds in datasets_facets:
#				len_total = len_total + len(ds['y'])
#		return len_total
#
#	def _create_global_indices(self,boundary,datasets_facets):
#		if boundary in ['north','south']:
#			offset=0
#			for dataset in datasets_facets:
#				dataset['xglo'] = dataset['x'] + offset
#				offset = offset + len(dataset['x'])
#		elif boundary in ['east','west']:
#			offset=0
#			for dataset in datasets_facets:
#				dataset['yglo'] = dataset['y'] + offset
#				offset = offset + len(dataset['y'])
#
#		return datasets_facets
#
#	def _extract_obc_from_datasets(self,variable,boundary,datasets_facets,namelist_info):
#		ndims = len(datasets_facets[0][variable].dims) - 2
#		if ndims == 2:
#			n0 = datasets_facets[0][variable][:].shape[0] # ugly
#			n1 = datasets_facets[0][variable][:].shape[1]
#			output = xr.DataArray(9999 * np.ones((n0,n1,namelist_info['ntotal'])),dims=('time','z','obc'))
#
#			if boundary in ['north','south']:
#				for kdata in np.arange(len(datasets_facets)):
#					tmp = datasets_facets[kdata].copy()
#					tmp['x'] = tmp['xglo']
#					for ji in tmp['x'].values:
#						if namelist_info['boundary_on_off_global'][ji] == 0: # no obc value
#							output[:,:,ji] = 0
#						else:
#							output[:,:,ji] = tmp[variable].sel(x=ji,y=namelist_info['boundary_on_off_global'][ji]-1)
#			elif boundary in ['east','west']:
#				for kdata in np.arange(len(datasets_facets)):
#					tmp = datasets_facets[kdata].copy()
#					tmp['y'] = tmp['yglo']
#					for jj in tmp['y'].values:
#						if namelist_info['boundary_on_off_global'][jj] == 0: # no obc value
#							output[:,:,jj] = 0
#						else:
#							output[:,:,jj] = tmp[variable].sel(y=jj,x=namelist_info['boundary_on_off_global'][jj]-1)
#							pass
#
#		return output

	def _write_ic_to_binary(self,ic,fileout,precision='single'):
		# write data to binary files
		fid   = open(fileout, "wb")
		flatdata = ic.flatten()
		for kk in np.arange(len(flatdata)):
			if precision == 'single':
				tmp = struct.pack('>f',flatdata[kk])
			elif precision == 'double':
				tmp = struct.pack('>d',flatdata[kk])
			fid.write(tmp)
		fid.close()
		return None
