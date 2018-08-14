import xarray as xr
import numpy as np
import sys
import struct

class create_ic():

	def __init__(self):
		'''  '''
		return None

	def create_ic_regional_llc(self,variable,datasets_facets,fileout,timestep=None,nx=270,clim=False,\
		background=None,precision='single'):
		''' create the IC file from list of datasets '''
		# in LLC configuration, we need to set nx
		self.nx = nx

		# compute global dimensions
		self.nt,self.nz,self.ny = self._compute_dimensions(datasets_facets)

		self.output = self._create_output_array(clim=clim)
		self._fill_array(datasets_facets,variable,timestep=timestep,clim=clim)

		if background is not None:
			self.output = self.output + background

		self._write_ic_to_binary(self.output,fileout,precision=precision)

		return None


	def _compute_dimensions(self,datasets_facets):
		''' compute the total size of domain in y '''
		ny_total=0
		self.nylist = []
		for ds in datasets_facets:
			nx_facet = len(ds['x'])
			ny_facet = len(ds['y'])
			dims = [nx_facet,ny_facet]
			dims.remove(self.nx)
			ny_total = ny_total + dims[0]
			self.nylist.append(dims[0])

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
				self.output[0,:,:self.nyend[0],:]              = \
				datasets_facets[0][variable][timestep,:,:,:].values
				self.output[0,:,self.nyend[0]:self.nyend[1],:] = \
				datasets_facets[1][variable][timestep,:,:,:].values
				self.output[0,:,self.nyend[1]:self.nyend[2],:] = \
				self._reorder_array(datasets_facets[2][variable][timestep,:,:,:].values)
				self.output[0,:,self.nyend[2]:self.nyend[3],:] = \
				self._reorder_array(datasets_facets[3][variable][timestep,:,:,:].values)
			if self.nt == 1:
				self.output[0,:,:self.nyend[0],:]              = \
				datasets_facets[0][variable][0,:,:,:].values
				self.output[0,:,self.nyend[0]:self.nyend[1],:] = \
				datasets_facets[1][variable][0,:,:,:].values
				self.output[0,:,self.nyend[1]:self.nyend[2],:] = \
				self._reorder_array(datasets_facets[2][variable][0,:,:,:].values)
				self.output[0,:,self.nyend[2]:self.nyend[3],:] = \
				self._reorder_array(datasets_facets[3][variable][0,:,:,:].values)
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
