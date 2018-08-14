import xarray as xr
import numpy as np
import sys
import struct
import extract_obc

class create_obc():

	def __init__(self,namelist='data.obcs'):
		''' args : full path to data.obcs used '''
		self.namelist = namelist
		return None

	def create_obc_regional_llc(self,variable,boundary,datasets_facets,fileout,precision='single',repeats=0):
		''' create the OBC file for one boundary from list of datasets '''
		# Parse data from data.obcs
		namelist_info = self._parse_namelist(boundary)

		# Compute length of OBC vector from N datasets combined
		ntotal = self._verif_size_combined_facets(boundary,datasets_facets)

		# Sanity Check
		if ntotal != namelist_info['ntotal']:
			message = 'datasets combined length do not match data.obcs \ntotal found = ' + \
			str(ntotal) + ' instead of ' + str(namelist_info['ntotal'])
			exit(message)
		else:
			print('working on boundary', boundary, 'of length', ntotal)
			pass

		datasets_facets = self._create_global_indices(boundary,datasets_facets)

		obc = self._extract_obc_from_datasets(variable,boundary,datasets_facets,namelist_info)

		self._write_obc_to_binary(obc,fileout,precision=precision,repeats=repeats)

		return None


	def _parse_namelist(self,boundary):
		if boundary == 'north':
			keyword = 'OB_Jnorth'
			#idx_length=0 ; idx_value=1
		elif boundary == 'south':
			keyword = 'OB_Jsouth'
		elif boundary == 'east':
			keyword = 'OB_Ieast'
		elif boundary == 'west':
			keyword = 'OB_Iwest'
		else:
			print('Unkown boundary condition')
			pass

		fid = open(self.namelist,'r')
		lines = fid.readlines()
		fid.close()
		for line in lines:
			if line.find(keyword) != -1:
				useful_line = line

		#print(useful_line)
		segment_length = []
		segment_value  = []
		segments = useful_line.replace(keyword,'').replace('=','').replace(',','').split()
		#print(segments)
		for segment in segments:
			segment_length.append(int(segment.replace('*',' ').split()[0]))
			segment_value.append(int(segment.replace('*',' ').split()[1]))

		segment_length = np.array(segment_length) ; segment_value = np.array(segment_value)
		ntotal = segment_length.sum() ; nsegments = len(segment_length)
		boundary_vector_global = np.arange(ntotal)
		boundary_on_off_global = 9999 * np.ones(ntotal,dtype='i')

		seg_start=0
		for k in np.arange(nsegments):
			seg_end = seg_start + segment_length[k]
			boundary_on_off_global[seg_start:seg_end] = segment_value[k]
			seg_start = seg_end

		namelist_info = {'boundary_vector_global':boundary_vector_global,'boundary_on_off_global':boundary_on_off_global,\
		'ntotal':ntotal}
		return namelist_info

	def _verif_size_combined_facets(self,boundary,datasets_facets):
		len_total=0
		if boundary in ['north','south']:
			for ds in datasets_facets:
				len_total = len_total + len(ds['x'])
		elif boundary in ['east','west']:
			for ds in datasets_facets:
				len_total = len_total + len(ds['y'])
		return len_total

	def _create_global_indices(self,boundary,datasets_facets):
		if boundary in ['north','south']:
			offset=0
			for dataset in datasets_facets:
				dataset['xglo'] = dataset['x'] + offset
				offset = offset + len(dataset['x'])
		elif boundary in ['east','west']:
			offset=0
			for dataset in datasets_facets:
				dataset['yglo'] = dataset['y'] + offset
				offset = offset + len(dataset['y'])

		return datasets_facets

	def _extract_obc_from_datasets(self,variable,boundary,datasets_facets,namelist_info):
		ndims = len(datasets_facets[0][variable].dims) - 2
		if ndims == 2:
			n0 = datasets_facets[0][variable][:].shape[0] # ugly
			n1 = datasets_facets[0][variable][:].shape[1]
			output = xr.DataArray(9999 * np.ones((n0,n1,namelist_info['ntotal'])),dims=('time','z','obc'))

			if boundary in ['north','south']:
				for kdata in np.arange(len(datasets_facets)):
					tmp = datasets_facets[kdata].copy()
					tmp['x'] = tmp['xglo']
					for ji in tmp['x'].values:
						if namelist_info['boundary_on_off_global'][ji] == 0: # no obc value
							output[:,:,ji] = 0
						else:
							output[:,:,ji] = tmp[variable].sel(x=ji,y=namelist_info['boundary_on_off_global'][ji]-1)
			elif boundary in ['east','west']:
				for kdata in np.arange(len(datasets_facets)):
					#print('facet #', kdata)
					tmp = datasets_facets[kdata].copy()
					tmp['y'] = tmp['yglo']
					jstart = tmp['yglo'].values.min()
					jend = tmp['yglo'].values.max() + 1
					#print(jstart,jend)
					#print(namelist_info['boundary_on_off_global'][jstart:jend])
					#print(type(namelist_info['boundary_on_off_global'][jstart:jend]))
					output[:,:,jstart:jend] = extract_obc.meridional_boundary_tz(tmp[variable].values,namelist_info['boundary_on_off_global'][jstart:jend])
#					for jj in tmp['y'].values:
#						print(jj)
#						if namelist_info['boundary_on_off_global'][jj] == 0: # no obc value
#							output[:,:,jj] = 0
#							print(' set to zero')
#						else:
#							output[:,:,jj] = tmp[variable].sel(y=jj,x=namelist_info['boundary_on_off_global'][jj]-1)
#							print(tmp[variable].dims)
#							print('take data at',   namelist_info['boundary_on_off_global'][jj]-1)

		return output

	def _write_obc_to_binary(self,obc,fileout,precision='single',repeats=0):
		# write data to binary files
		fid   = open(fileout, "wb")
		flatdata = obc.values.flatten()
		for kt in np.arange(1+repeats):
			for kk in np.arange(len(flatdata)):
				if precision == 'single':
					tmp = struct.pack('>f',flatdata[kk])
				elif precision == 'double':
					tmp = struct.pack('>d',flatdata[kk])
				fid.write(tmp)
		fid.close()
		return None
