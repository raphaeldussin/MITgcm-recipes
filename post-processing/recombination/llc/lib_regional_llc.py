import numpy as np
import matplotlib.pylab as plt
import cartopy as cart
from matplotlib import cm, colors

class regional_llc():

	def __init__(self,domain='ASTE'):
		''' Building class object '''
		if domain == 'ASTE':
			self.domain = domain
			self.nx = 270
			self.nx_fn = np.array([self.nx,0,self.nx,self.nx,self.nx])
			self.ny_fn = np.array([450,0,270,180,450])
			self.ny_total = self.ny_fn.sum()
			self.ny_end = self.ny_fn.cumsum()
			self.nz = 50
		elif domain == 'ASTE_hr':
			self.nx = 1080
			self.nx_fn = np.array([self.nx,0,self.nx,self.nx,self.nx])
			self.ny_fn = np.array([1260,0,1080,540,1260])
			self.ny_total = self.ny_fn.sum()
			self.ny_end = self.ny_fn.cumsum()
			self.nz = 80
		else:
			pass

	def read_2d_field(self,dirin,variable,timestep=None,precision='single',filealt=None):
		''' read 2d field from .data file '''
		if timestep is None:
			filein = dirin + '/' + variable + '.data'
		else:
			filein = dirin + '/' + variable + '.' + str(timestep).zfill(10) + '.data'
		# overide filein if filealt exists (usefull for input binaries,...)
		if filealt is not None:
			filein = dirin + '/' + filealt
		# load data
		if precision == 'single':
			data_raw = np.fromfile(filein,'>f')
		elif precision == 'double':
			data_raw = np.fromfile(filein,'>d')

		# reshape on dims
		data_compact = np.reshape(data_raw,[self.ny_total,self.nx])
		del data_raw
		# split in faces
		data_f1 = data_compact[:self.ny_end[0],:]
		data_f2 = data_compact[self.ny_end[0]:self.ny_end[1],:]
		data_f3 = data_compact[self.ny_end[1]:self.ny_end[2],:]
		data_f4 = np.reshape(data_compact[self.ny_end[2]:self.ny_end[3],:],[self.nx_fn[3],self.ny_fn[3]])
		data_f5 = np.reshape(data_compact[self.ny_end[3]:self.ny_end[4],:],[self.nx_fn[4],self.ny_fn[4]])

		return [data_f1,data_f2,data_f3,data_f4,data_f5]

	def read_3d_field(self,dirin,variable,timestep=None,precision='single',filealt=None):
		''' read 3d field from .data file '''
		if timestep is None:
			filein = dirin + '/' + variable + '.data'
		else:
			filein = dirin + '/' + variable + '.' + str(timestep).zfill(10) + '.data'
		if filealt is not None:
			filein = dirin + '/' + filealt
		# load data
		if precision == 'single':
			data_raw = np.fromfile(filein,'>f')
		elif precision == 'double':
			data_raw = np.fromfile(filein,'>d')

		# reshape on dims
		data_compact = np.reshape(data_raw,[self.nz,self.ny_total,self.nx])
		del data_raw
		# split in faces
		data_f1 = data_compact[:,:self.ny_end[0],:]
		data_f2 = data_compact[:,self.ny_end[0]:self.ny_end[1],:]
		data_f3 = data_compact[:,self.ny_end[1]:self.ny_end[2],:]
		data_f4 = np.reshape(data_compact[:,self.ny_end[2]:self.ny_end[3],:],[self.nz,self.nx_fn[3],self.ny_fn[3]])
		data_f5 = np.reshape(data_compact[:,self.ny_end[3]:self.ny_end[4],:],[self.nz,self.nx_fn[4],self.ny_fn[4]])

		return [data_f1,data_f2,data_f3,data_f4,data_f5]

	def read_4d_field(self,dirin,variable,ntimes,timestep=None,precision='single',filealt=None):
		''' read 4d field from .data file '''
		if timestep is None:
			filein = dirin + '/' + variable + '.data'
		else:
			filein = dirin + '/' + variable + '.' + str(timestep).zfill(10) + '.data'
		if filealt is not None:
			filein = dirin + '/' + filealt
		# load data
		if precision == 'single':
			data_raw = np.fromfile(filein,'>f')
		elif precision == 'double':
			data_raw = np.fromfile(filein,'>d')

		# reshape on dims
		data_compact = np.reshape(data_raw,[ntimes,self.nz,self.ny_total,self.nx])
		del data_raw
		# split in faces
		data_f1 = data_compact[:,:,:self.ny_end[0],:]
		data_f2 = data_compact[:,:,self.ny_end[0]:self.ny_end[1],:]
		data_f3 = data_compact[:,:,self.ny_end[1]:self.ny_end[2],:]
		data_f4 = np.reshape(data_compact[:,:,self.ny_end[2]:self.ny_end[3],:],\
		[ntimes,self.nz,self.nx_fn[3],self.ny_fn[3]])
		data_f5 = np.reshape(data_compact[:,:,self.ny_end[3]:self.ny_end[4],:],\
		[ntimes,self.nz,self.nx_fn[4],self.ny_fn[4]])

		return [data_f1,data_f2,data_f3,data_f4,data_f5]

	def plot_map_field(self,dirin,variable,dict_plt,frame=None,level=None,timestep=None,\
		                precision='single',precision_xy='single',filealt=None,dirxy=None):
		''' plot map of 2d or 3d field '''
		if dirxy is None:
			dirxy = dirin
		XCs = self.read_2d_field(dirxy,'XC',precision=precision_xy)
		YCs = self.read_2d_field(dirxy,'YC',precision=precision_xy)
		Depth = self.read_2d_field(dirxy,'Depth',precision=precision_xy)
		if frame is None:
			if level is None:
				FIELDs = self.read_2d_field(dirin,variable,timestep=timestep,\
				precision=precision,filealt=filealt)
			else:
				tmp = self.read_3d_field(dirin,variable,timestep=timestep,\
				precision=precision,filealt=filealt)
				FIELDs = []
				for kk in np.arange(5):
					FIELDs.append(tmp[kk][level,:,:])

		# do the actual plot
		plt.figure(figsize=[10,8])
		ncontours=45 # number of contours could be moved to input dict_plt
		norm = colors.Normalize(vmin=dict_plt['vmin'], vmax=dict_plt['vmax'])
		if 'mult_fact' in dict_plt:
			mult_fact = dict_plt['mult_fact']
		else:
			mult_fact = 1

		if self.domain in ['ASTE','ASTE_hr']:
			m = plt.axes(projection=cart.crs.Orthographic(central_longitude=-45, central_latitude=60))
			# facet1
			fieldplt1 = self.mask_invalid_data(Depth[0],FIELDs[0])
			C = m.pcolormesh(XCs[0],YCs[0],mult_fact*fieldplt1,norm=norm,\
			cmap=dict_plt['colorbar'],transform=cart.crs.PlateCarree())
			plt.colorbar(C, norm=norm)
			# facet2 is empty
			# facet3
			fieldplt3 = self.mask_invalid_data(Depth[2],FIELDs[2])
			m.pcolormesh(XCs[2],YCs[2],mult_fact*fieldplt3,norm=norm,\
			cmap=dict_plt['colorbar'],transform=cart.crs.PlateCarree())
			# facet4
			fieldplt4 = self.mask_invalid_data(Depth[3],FIELDs[3])
			m.pcolormesh(XCs[3],YCs[3],mult_fact*fieldplt4,norm=norm,\
			cmap=dict_plt['colorbar'],transform=cart.crs.PlateCarree())
			# facet5
			fieldplt5 = self.mask_invalid_data(Depth[4],FIELDs[4])
			m.pcolormesh(XCs[4],YCs[4],mult_fact*fieldplt5,norm=norm,\
			cmap=dict_plt['colorbar'],transform=cart.crs.PlateCarree())
			m.coastlines()
			m.add_feature(cart.feature.LAND, facecolor='0.75')
			gl = m.gridlines(draw_labels=False)
			plt.show()

		return None

	def mask_invalid_data(self,depth,field):
		mask = np.full(depth.shape, False)
		mask[np.where(depth == 0)] = True
		fieldout = np.ma.masked_values(field,0)
		fieldout = np.ma.masked_array(data=field,mask=mask)
		return fieldout
