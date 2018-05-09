import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
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
			filein = dirin + filealt
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
			filein = dirin + filealt
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
			filein = dirin + filealt
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

	def plot_map_field(self,dirin,variable,dict_plt,frame=None,level=None,timestep=None,precision='single',filealt=None):
		''' plot map of 2d or 3d field '''
		XCs = self.read_2d_field(dirin,'XC')
		YCs = self.read_2d_field(dirin,'YC')
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
		if self.domain == 'ASTE':
			m = Basemap(projection='npaeqd',boundinglat=-0,lon_0=320,resolution='l')
			xx1,yy1 = m(XCs[0],YCs[0])
			C=m.contourf(xx1,yy1,np.ma.masked_values(FIELDs[0][:,:],0),ncontours,norm=norm,\
			cmap=dict_plt['colorbar'])
			plt.colorbar(C, norm=norm)
			xx3,yy3 = m(XCs[2],YCs[2])
			m.contourf(xx3,yy3,np.ma.masked_values(FIELDs[2][:,:],0),ncontours,norm=norm,\
			cmap=dict_plt['colorbar'])
			xx4,yy4 = m(XCs[3],YCs[3])
			m.contourf(xx4,yy4,np.ma.masked_values(FIELDs[3][:,:],0),ncontours,norm=norm,\
			cmap=dict_plt['colorbar'])
			xx5,yy5 = m(XCs[4],YCs[4])
			# occasionally we can run into issues with input files
			field4 = np.ma.masked_values(FIELDs[4][:,:],0)
			field4[np.where(XCs[4] == 0)] = 0
			field4 = np.ma.masked_values(field4,0)
			m.contourf(xx5,yy5,field4,ncontours,norm=norm,\
			cmap=dict_plt['colorbar'])
			m.fillcontinents(color='grey',lake_color='white')
			m.drawcoastlines()
			plt.show()

		return None

