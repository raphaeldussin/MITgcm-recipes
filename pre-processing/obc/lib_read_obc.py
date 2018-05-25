import numpy as np
import matplotlib.pylab as plt
#from mpl_toolkits.basemap import Basemap
from matplotlib import cm, colors

class obc_regional_llc():

	def __init__(self,domain='ASTE'):
		''' Building class object '''
		if domain == 'ASTE':
			self.domain = domain
			self.nx_fn = np.array([270,0,270,270,270])
			self.nyEW_fn = np.array([450,0,270,270,270])
			self.nyNS_fn = np.array([270,0,270,180,450])
			self.nyEW_total = self.nyEW_fn.sum()
			self.nyNS_total = self.nyNS_fn.sum()
			self.nyEW_end = self.nyEW_fn.cumsum()
			self.nyNS_end = self.nyNS_fn.cumsum()
			self.nz = 50
		else:
			pass

	def read_2d_obc(self,dirin,obcfile,step=0,orientation='NS',precision='single',ntimes=12):
		filein = dirin + obcfile
		self.ntime = ntimes
		# load data
		if precision == 'single':
			data_raw = np.fromfile(filein,'>f')
		elif precision == 'double':
			data_raw = np.fromfile(filein,'>d')

		# reshape on dims
		if orientation == 'NS':
			data_compact = np.reshape(data_raw,[self.ntime,self.nz,self.nyNS_total])
		elif orientation == 'EW':
			data_compact = np.reshape(data_raw,[self.ntime,self.nz,self.nyEW_total])

		del data_raw
		## split in faces
		if orientation == 'NS':
			data_f1 = data_compact[:,:,:self.nyNS_end[0]]
			data_f2 = data_compact[:,:,self.nyNS_end[0]:self.nyNS_end[1]]
			data_f3 = data_compact[:,:,self.nyNS_end[1]:self.nyNS_end[2]]
			data_f4 = data_compact[:,:,self.nyNS_end[2]:self.nyNS_end[3]]
			data_f5 = data_compact[:,:,self.nyNS_end[3]:self.nyNS_end[4]]
		elif orientation == 'EW':
			data_f1 = data_compact[:,:,:self.nyEW_end[0]]
			data_f2 = data_compact[:,:,self.nyEW_end[0]:self.nyEW_end[1]]
			data_f3 = data_compact[:,:,self.nyEW_end[1]:self.nyEW_end[2]]
			data_f4 = data_compact[:,:,self.nyEW_end[2]:self.nyEW_end[3]]
			data_f5 = data_compact[:,:,self.nyEW_end[3]:self.nyEW_end[4]]

		del data_compact

		#plt.figure()
		#plt.contourf(np.ma.masked_values(data_compact[0,:,:],0))
		plt.figure()
		plt.subplot(151)
		if data_f1.shape[-1] != 0:
			plt.contourf((np.ma.masked_values(data_f1[step,:],0)))
		plt.subplot(152)
		if data_f2.shape[-1] != 0:
			plt.contourf((np.ma.masked_values(data_f2[step,:],0)))
		plt.subplot(153)
		if data_f3.shape[-1] != 0:
			plt.contourf((np.ma.masked_values(data_f3[step,:],0)))
		plt.subplot(154)
		if data_f4.shape[-1] != 0:
			plt.contourf((np.ma.masked_values(data_f4[step,:],0)))
		plt.subplot(155)
		if data_f5.shape[-1] != 0:
			plt.contourf((np.ma.masked_values(data_f5[step,:],0)))
		plt.show()

		return None #[data_f1,data_f2,data_f3,data_f4,data_f5]
