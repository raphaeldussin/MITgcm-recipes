import numpy as np
import xarray as xr
import matplotlib.pylab as plt

class mitgcm_grid():

        def __init__(self,domain='ASTE',tile=1):
                ''' Building class object '''

		self.tile=tile
		self.grid_var_list = ['XC','YC','DXF','DYF','RAC',\
		                      'XG','YG','DXV','DYU','RAZ',\
                                      'DXC','DYC','RAW','RAS','DXG','DYG']

                if domain == 'ASTE':
                        self.domain = domain
                        self.nx = 270
                        self.nx_fn = np.array([self.nx,self.nx,self.nx,180,450]) + 1
                        self.ny_fn = np.array([450,180,270,self.nx,self.nx]) + 1
                        self.nx_fn_interior = np.array([self.nx,self.nx,self.nx,180,450])
                        self.ny_fn_interior = np.array([450,180,270,self.nx,self.nx])

                        self.ny_total = self.ny_fn.sum()
                        self.ny_end = self.ny_fn.cumsum()
			self.ny_total_interior = self.ny_fn_interior.sum()
			self.ny_end_interior = self.ny_fn_interior.cumsum()
			self.nxy = self.nx * self.ny_total

			self.delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,\
			                      10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,\
			                      31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,\
			                      93.96, 96.58, 98.25, 99.25,100.01,101.33,104.56,111.33,122.83,\
			                      139.09,158.94,180.83,203.55,226.50,249.50,272.50,295.50,318.50,\
			                      341.50,364.50,387.50,410.50,433.50,456.50])
                        self.nz = len(self.delR)
			self.Z = self.delR.cumsum()

                elif domain == 'ASTE_hr':
                        self.nx = 1080
                        self.nx_fn = np.array([self.nx,0,self.nx,540,1260]) + 1
                        self.ny_fn = np.array([1260,0,1080,self.nx,self.nx]) + 1
                        self.ny_total = self.ny_fn.sum()
                        self.ny_end = self.ny_fn.cumsum()
                        self.nz = 80
                else:
                        pass

		Z = xr.DataArray(self.Z,coords=[np.arange(self.nz)],dims=['z'])
		dz = xr.DataArray(self.delR,coords=[self.Z],dims=['depth'])
		self.grid = xr.Dataset({'Z':Z,'dz':dz})

		return None

	def load_grid_from_bin(self,gridfile='tileNTILE.mitgrid',precision='double'):
		''' populate grid object from data read in binary files tile000*.mitgrid '''
		# read all grid file
		data_raw = self._read_raw_data_from_bin(gridfile=gridfile,precision=precision)
		# extract individual variables
		#for varname in self.grid_var_list:
		#	tmp = self._extract_variable_from_raw_data(varname,data_raw,tile=tile)
		#	exec('self.' + varname + ' = tmp' )
		rawXC = self._extract_variable_from_raw_data('XC',data_raw)
		rawYC = self._extract_variable_from_raw_data('YC',data_raw)
		rawDXC = self._extract_variable_from_raw_data('DXC',data_raw)
		rawDYC = self._extract_variable_from_raw_data('DYC',data_raw)

		# grid is asymetric but stored symetric with a padding of one
		XC = xr.DataArray(rawXC[:-1,:-1], coords=[np.arange(self.ny_fn[self.tile-1]-1),np.arange(self.nx_fn[self.tile-1]-1)], dims=['y', 'x'])
		YC = xr.DataArray(rawYC[:-1,:-1], coords=[np.arange(self.ny_fn[self.tile-1]-1),np.arange(self.nx_fn[self.tile-1]-1)], dims=['y', 'x'])

		self.grid.update({'XC': XC,'YC': YC})

		return None

	def _read_raw_data_from_bin(self,gridfile='tileNTILE.mitgrid',precision='double'):
		''' load all data from binary file '''
		gridfile = gridfile.replace('NTILE',str(self.tile).zfill(3))
		print('loading grid from file' + gridfile)
                if precision == 'single':
                        data_raw = np.fromfile(gridfile,'>f')
                elif precision == 'double':
                        data_raw = np.fromfile(gridfile,'>d')
		return data_raw

	def _extract_variable_from_raw_data(self,varname,data_raw):
		''' extract the chosen variable from the sequence of variables '''
		index_var = self.grid_var_list.index(varname)
		if index_var == -1:
			print('variable not found')
			pass
		else:
			#print('nx =' , self.nx_fn[self.tile-1] ,', ny =' , self.ny_fn[self.tile-1])
			nxy= self.nx_fn[self.tile-1] * self.ny_fn[self.tile-1]
			kstart=index_var*nxy
			kend=kstart+nxy
			data_extract = data_raw[kstart:kend]
		# reshape the variable
		data_out = np.reshape(data_extract,[self.ny_fn[self.tile-1],self.nx_fn[self.tile-1]])
		return data_out

	def write_grid_to_nc(self,filename):
		self.grid.to_netcdf(filename)
		return None

        def infer_mask_from_output_file(self,filein,precision='single',spval=0,skiptiles=[2]):
                ''' read 3d output .data file and infers mask from spval '''
                # load data
                if precision == 'single':
                        data_raw = np.fromfile(filein,'>f')
                elif precision == 'double':
                        data_raw = np.fromfile(filein,'>d')

		nx_fn = self.nx_fn_interior.copy()
		ny_fn = self.ny_fn_interior.copy()
		for tile in skiptiles:
			nx_fn[tile-1] = 0
			ny_fn[tile-1] = 0

		nxy = nx_fn * ny_fn
		nxy_end = nxy.cumsum()
		nxy_start = nxy_end - nxy
                # reshape on vertical dimension
                data_compact = np.reshape(data_raw,[self.nz,-1])
                del data_raw

		n1 = nx_fn[self.tile-1]
		n2 = ny_fn[self.tile-1]

               	data = data_compact[:,nxy_start[self.tile-1]:nxy_end[self.tile-1]]
		data_r =np.reshape(data,[self.nz,n2,n1])

		mask = np.ones(data_r.shape)
		mask[np.where(data_r == spval)] = 0

		lsm = xr.DataArray(mask, coords=[np.arange(self.nz),np.arange(n2),np.arange(n1)], dims=['z','y','x'])
		self.grid.update({'lsm': lsm})

		#plt.figure()
		#plt.pcolormesh(mask[0,:,:])
		#plt.show()

		return None

