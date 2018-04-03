import numpy as np
import subprocess as sp
import matplotlib.pylab as plt
import struct

class joinmds():

	def __init__(self,ddir,dirout,var,step=None):
		#self.ddir = '/Users/raphael/WORK/RUNS_MITGCM/tutorial_global_oce_latlon/results/'
		#self.dirout = '/Users/raphael/WORK/RUNS_MITGCM/tutorial_global_oce_latlon/results_joined/'
		self.ddir = ddir
		self.dirout = dirout
		self.var = var
		self.step = step
		#self.var = 'V'
		#self.step = 20
		if self.step is not None:
			self.basename = self.ddir + self.var + '.' + str(self.step).zfill(10)
		else:
			self.basename = self.ddir + self.var
		# Assuming that there will always be a first one
		metadata_1 = self.read_meta(1,1)
		# copy time and space dims into object attributes
		self.nrecords = int(metadata_1['nrecords'])
		self.ndims = metadata_1['nDims']
		for kdim in np.arange(self.ndims):
			exec("self.n" + str(kdim) + " = int(metadata_1['dimList'][" + str(3*kdim) + "])")
		return None

	def __call__(self):
		self.allocate_global_array()
		ntiles, ntiles0, ntiles1 = self.how_many_tiles()
		for ntile0 in np.arange(ntiles0):
			for ntile1 in np.arange(ntiles1):
				metadata = self.read_meta(ntile0+1,ntile1+1)
				for kdim in np.arange(self.ndims):
					exec("self.start" + str(kdim) + " = int(metadata['dimList'][" + str(3*kdim+1) + "])")
					exec("self.start" + str(kdim) + " = self.start" + str(kdim) + " - 1" ) # python indexing
					exec("self.end" + str(kdim) + " = int(metadata['dimList'][" + str(3*kdim+2) + "])")
					exec("self.end" + str(kdim) + " = self.end" + str(kdim) + "- 1" ) # python indexing

				data_raw = self.read_data(ntile0+1,ntile1+1)
				if self.ndims == 3:
					data = np.reshape(data_raw,(1,self.n2,self.end1-self.start1+1,self.end0-self.start0+1))
					self.joindata[:,self.start2:self.end2+1,self.start1:self.end1+1,self.start0:self.end0+1] = data
				elif self.ndims == 2:
					data = np.reshape(data_raw,(1,self.end1-self.start1+1,self.end0-self.start0+1))
					self.joindata[:,self.start1:self.end1+1,self.start0:self.end0+1] = data
		self.write_new_data()
		self.write_new_metadata()
		return None

	def read_data(self,ntilei=1,ntilej=1):
		''' Read meta file and copy its content into a dictionary '''
		# define filename
		fdata = self.basename + '.' + str(ntilei).zfill(3) + '.' + str(ntilej).zfill(3) + '.data'
		data_raw = np.fromfile(fdata, dtype='>f')
		return data_raw

	def read_meta(self,ntilei=1,ntilej=1):
		''' Read meta file and copy its content into a dictionary '''
		# define filename
		fmetadata = self.basename + '.' + str(ntilei).zfill(3) + '.' + str(ntilej).zfill(3) + '.meta'
		# load the raw metadata
		raw_metadata = self.load_txtfile(fmetadata)
		# transfer content into dictionary
		metadata = {}	
		for line in raw_metadata:
			# clean up symbols
			items2remove = ['=','{','[',',','}',']']
			infoline=line
			for item in items2remove:
				infoline = infoline.replace(item,' ')
			# split line into items
			items = infoline.split()
			# special case for dimList to allow to pass a list of values
			# else get one value
			if items[0] == 'dimList':
				metadata['dimList'] = []
				for item in items[1:]:
					metadata['dimList'].append(item)
			else:
				exec("metadata['" + items[0] + "'] = " + items[1] )
		return metadata

	def load_txtfile(self,fname):
		''' Load meta file and reformat its content into one attribute per line '''
		# open metadata file in read mode
		fid = open(fname,'r')
		# dump content
		lines = fid.readlines()
		# close file
		fid.close()
		# reformat into one long string
		raw_metadata = ''
		for line in lines:
			raw_metadata += line
		# reformat the metadata into lines separated by the semi-colons
		# so that each line is one attribute
		raw_metadata = raw_metadata.replace('\n','').replace(';','\n')
		raw_metadata = raw_metadata.splitlines()
		return raw_metadata

	def allocate_global_array(self):
		''' Allocate the global array '''
		if self.ndims == 3:
			self.joindata = np.empty((self.nrecords,self.n2,self.n1,self.n0),dtype='float32')
		elif self.ndims == 2:
			self.joindata = np.empty((self.nrecords,self.n1,self.n0),dtype='float32')
		print 'allocated', self.joindata.shape
		return None
			
	def how_many_tiles(self):
		''' find out how many tiles there is to join '''
		# how many tiles ?
		cmd = 'ls ' + self.basename + '*meta' + ' | wc -l'
		ntiles = sp.check_output(cmd,shell=True).replace('\n',' ')
		# layout
		cmd = 'ls ' + self.basename + '*meta'
		listtiles = sp.check_output(cmd,shell=True).replace('\n',' ').split()
		n0max = 1
		n1max = 1
		for tile in listtiles:
			froot = tile.replace('/',' ').split()[-1]
			if self.step is not None:
				n0,n1 = froot.replace('.',' ').split()[2:4]
			else:
				n0,n1 = froot.replace('.',' ').split()[1:3]
			n0max = max(n0max,int(n0))
			n1max = max(n1max,int(n1))
		return ntiles, n0max, n1max

	def write_new_data(self):
		froot = self.basename.replace('/',' ').split()[-1]
		fileout = self.dirout + froot + '.data'
		#self.joindata.tofile(fileout,format='>f')
		fid   = open(fileout, "wb")
		flatdata = self.joindata.flatten()
		for kk in np.arange(len(flatdata)):
			tmp = struct.pack('>f',flatdata[kk])
			fid.write(tmp)
		fid.close()
		return None

	def write_new_metadata(self):
		froot = self.basename.replace('/',' ').split()[-1]
		fileout = self.dirout + froot + '.meta'
		metadata_1 = self.read_meta(1,1)
		metadata_out = metadata_1
		for kdim in np.arange(self.ndims):
			metadata_out['dimList'][3*kdim+1] = 1
			metadata_out['dimList'][3*kdim+2] = metadata_out['dimList'][3*kdim]

		fid = open(fileout,'w')
		fid.write(" simulation = {'" + metadata_out['simulation'] + "'}; \n")
                fid.write(" nDims = [   " + str(metadata_out['nDims']) + " ]; \n")
		fid.write(" dimList = [ \n")
		fid.write("     " + str(metadata_out['dimList'][0]) + "," + str(metadata_out['dimList'][1]) + ',' + str(metadata_out['dimList'][2]) + ",\n")
		if self.ndims == 2:
			fid.write("     " + str(metadata_out['dimList'][3]) + "," + str(metadata_out['dimList'][4]) + ',' + str(metadata_out['dimList'][5]) + "\n")
		if self.ndims == 3:
			fid.write("     " + str(metadata_out['dimList'][3]) + "," + str(metadata_out['dimList'][4]) + ',' + str(metadata_out['dimList'][5]) + ",\n")
			fid.write("     " + str(metadata_out['dimList'][6]) + "," + str(metadata_out['dimList'][7]) + ',' + str(metadata_out['dimList'][8]) + "\n")
		fid.write(" ]; \n")
		fid.write(" dataprec = [ 'float32' ]; \n")
		fid.write(" nrecords = [     1 ]; \n")
		fid.close()
		return None

ddir = '/Users/raphael/WORK/RUNS_MITGCM/tutorial_global_oce_latlon/results/'
dirout = '/Users/raphael/WORK/RUNS_MITGCM/tutorial_global_oce_latlon/results_joined/'

test = joinmds(ddir,dirout,'U',20)
test()
test = joinmds(ddir,dirout,'V',20)
test()
		

for var in ['XC','YC','XG','YG','DXC','DYC','DXG','DYG','RAC','RAZ','RAW','Depth','hFacC','hFacS','hFacW']:
	metric = joinmds(ddir,dirout,var)
	metric()
	#plt.figure()
	#if metric.ndims == 3:
	#	plt.pcolormesh(metric.joindata[0,0,:,:])
	#elif metric.ndims == 2:
	#	plt.pcolormesh(metric.joindata[0,:,:])
	#plt.colorbar()
	#plt.title(var)
	#plt.show()



#plt.figure()
#plt.pcolormesh(test.joindata[0,0,:,:])
#plt.show()
