from lib_read_obc import *

aste = obc_regional_llc()
aste.read_2d_obc('./','alk_south_bling.bin',\
orientation='NS',precision='single')

aste.read_2d_obc('./','alk_east_bling.bin',\
orientation='EW',precision='single')

