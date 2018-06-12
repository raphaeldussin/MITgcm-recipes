from lib_read_obc import *

dir_obc_bling = '/Volumes/L2/ASTE/INPUTS/inputs_bgc_from_global_BLINGv2/'

aste = obc_regional_llc()

aste.read_2d_obc(dir_obc_bling,'o2_obc_south_monthly_global_bling_v2.bin',\
orientation='NS',precision='single')

aste.read_2d_obc(dir_obc_bling,'o2_obc_east_monthly_global_bling_v2.bin',\
orientation='EW',precision='single')
