from lib_read_obc import *

dir_obc_bling = '/Volumes/L2/ASTE/INPUTS/inputs_cfc_from_glodap+carmack/'

aste = obc_regional_llc()

aste.read_2d_obc(dir_obc_bling,'CFC11_obc_south_monthly_glodap+carmack.bin',\
orientation='NS',precision='single',ntimes=1)

aste.read_2d_obc(dir_obc_bling,'CFC11_obc_east_monthly_glodap+carmack.bin',\
orientation='EW',precision='single',ntimes=1)

aste.read_2d_obc(dir_obc_bling,'CFC12_obc_south_monthly_glodap+carmack.bin',\
orientation='NS',precision='single',ntimes=1)

aste.read_2d_obc(dir_obc_bling,'CFC12_obc_east_monthly_glodap+carmack.bin',\
orientation='EW',precision='single',ntimes=1)
