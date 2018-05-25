from lib_read_obc import *

aste = obc_regional_llc()
aste.read_2d_obc('/Users/raphael/WORK/ASTE/obc_ASTE/physics/','OBSt_1170x50x290_10Jul2017_bl.bin',\
orientation='NS',precision='single')

#aste.read_2d_obc('/Users/raphael/WORK/ASTE/obc_ASTE/physics/','OBNt_1170x50x290_10Jul2017_bl.bin',\
#orientation='NS',precision='single')

aste.read_2d_obc('/Users/raphael/WORK/ASTE/obc_ASTE/physics/','OBEt_1260x50x290_10Jul2017_bl.bin',\
orientation='EW',precision='single')

aste.read_2d_obc('/Users/raphael/WORK/ASTE/obc_ASTE/physics/','OBWt_1260x50x290_10Jul2017_bl.bin',\
orientation='EW',precision='single')
