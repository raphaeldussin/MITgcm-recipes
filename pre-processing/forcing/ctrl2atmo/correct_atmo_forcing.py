from lib_correct_atmo_forcing import *
import pandas as pd

jragrid='./jra_grid.nc'

dirin='/t0/scratch/raphael/ASTE/config_AnTNguyen.v2/ADXXfiles/'
astegrid='/t0/scratch/raphael/ASTE/tmpdir_ASTE-Phy-REF01/'

lon_atm, lat_atm = read_forcing_grid(jragrid, lonvar='longitude', latvar='latitude')
data_ctrl, grid_ctrl = read_ctrl_file(dirin + 'xx_aqh', 55, 419, astegrid)

data = regrid_to_forcing(lon_atm, lat_atm, grid_ctrl, data_ctrl)

ds = xr.Dataset({'temperature': (['time', 'lat', 'lon'], data)},
                 coords={'lon': (['lon'], lon_atm),
                         'lat': (['lat'], lat_atm),
                         'time': pd.date_range('2002-01-01', freq='14D', periods=419),
                         'reference_time': pd.Timestamp('2002-01-01')})

ds.to_netcdf('test_spechum.nc')
