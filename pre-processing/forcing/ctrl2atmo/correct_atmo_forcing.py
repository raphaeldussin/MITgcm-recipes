from lib_correct_atmo_forcing import *

jradir='/t0/scratch/raphael/ASTE/jra55/'
jragrid='./jra_grid.nc'
dirxx='/t0/scratch/raphael/ASTE/config_AnTNguyen.v2/ADXXfiles/'
astegrid='/t0/scratch/raphael/ASTE/tmpdir_ASTE-Phy-REF01/'
outdir='/t0/scratch/raphael/ASTE/workdir_jra55+ctrl/'

lon_atm, lat_atm = read_forcing_grid(jragrid, lonvar='longitude', latvar='latitude')

# Spec humidity test
year=2002
nrecsxx=419
iter=55
dateref='2002-01-01 0:00:00'
startdate=str(year) + '-01-01 0:00:00' # 1:30:00 for precips
mult_fact=1.

# read control file for full period
ds_ctrl, grid_ctrl = read_ctrl_file(dirxx + 'xx_aqh', iter, nrecsxx, astegrid)

# regrid control onto JRA grid
ds_ctrl_regridded = regrid_to_forcing(lon_atm, lat_atm, grid_ctrl, ds_ctrl,'Qair')

# verif: write to file for visual check
ds_ctrl_regridded.to_netcdf(outdir + 'xx_spechum.nc')

# resample one year to 3hourly data/30 min for precips
ds_ctrl_regridded_resampled = resample_ctrl_to_forcing_dates(ds_ctrl_regridded, year, freq='3H')

# verif: write to file for visual check
ds_ctrl_regridded_resampled.to_netcdf(outdir + 'xx_spechum_resampled.nc')

# read original forcing
jrafile=jradir + 'jra55_spfh2m_2005'

ds_jra = read_forcing_from_binary(jrafile, 'Qair', lon_atm, lat_atm, 
                                  dateref=dateref, startdate=startdate,
                                  freq='3H', precision='single')

# verif: write to file for visual check
ds_jra.to_netcdf(outdir + 'jra55_spfh2m_2005.nc')

# shortcut
#import xarray as xr
#ds_jra = xr.open_dataset(outdir + 'jra55_spfh2m_2005.nc')
#ds_ctrl_regridded_resampled=xr.open_dataset(outdir + 'xx_spechum_resampled.nc')

# make correction
ds_jra_corrected = add_correction(ds_jra, ds_ctrl_regridded_resampled, mult_fact)

# verif: write to file for visual check
ds_jra_corrected.to_netcdf(outdir + 'jra55_spfh2m_2005_corrected.nc')


