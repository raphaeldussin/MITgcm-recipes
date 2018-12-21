from lib_correct_atmo_forcing import *

jradir='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/'
jragrid='./jra_grid.nc'
dirxx='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/xx_effective/'
astegrid='/local/data/artemis/workspace/rdussin/ASTE/RUNS/ASTE-Phy-Release1/outputs/'
outdir='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/corrected/'

lon_atm, lat_atm = read_forcing_grid(jragrid, lonvar='longitude', latvar='latitude')

firstyear=2017
lastyear=2017
dateref='2002-01-01 0:00:00'
nrecsxx=419   # number of time records from xx_var.effective, infered from file size / 1350 / 270 / 4
optimcycle=55 # optimization cycle number
debug=False

# define dictionaries for atmo variables
temp = {'ctrl_name': 'xx_atemp.effective', 'jra_name': 'jra55_tmp2m_degC', 'ncname': 't2', 'sf': 1, 'freq':'3H'}
humi = {'ctrl_name': 'xx_aqh.effective', 'jra_name': 'jra55_spfh2m', 'ncname': 'q2', 'sf': 1, 'freq':'3H'}
radlw = {'ctrl_name': 'xx_lwdown.effective', 'jra_name': 'jra55_dlw', 'ncname': 'radlw', 'sf': -1, 'freq':'30min'}
radsw = {'ctrl_name': 'xx_swdown.effective', 'jra_name': 'jra55_dsw', 'ncname': 'radsw', 'sf': -1, 'freq':'30min'}
precip = {'ctrl_name': 'xx_precip.effective', 'jra_name': 'jra55_rain', 'ncname': 'precip', 'sf': 1, 'freq':'30min'}
uwind = {'ctrl_name': 'xx_uwind.effective', 'jra_name': 'jra55_u10m', 'ncname': 'u10', 'sf': 1, 'freq':'3H'}
vwind = {'ctrl_name': 'xx_vwind.effective', 'jra_name': 'jra55_v10m', 'ncname': 'v10', 'sf': 1, 'freq':'3H'}

#vars_to_correct = [temp, humi, radlw, radsw, precip, uwind, vwind]
vars_to_correct = [temp, radsw]

for var in vars_to_correct:
    print('working on variable', var)

    # read control file for full period
    print('read control (eager)')
    ds_ctrl, grid_ctrl = read_ctrl_file(dirxx + var['ctrl_name'], optimcycle, nrecsxx, astegrid)
    
    # regrid control onto JRA grid
    print('regrid control (eager/parallel)')
    ds_ctrl_regridded = regrid_to_forcing(lon_atm, lat_atm, grid_ctrl, ds_ctrl, var['ncname'])
    
    # debug: write to file for visual check
    if debug:
        ds_ctrl_regridded.to_netcdf(outdir + var['ctrl_name'] + '_regridded.nc')

    # loop on years
    for year in np.arange(firstyear, lastyear+1):
        print('working on year', year)
        if var['ncname'] in ['precip', 'radlw', 'radsw']:
            startdate=str(year) + '-01-01 1:30:00' # 1:30:00 for precips, radiative
        else:
            startdate=str(year) + '-01-01 0:00:00' # for all the rest

        # resample one year to 3hourly data/30 min for precips and rads
        print('resample control (lazy)')
        ds_ctrl_regridded_resampled = resample_ctrl_to_forcing_dates(ds_ctrl_regridded, year, freq=var['freq'])
        
        # debug: write to file for visual check
        if debug:
            ds_ctrl_regridded_resampled.to_netcdf(outdir + var['ctrl_name'] + '_resampled.nc')
        
        # read original forcing
        print('read forcing (eager)')
        jrafile = jradir + var['jra_name'] + '_' + str(year)
        
        ds_jra = read_forcing_from_binary(jrafile, var['ncname'], lon_atm, lat_atm, 
                                          dateref=dateref, startdate=startdate,
                                          freq='3H', precision='single')
        
        # debug: write to file for visual check
        if debug:
            ds_jra.to_netcdf(outdir + var['jra_name'] + '_' + str(year) + '.nc')
        
        # make correction
        print('make corrected forcing (lazy)')
        ds_jra_corrected = add_correction(ds_jra, ds_ctrl_regridded_resampled, var['sf'])
        
        # debug: write to file for visual check
        if debug:
            ds_jra_corrected.to_netcdf(outdir + var['jra_name'] + '_it' + str(optimcycle) + 'xx_' + str(year) +'.nc' )

        # write to binary file
        print('write corrected forcing (eager)')
        fileout = var['jra_name'] + '_it' + str(optimcycle) + 'xx_' + str(year)
        write_to_binary(ds_jra_corrected, var['ncname'], outdir + fileout, precision='single')

