from lib_correct_atmo_forcing import *
import numpy as np
import xmitgcm
import datetime

jradir='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/'
jragrid='./jra_grid.nc'
dirxx='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/xx_effective/'
astegrid='/local/data/artemis/workspace/rdussin/ASTE/RUNS/ASTE-Phy-Release1/outputs/'
astegridnc='/local/data/artemis/workspace/rdussin/ASTE/GRID/nc/aste_grid.nc'
outdir='/local/data/artemis/workspace/rdussin/ASTE/JRA55/jra55/corrected_modelgrid/'

lon_atm, lat_atm = read_forcing_grid(jragrid, lonvar='longitude', latvar='latitude')

firstyear=2002
lastyear=2017
nrecsxx=419   # number of time records from xx_var.effective, infered from file size / 1350 / 270 / 4
optimcycle=55 # optimization cycle number
debug=False

# define dictionaries for atmo variables
temp = {'ctrl_name': 'xx_atemp.effective', 'jra_name': 'jra55_tmp2m_degC', 'ncname': 't2', 'sf': 1, 'freq':'3H'}
humi = {'ctrl_name': 'xx_aqh.effective', 'jra_name': 'jra55_spfh2m', 'ncname': 'q2', 'sf': 1, 'freq':'3H'}
radlw = {'ctrl_name': 'xx_lwdown.effective', 'jra_name': 'jra55_dlw', 'ncname': 'radlw', 'sf': -1, 'freq':'1.5H'}
radsw = {'ctrl_name': 'xx_swdown.effective', 'jra_name': 'jra55_dsw', 'ncname': 'radsw', 'sf': -1, 'freq':'1.5H'}
precip = {'ctrl_name': 'xx_precip.effective', 'jra_name': 'jra55_rain', 'ncname': 'precip', 'sf': 1, 'freq':'1.5H'}
uwind = {'ctrl_name': 'xx_uwind.effective', 'jra_name': 'jra55_u10m', 'ncname': 'u10', 'sf': 1, 'freq':'3H'}
vwind = {'ctrl_name': 'xx_vwind.effective', 'jra_name': 'jra55_v10m', 'ncname': 'v10', 'sf': 1, 'freq':'3H'}

vars_to_correct = [temp, humi, radlw, radsw, precip, uwind, vwind]

for var in vars_to_correct:
    print('working on variable', var)

    # read control file for full period
    print('read control (eager)')
    if var['ncname'] in ['precip', 'radlw', 'radsw']:
        dateref='2002-01-01 0:00:00'
    else:
        dateref='2002-01-01 0:00:00'

    ds_ctrl, grid_ctrl = read_ctrl_file(dirxx + var['ctrl_name'], optimcycle, nrecsxx, astegrid, 
                                        refdate=dateref, ncname=var['ncname'])
    
    # loop on years
    for year in np.arange(firstyear, lastyear+1):
        now = datetime.datetime.now().isoformat()
        print(now, ' working on year', year)
        if var['ncname'] in ['precip', 'radlw', 'radsw']:
            startdate=str(year) + '-01-01 1:30:00' # 1:30:00 for precips, radiative
        else:
            startdate=str(year) + '-01-01 0:00:00' # for all the rest

        # resample one year to 3hourly data/30 min for precips and rads
        print('resample control (eager)')
        ds_ctrl.load()
        ds_ctrl_resampled = resample_ctrl_to_forcing_dates(ds_ctrl, year, freq=var['freq'])
        
        # debug: write to file for visual check
        if debug:
            ds_ctrl_resampled.to_netcdf(outdir + var['ctrl_name'] + '_resampled.nc')
        
        # read original forcing
        print('read forcing (eager)')
        jrafile = jradir + var['jra_name'] + '_' + str(year)
        
        ds_jra = read_forcing_from_binary(jrafile, var['ncname'], lon_atm, lat_atm, 
                                          dateref=dateref, startdate=startdate,
                                          freq='3H', precision='single')
        # move to -180/180
        ds_jra_180 = MITgcm_recipes.regridding.geo_roll_360_to_180(ds_jra, 'lon')

        # debug: write to file for visual check
        if debug:
            ds_jra_180.to_netcdf(outdir + var['jra_name'] + '_180_' + str(year) + '.nc')

        # project onto ASTE grid
        ds_jra_modelgrid = regrid_to_model(ds_jra_180, astegridnc, var['ncname'])

        # debug: write to file for visual check
        if debug:
            ds_jra_modelgrid.to_netcdf(outdir + var['jra_name'] + '_ASTEgrid_' + str(year) + '.nc')
        
        # make correction
        print('make corrected forcing (lazy)')
        ds_jra_corrected = add_correction(ds_jra_modelgrid, ds_ctrl_resampled, var['sf'])
        
        # write to file for visual check
        ds_jra_corrected.to_netcdf(outdir + var['jra_name'] + '_it' + str(optimcycle) + 'xx_' + str(year) +'.nc' )

        # put in compact form
        # write to binary file
        print('write corrected forcing (eager)')
        fileout = outdir + var['jra_name'] + '_it' + str(optimcycle) + 'xx_ASTEgrid_' + str(year)
        md = xmitgcm.utils.get_extra_metadata(domain='aste', nx=270)
        # cut into facets
        print('cut into facets')
        facets = xmitgcm.utils.rebuild_llc_facets(ds_jra_corrected[var['ncname']], md)
        # load facets
        # create compact
        print('create compact')
        compact = xmitgcm.utils.llc_facets_2dtime_to_compact(facets, md)
        # write to binary
        print('write binary')
        xmitgcm.utils.write_to_binary(compact, fileout, dtype=np.dtype('f'))

        # clean up
        del ds_ctrl_resampled
        del ds_jra 
        del ds_jra_180 
        del ds_jra_modelgrid
        del ds_jra_corrected

