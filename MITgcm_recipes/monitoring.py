import pandas as pd
import glob
import subprocess as sp
import datetime as dt
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as mdates

def make_monitor_file(dirlogs):
    ''' create a csv file with content from monitor package that
    pandas can read '''
    filenames = glob.glob(dirlogs + 'STDOUT.0000.*')
    sp.call('rm monitor.log', shell=True)
    for filename in filenames:
        sp.call("cat " + filename + " | grep %MON | awk '{print $4,$6}' >> monitor.log",shell=True)
    return None

def read_monitoring(dirlogs):
    ''' create monitor.log and read with pandas '''
    make_monitor_file(dirlogs)
    monitoring = pd.read_csv('monitor.log',sep=' ',names=['gcmvariable','gcmvalue'])
    # create the time axis
    timesteps = monitoring[monitoring.gcmvariable=='time_tsnumber'].gcmvalue.values
    seconds_from = monitoring[monitoring.gcmvariable=='time_secondsf'].gcmvalue.values
    # start date can be found in STDOUT looking for startdate
    datestr=sp.check_output('grep startDate_1 STDOUT.0000.1', shell=True).replace('startDate_1=',' ').split()[-1].replace(',','')
    timestr=sp.check_output('grep startDate_2 STDOUT.0000.1', shell=True).replace('startDate_2=',' ').split()[-1].replace(',','')

    year=datestr[:4]
    month=datestr[4:6]
    day=datestr[6:8]
    hour=timestr[:2]
    mins=timestr[2:4]
    secs=timestr[4:6]

    dates= []
    for kt in np.arange(len(seconds_from)):
        dates.append(dt.datetime(year,month,day,hour,mins,secs) + dt.timedelta(seconds=seconds_from[kt]))

    monitoring['dates'] = dates

    return monitoring

# monitoring plot
def plot_statistics_of(cvar,monitoring):
    ''' plots min/max/mean timeseries,
    takes variable without suffix as string (e.g. dynstat_sst)'''
    meanvar = monitoring[monitoring.gcmvariable==cvar+'_mean'].gcmvalue.values
    minvar  = monitoring[monitoring.gcmvariable==cvar+'_min' ].gcmvalue.values
    maxvar  = monitoring[monitoring.gcmvariable==cvar+'_max' ].gcmvalue.values

    fig = plt.figure(figsize=[12,8])

    plt.subplot(121)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator())
    plt.fill_between(monitoring['dates'],minvar,maxvar,color='grey')
    plt.plot(monitoring['dates'],meanvar,'r.')
    plt.title('range and mean '+cvar)
    plt.gcf().autofmt_xdate()

    plt.subplot(122)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
    plt.gca().xaxis.set_minor_locator(mdates.DayLocator())
    plt.plot(monitoring['dates'],meanvar,'r.')
    plt.title('zoom on mean '+cvar)
    plt.gcf().autofmt_xdate()
    plt.show()
