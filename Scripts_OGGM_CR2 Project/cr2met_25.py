#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import numpy as np
import pandas as pd
from scipy import stats

# Optional libs
try:
    import salem
except ImportError:
    pass

# Locals
from oggm import cfg
from oggm import utils
from oggm import entity_task
from oggm.exceptions import InvalidParamsError

# Module logger
log = logging.getLogger(__name__)

CR2MET_SERVER = '/home/caro/00_climate_input/'


def set_cr2met_url(url):
    """If you want to use a different server for CR2MET (for testing, etc)."""
    global CR2MET_SERVER
    CR2MET_SERVER = url


@utils.locked_func
def get_cr2met_25_file(var=None):
    """Returns a path to the desired CR2MET baseline climate file.

    If the file is not present, download it.

    Parameters
    ----------
    var : str
        'tmp' for temperature
        'pre' for precipitation

    Returns
    -------
    str
        path to the file
    """

    # Be sure input makes sense
    if var not in ['tmp', 'pre']:
        raise InvalidParamsError('cr2met variable {} '
                                 'does not exist!'.format(var))

    # File to look for
    if var == 'tmp':
        bname = 'CR2met_t2m_hgt_2022_1960_dic_2021_2.5.nc'
    else:
        bname = 'CR2met_pr_2022_1960_dic_2021_2.5.nc'
                 
    h_url = CR2MET_SERVER + bname #+ '.bz2'
    #return utils.file_extractor(utils.file_downloader(h_url))
    return h_url


@entity_task(log, writes=['climate_historical'])
def process_cr2met_25_data(gdir, y0=None, y1=None, output_filesuffix=None):
    """Processes and writes the CR2MET baseline climate data for this glacier.

    Extracts the nearest timeseries and writes everything to a NetCDF file.

    Parameters
    ----------
    gdir : :py:class:`oggm.GlacierDirectory`
        the glacier directory to process
    y0 : int
        the starting year of the timeseries to write. The default is to take
        1850 (because the data is quite bad before that)
    y1 : int
        the ending year of the timeseries to write. The default is to take
        the entire time period available in the file, but with this kwarg
        you can shorten it (to save space or to crop bad data)
    output_filesuffix : str
        this add a suffix to the output file (useful to avoid overwriting
        previous experiments)
    """

    if cfg.PARAMS['baseline_climate'] != 'CR2MET25': # change HISTALP por cr2met
        raise InvalidParamsError("cfg.PARAMS['baseline_climate'] should be "
                                 "set to CR2MET25.")


    # read the time out of the pure netcdf file
    ft = get_cr2met_25_file('tmp')
    fp = get_cr2met_25_file('pre')
    
    with utils.ncDataset(ft) as nc:
        vt = nc.variables['time']
        assert vt[0] == 1
        assert vt[-1] == vt.shape[0]
        t0 = vt.units.split(' since ')[1][:7] # comienzo '1960-01'
        time_t = pd.date_range(start=t0, periods=vt.shape[0], freq='MS')
        
    with utils.ncDataset(fp) as nc:
        vt = nc.variables['time']
        #vt = utils.ncDataset(fp).variables['time']
        assert vt[0] == 1
        assert vt[-1] == vt.shape[0]
        t0 = vt.units.split(' since ')[1][:7]
        time_p = pd.date_range(start=t0, periods=vt.shape[0], freq='MS')

    # now open with salem
    nc_ts_tmp = salem.GeoNetcdf(ft, time=time_t)
    nc_ts_pre = salem.GeoNetcdf(fp, time=time_p)

    # some default
    if y0 is None:
        y0 = 1960 # comienzo cr2met
        
    # set temporal subset for the ts data (hydro years)
    # the reference time is given by precip, which is shorter
    sm = cfg.PARAMS['hydro_month_' + gdir.hemisphere]
    em = sm - 1 if (sm > 1) else 12
    yrs = nc_ts_pre.time.year
    y0 = yrs[0] if y0 is None else y0
    y1 = yrs[-1] if y1 is None else y1

    nc_ts_tmp.set_period(t0='{}-{:02d}-01'.format(y0, sm),
                         t1='{}-{:02d}-01'.format(y1, em))
    nc_ts_pre.set_period(t0='{}-{:02d}-01'.format(y0, sm),
                         t1='{}-{:02d}-01'.format(y1, em))
    time = nc_ts_pre.time
    ny, r = divmod(len(time), 12)
    assert r == 0
    
    # units
    assert nc_ts_tmp._nc.variables['hgt'].units.lower() in ['m', 'meters',
                                                              'meter',
                                                              'metres',
                                                              'metre']
    assert nc_ts_tmp._nc.variables['temp'].units.lower() in ['degc', 'degrees',
                                                             'degrees celcius',
                                                             'degree', 'c']
    assert nc_ts_pre._nc.variables['prcp'].units.lower() in ['kg m-2',
                                                                 'l m-2', 'mm',
                                                                 'millimeters',
                                                                 'millimeter']
  
    # geoloc
    lon = gdir.cenlon
    lat = gdir.cenlat
    nc_ts_tmp.set_subset(corners=((lon, lat), (lon, lat)), margin=1)
    nc_ts_pre.set_subset(corners=((lon, lat), (lon, lat)), margin=1)   
    
    # read the data
    temp = nc_ts_tmp.get_vardata('temp')
    prcp = nc_ts_pre.get_vardata('prcp')
    hgt = nc_ts_tmp.get_vardata('hgt')
    ref_lon = nc_ts_tmp.get_vardata('lon')
    ref_lat = nc_ts_tmp.get_vardata('lat')
    source = 'CR2MET25'
    nc_ts_tmp._nc.close()
    nc_ts_pre._nc.close()

    # should we compute the gradient?
    use_grad = cfg.PARAMS['temp_use_local_gradient']
    igrad = None
    if use_grad:
        igrad = np.zeros(len(time)) * np.NaN
        for t, loct in enumerate(temp):
            slope, _, _, p_val, _ = stats.linregress(hgt.flatten(),
                                                     loct.flatten())
            igrad[t] = slope if (p_val < 0.01) else np.NaN

    gdir.write_monthly_climate_file(time, prcp[:, 1, 1], temp[:, 1, 1],
                                    hgt[1, 1], ref_lon[1], ref_lat[1],
                                    filesuffix=output_filesuffix,
                                    source=source)

