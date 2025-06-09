#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:09:40 2022

@author: caro
"""

##########################################################################
# global initilization & pre-processed tasks #############################

# python module import ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
from oggm import cfg, utils, workflow, tasks
import logging
import numpy as np
import os
import pandas

# oggm initialization ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cfg.initialize(logging_level='WARNING')
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['continue_on_error'] = True
cfg.PARAMS['run_mb_calibration'] = True
cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['border'] = 80
cfg.PARAMS['store_fl_diagnostics'] = True

cfg.PARAMS["climate_qc_months"] = 3 # hugonnet 
cfg.PARAMS["hydro_month_nh"] = 1 # hugonnet
cfg.PARAMS["hydro_month_sh"] = 1 # hugonnet
cfg.PARAMS["max_mu_star"] = 600

cfg.PARAMS["temp_all_solid"] = 0.0
cfg.PARAMS["temp_all_liq"] = 2.0 
cfg.PARAMS["temp_melt"] = 0

cfg.PARAMS['cfl_min_dt'] = 10 # para simular glaciares con problemas

# # glacier call #############################
# datos_rgi = pandas.read_csv('/home/caro/04_Andes_2000_2019_CR2_2023/datos/RGI_BNA_Clusters.csv')
# # call parameters
# datos_param = pandas.read_csv('/home/caro/04_Andes_2000_2019_CR2_2023/datos/LR_Pf.csv')
# glacier call #############################
datos_rgi = pandas.read_csv('/Users/milliespencer/Desktop/files_chile_OGGM_climate_comparison/cautin_RGI_BNA_Clusters.csv')
# call parameters
datos_param = pandas.read_csv('/Users/milliespencer/Desktop/files_chile_OGGM_climate_comparison/LR_Pf.csv')

# # INICIO ITERACION por region GLACIO +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# listas de llamado
list_region = ['OT3','DA1','DA2','DA3','WA1','WA2','WA3','WA4','WA5','WA6']


for zona in list_region:
    
# directorio
#    zona = 'OT3'
    salida = '/Users/milliespencer/Desktop/files_chile_OGGM_climate_comparison' + zona + '/'
    cfg.PATHS['working_dir'] = salida                                 

# extraccion de glaciares por zona
    gdf_sel = datos_rgi[datos_rgi['Cluster'] == zona]
    gdf_sel = list(gdf_sel.RGIId)
    
# se extrae LT usando id glaciar
    param = datos_param[datos_param['Cluster'] == zona]
    LT_n = float(param.LR)
# se extrae Pf 
    Pf_n = float(param.Pf)

    # definion LT en modelo
    cfg.PARAMS["temp_default_gradient"] = LT_n
    cfg.PARAMS["prcp_scaling_factor"] = Pf_n
    cfg.PARAMS['use_winter_prcp_factor'] = False
    
        
###
### LLAMADO PERIODO 2000-2019 +++++++++++++++++++++++++++++++++++++++++++++++++
###
    # Guardado statistics
    output_folder= salida
    rgi_version  = '_'+ zona
    border = 80
    output_base_dir = os.path.join(output_folder,'RGI{}'.format(rgi_version), 'b_{:03d}'.format(border))
    sum_dir = os.path.join(output_base_dir, 'L3', 'summary')
    utils.mkdir(sum_dir)

    # Level = 2 2000-2019+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#    gdirs = workflow.init_glacier_directories(gdf_sel, from_prepro_level=2)
    
    gdirs = workflow.init_glacier_directories(gdf_sel, from_prepro_level=2,
                                              prepro_base_url='https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L1-L2_files/centerlines/',
                                              prepro_border=80)
    
       # Climate
    cfg.PARAMS['baseline_climate'] = 'CR2MET25'
    from oggm.shop import cr2met_25
    workflow.execute_entity_task(cr2met_25.process_cr2met_25_data, gdirs) # ,output_filesuffix='_cr2met'
    utils.get_geodetic_mb_dataframe()  # Small optim to avoid concurrency
    workflow.execute_entity_task(tasks.mu_star_calibration_from_geodetic_mb, gdirs,ref_period='2000-01-01_2020-01-01') # edit climate.py Between lines 539 and 554
    workflow.execute_entity_task(tasks.apparent_mb_from_any_mb, gdirs)
    
    filter = border >= 20
    workflow.calibrate_inversion_from_consensus(gdirs,apply_fs_on_mismatch=True,error_on_mismatch=False,filter_inversion_output=filter)
         
        # We get ready for modelling
    log = logging.getLogger(__name__)     
    if border >= 20:
                workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)
    else:
        log.workflow('L3: for map border values < 20, wont initialize glaciers '
                      'for the run.')
    
    compilacion = '_2000_2019_hydro_TC_'+ zona 
    # save outputs 2000-2019
    workflow.execute_entity_task (tasks.run_with_hydro, gdirs,
                         ys=1999,     
                         run_task = tasks.run_from_climate_data,          # The task to run with hydro
                         store_monthly_hydro = True,                      # compute monthly hydro diagnostics
                         ref_area_from_y0 = True,                           # Even if the glacier may grow, keep the reference area as the year 0 of the simulation
                         output_filesuffix = compilacion,                # Where to write the output - this is needed to stitch the runs together afterwards
                         );

    # Guardado glacier statistics 
    rgi_reg = '_' + zona
    opath = os.path.join(sum_dir, 'glacier_statistics_{}.csv'.format(rgi_reg))
    utils.compile_glacier_statistics(gdirs, path=opath)
    mb_dif = utils.compile_glacier_statistics(gdirs, path=opath)
    opath = os.path.join(sum_dir, 'climate_statistics_{}.csv'.format(rgi_reg))
    utils.compile_climate_statistics(gdirs, path=opath)
    opath = os.path.join(sum_dir, 'fixed_geometry_mass_balance_{}.csv'.format(rgi_reg))
    utils.compile_fixed_geometry_mass_balance(gdirs, path=opath)
    
    # compilación datos en un archivo
    ds2000 = utils.compile_run_output(gdirs, input_filesuffix= compilacion )
