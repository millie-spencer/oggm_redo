#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:09:40 2022

@author: caro
"""

##########################################################################
# global initilization & pre-processed tasks #############################

# python module import ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
from oggm import cfg, utils, workflow, tasks, graphics
import matplotlib.pyplot as plt
import logging
import oggm
from oggm import cfg, utils, workflow, tasks, graphics, global_tasks
from oggm.core import massbalance, flowline
from oggm.shop import gcm_climate
from oggm.core.massbalance import ScalarMassBalance
from oggm.tests.funcs import bu_tidewater_bed
from oggm.core.flowline import FluxBasedModel

import pandas as pd
from oggm.core.massbalance import MultipleFlowlineMassBalance # mass-balance at the glacier level

import xarray as xr
import numpy as np
import pandas as pd
import os
import geopandas as gpd

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

# glacier call #############################
datos_rgi = pd.read_csv('/Users/milliespencer/Desktop/CR2_OGGM_Paper/files_chile_OGGM_climate_comparison/RGI_BNA_Clusters.csv')
# call parameters
datos_param = pd.read_csv('/Users/milliespencer/Desktop/CR2_OGGM_Paper/files_chile_OGGM_climate_comparison/LR_Pf.csv')


# # INICIO ITERACION por region GLACIO +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# listas de llamado
list_region = ['OT3','DA1','DA2','DA3','WA1','WA2','WA3','WA4','WA5','WA6']


for zona in list_region:
    
# directorio
 #   zona = 'OT3'
    salida = '/Users/milliespencer/Desktop/CR2_OGGM_Paper/' + zona + '/'
    cfg.PATHS['working_dir'] = salida                                 
    compilacion_zona = salida + 'run_output_2000_2019_hydro_TC_' + zona + '.nc'

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

# entradas iniciales     
    zonax = xr.open_dataset(compilacion_zona)
    rgi_id = zonax['rgi_id'].values
    area = zonax['area'].values

# area total por clusters desde RGI 6

    rgi_area = datos_rgi[datos_rgi['Cluster'] == zona]
    rgi_area = rgi_area.Area.sum()
    
# llamado de simulaciones por zonas glacio, 13 en total

# zona = 'OT3' 
# directorio = '/home/caro/03_Andes_2000_2019_TC/'
# llamado_zona = directorio + zona + '/'
# compilacion_zona = llamado_zona + 'run_output_2000_2019_hydro_TC_' + zona + '.nc'

# zonax = xr.open_dataset(compilacion_zona)
# rgi_id = zonax['rgi_id'].values
# area = zonax['area'].values


    # calculo area usando simulacion
            # revisar numero de filas, deben ser 21.
    
    ID_area = pd.DataFrame(area)
    ID_rgi = pd.DataFrame(rgi_id)
    ID_rgi.rename(columns={ID_rgi.columns[0]:'id' }, inplace = True)
    
    # ver desde que año comienza sim area,
    # debería ser 2000, pero algoa pasa
    ID_area = pd.DataFrame(ID_area.iloc[10])#.transpose() # sacar area
    
    ID_area.rename(columns={ID_area.columns[0]:'area' }, inplace = True)
    ID_area = pd.concat([ID_rgi.reset_index(drop=True), ID_area], axis=1) # merge two dataframes
    ID_area = ID_area.dropna(axis=0) # delere na rows
    max_n = ID_area[ID_area.columns[0]].count()
    ID_area['number'] = range(0,max_n)
    ID_data = ID_area.number.values.tolist()
    ID_rgi_id = ID_area.id.values.tolist()
    
    #     # calculo area usando RGI
        
    # ID_rgi = pd.DataFrame(rgi_id)
    # ID_rgi.rename(columns={ID_rgi.columns[0]:'id' }, inplace = True)
    
    # # llamo hugonnet areas mismas de rgi
    # gmb = utils.get_geodetic_mb_dataframe()
    # gmb.to_csv('/home/caro/OGGM_Maipo_terraclimate_2_2000_2100/gmb_hugonnet.csv')
    # gmb = gmb.loc[gmb['period'] == '2000-01-01_2020-01-01'] # extraigo solo '2000-01-01_2020-01-01'
    # gmb['rgiid2'] = gmb.index 
    # gmb = gmb[gmb['rgiid2'].isin(ID_rgi['id'] )]
    
    # ID_area = pd.DataFrame(gmb['area'])
    # ID_area['rgi_id'] = ID_area.index 
    
    # max_n = ID_area[ID_area.columns[0]].count()
    # ID_area['number'] = range(0,max_n)
    # ID_data = ID_area.number.values.tolist()
    # ID_rgi_id = ID_area.rgi_id.values.tolist()
    
    
    
        # mass balance ++++++++++++++++++++++++
    
    # cfg.PATHS['working_dir'] = llamado_zona 
    gdirs = workflow.init_glacier_directories(ID_rgi_id)
    
    # # parametros
    # cfg.PARAMS['prcp_scaling_factor']=2.5
    # cfg.PARAMS['temp_default_gradient'] = -0.0066
    
    # probando si se puesde sacar mb en base a todos los glaciare
    # y no solo los que tienen area sim el año 2000
    
              # SMB MAIPO
              
    lista = []
    for i in ID_data:
                x = MultipleFlowlineMassBalance(gdirs[i],use_inversion_flowlines=True)
                lista.append(x)
    
    years = np.arange(2000, 2022)
    mb = []
    for i in lista:
                xx = i.get_specific_mb(year=years) # Specific mb for this year and a specific glacier geometry. Units: [mm w.e. yr-1]
                mb.append(xx)
    mb_per_glacier = pd.DataFrame(mb) # mm/year
    # asociacion con id
    mb_per_glacier_id = mb_per_glacier
    mb_per_glacier_id['rgi_id'] = ID_rgi_id
    c = salida + 'mb.csv' # zona
    mb_per_glacier_id.to_csv(c,index=False)
    
        
    # area
    areas = ID_area.area.values.tolist()     
    areas  = np.divide(areas, 1000000)
    # mb promedio poderada por area 
    mb_mean = pd.DataFrame()
    itera = np.arange(0,22).tolist()
        
    for i in itera:
            mb_list = mb_per_glacier[i].tolist()
            avg = np.average(mb_list, axis=0, weights = areas)
            data = {'itera': [i],
                  'mb': [avg]
                }
            df = pd.DataFrame(data)
            mb_mean = pd.concat([mb_mean, df], axis=0)
        
    mb_mean['year'] =  np.arange(2000, 2022)
    c = salida + 'mb_promedio.csv' # zona
    mb_mean.to_csv(c,index=False)
        
        # # geodetic mass balance
        # # geodetic mass balance
        # # geodetic mass balance
        
    gmb = utils.get_geodetic_mb_dataframe()
    gmb.to_csv('/home/caro/OGGM_Maipo_terraclimate_2_2000_2100/gmb_hugonnet.csv')
    gdf_sel_RGIId_list = ID_rgi_id #gdf_sel['RGIId'].tolist() # lista glaciares maipo
    lista = pd.DataFrame(gdf_sel_RGIId_list)
    c = salida + 'lista_rgi.csv' # zona
    lista.to_csv(c,index=False)
        
        # comparación GMB y SMB
        
    gmb_maipo = gmb.loc[gmb['period'] == '2000-01-01_2020-01-01'] # extraigo solo '2000-01-01_2020-01-01'
    gmb_maipo['dmdtda_mm'] = gmb_maipo['dmdtda']*1000 # m a mm
    gmb_maipo['err_dmdtda_mm'] = gmb_maipo['err_dmdtda']*1000 # m a mm
    gmb_maipo['rgiid2'] = gmb_maipo.index 
    gmb_maipo = gmb_maipo[gmb_maipo['rgiid2'].isin(ID_rgi_id)]
    gmb_maipo['dmdtda_mm_pon'] = gmb_maipo.dmdtda_mm*((gmb_maipo.area/1000000)/(gmb_maipo.area.sum()/1000000))
    gmb_maipo['err_dmdtda_mm_pon'] = gmb_maipo.err_dmdtda_mm*((gmb_maipo.area/1000000)/(gmb_maipo.area.sum()/1000000))
    
    gmb = (gmb_maipo.dmdtda_mm_pon.sum(),"GMB") # -271
    gmb_error = (gmb_maipo.err_dmdtda_mm_pon.sum(),"GMB_error") # -271
    
    # ponderar mb por area y promedir glaciares 
        # preparación mb sim 
    mb_sim = mb_per_glacier_id.drop(mb_per_glacier_id.columns[20], axis=1) # eliminar columna id 
    mb_sim = pd.DataFrame(mb_sim.mean(axis=1)) # mean por row
    mb_sim.rename(columns={mb_sim.columns[0]:'mb_sim'}, inplace = True)
    mb_sim['rgiid'] = mb_per_glacier_id['rgi_id'] # agregar columna id
    mb_sim['area_sim'] = areas
    mb_sim['mb_mm_pon'] = mb_sim.mb_sim*((mb_sim.area_sim)/(mb_sim.area_sim.sum()))
    
    smb = (mb_sim.mb_mm_pon.sum(),"SMB")             #    
    
    rgi_area_s = (round(rgi_area,3),'area_rgi') # area rgi
    rgi_oggm_s = (round(((ID_area.area.sum())/1000000),3),'area_rgi') # area rgi
    porcentaje_s = (round(((rgi_oggm_s[0]/rgi_area_s[0])*100),3),'por_area_oggm')
    n_g_rgi = (len(gdf_sel),'n_g_rgi')
    n_g_oggm = (max_n,'n_g_oggm')
    
    comp_gmb_smb = []
    comp_gmb_smb.append(gmb)
    comp_gmb_smb.append(gmb_error)
    comp_gmb_smb.append(smb)
    comp_gmb_smb.append(rgi_area_s)
    comp_gmb_smb.append(rgi_oggm_s)
    comp_gmb_smb.append(porcentaje_s)
    comp_gmb_smb.append(n_g_rgi)
    comp_gmb_smb.append(n_g_oggm)
    
    comp_gmb_smb = pd.DataFrame(comp_gmb_smb)
    
    c = salida + 'comparacion_gmb_smb_' + zona + '_.csv' # zona
    comp_gmb_smb.to_csv(c,index=False)
    
    
    #    plot mb maipo +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    a_gmb = np.array([(float(gmb[0]))/1000]*22)
    a_error = np.array(([(float(gmb_error[0]))/1000]*22)) 
    
    rgi_area = round(rgi_area,1) # area rgi
    oggm_area = round(((ID_area.area.sum())/1000000),1) # area procesada por oggm
    
    porcetaje = round(((oggm_area/rgi_area)*100),1) 
    #mean_smb = str(round(float(smb[0]),1))
    mean_gmb = str(round(float(gmb[0]),1))
    mean_smb = mean_gmb
    oggm_area = str(oggm_area)
    porcetaje = str(porcetaje)
    
    
    from pylab import figure, text
    plt.rcParams["figure.figsize"] = (7,4)
    name =  'Annual_mb_'+zona+'_2000-2021.png' # save the figure
    pathname = salida + name
    
    f = figure()
    ax = f.add_subplot(111)
    plt.plot(mb_mean.year,mb_mean.mb/1000, color="blue", label="OGGM MB")# set x-axis label
    plt.plot(years, ([(float(gmb[0]))/1000]*22), color="black", label="GMB")# set x-axis label
    plt.fill_between(years, a_gmb-a_error, a_gmb+a_error, alpha=0.3, color='black',label="GMB error");
    plt.legend(loc="lower left")
    #L.get_texts()[3].set_text('WGMS')
    plt.xlabel('Time [yr]',fontsize=15)
    plt.ylabel("smb [mm w.e./yr]",color="black",fontsize=15)# set y-axis label
    plt.yticks(fontsize=15) 
    plt.xticks(fontsize=15)#,rotation=90)
    plt.xticks(np.arange(2000, 2022, 5))
    plt.title(zona,fontsize=20)
    plt.tight_layout()
    
    text(0.7, 0.92,'A = '+ oggm_area+' '+"km\u00b2", ha='left', va='center', transform=ax.transAxes,fontsize=13)
    text(0.7, 0.84,'A sim = '+ porcetaje+'%', ha='left', va='center', transform=ax.transAxes,fontsize=13)
    text(0.7, 0.76,'GMB = '+ mean_gmb, ha='left', va='center', transform=ax.transAxes,fontsize=13)
    text(0.7, 0.68,'SMB = '+ mean_smb, ha='left', va='center', transform=ax.transAxes,fontsize=13)
        
    plt.savefig(pathname, format = 'png', dpi=300)
    plt.show()
    
    #
    # Temp and P
    #
    # iteracion extracción variable por 
    lista_temp = []
    lista_prcp = []
    lista_ele = []
    lista_temp = pd.DataFrame(lista_temp)
    lista_prcp = pd.DataFrame(lista_prcp)
    lista_ele = pd.DataFrame(lista_ele)
 
    itera = np.arange(0,max_n).tolist()

    for i in itera:
       clima = gdirs[i].get_filepath('climate_historical')
       ds = xr.open_dataset(clima)
       prcp_x = pd.DataFrame(ds.prcp) # star on 1950-01-01
       temp_x = pd.DataFrame(ds.temp)
       
       id_rgi = gdirs[i]
       id_rgi = id_rgi.rgi_id
       ele_x = {id_rgi:[ds.ref_hgt]}
       ele_x = pd.DataFrame(ele_x)
       
       lista_prcp = pd.concat([lista_prcp,prcp_x],axis=1) 
       lista_temp = pd.concat([lista_temp,temp_x],axis=1) 
       lista_ele = pd.concat([lista_ele,ele_x],axis=1) 
           
    ct = salida + 'lista_temp.csv' # zona
    cp = salida + 'lista_prec.csv' # zona
    ce = salida + 'lista_ele.csv' # zona
    
    lista_temp.to_csv(ct,index=False)       
    lista_prcp.to_csv(cp,index=False)
    lista_ele.to_csv(ce,index=False)
