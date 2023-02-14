#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:58:54 2021
Create Pandas dataframes out of TM5-4DVAR flux data
@author: eschoema
"""

datapath = "."
savepath = "."

import datetime
#from datetime import timedelta
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas
import xarray as xr
import time

def CreateDF(DS, filepath):
    '''
    # function returning a pandas dataframe with the fluxes in the given xarray dataset
    # arguments:
    #           DS: xarray dataset
    #           filepath: where to save the dataframe
    '''
    Flux_ocn = []
    Flux_bio = []
    Flux_ff = []
    Flux_fire = []
    Flux = []
    region = []
    Month = []
    Year = []
    Date =[]
    
    for i in range(len(DS.field.values[0])):# loop through all timesteps (dimensions [category, times, region])
        #check for monthly fluxes 
        # get number of days of current timestep by subtracting starttime from endtime
        dnum = (datetime.date(DS.time_end.values[i][0],DS.time_end.values[i][1],DS.time_end.values[i][2]) 
            - datetime.date(DS.time_start.values[i][0],DS.time_start.values[i][1],DS.time_start.values[i][2])).days
        if dnum <=31 and dnum >= 28:
            for x,j in enumerate(list(range(16))+[2,5,7,11,13,43]):
                if x <= 15 or j == 43:
                    Flux_ocn.append(DS.field.values[0][i][j])
                    Flux_bio.append(DS.field.values[1][i][j])
                    Flux_ff.append(DS.field.values[2][i][j])
                    Flux_fire.append(DS.field.values[3][i][j])
                    Flux.append(DS.field.values[4][i][j])
                else:
                    Flux_ocn.append(DS.field.values[0][i][j]+DS.field.values[0][i][j+1])
                    Flux_bio.append(DS.field.values[1][i][j]+DS.field.values[1][i][j+1])
                    Flux_ff.append(DS.field.values[2][i][j]+DS.field.values[2][i][j+1])
                    Flux_fire.append(DS.field.values[3][i][j]+DS.field.values[3][i][j+1])
                    Flux.append(DS.field.values[4][i][j]+DS.field.values[4][i][j+1])
                    
                Date.append(datetime.date(DS.time_start.values[i][0],DS.time_start.values[i][1],15))
                region.append(x)
                Month.append(DS.time_start.values[i][1])
                Year.append(DS.time_start.values[i][0])

    # if current file contains the prior values, add _apri to parameter name in dataframe
    kind = ''            
    if 'APRI' in filepath.split('_')[-1]:
        kind = '_apri'

    d = {'Flux_ocn' + kind: Flux_ocn, 
         'Flux_bio' + kind: Flux_bio, 
         'Flux_ff' + kind:Flux_ff,
         'Flux_fire' + kind:Flux_fire,
         'Flux' + kind:Flux,
         'Region': region, 
         'Year': Year, 
         'Month': Month, 
         'Date': Date}
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameTM5Flux(dataset,res):
    '''
    # function saving a pandas dataframe with the TM5-4DVar fluxes for all chosen regions in CreateDF() for the given dataset
    # arguments:
    #           dataset: which TM5-4DVar Flux dataset (which data is assimilated by TM5-4DVar) e.g. 'RemoTeCISloc_Flux','IS_Flux','ACOSIS_Flux'
    #           res: 'glb3x2' or 'glb6x4'
    '''

    print("Start reading data:" + dataset)
    savepath = datapath + "/TM5Inversion/dataframes_flux/"
    if 'RemoTeCISloc_Flux' in dataset:
        sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/RemoTeC+IS-land_ocean_bc/2009010100/2019070100/output/"        
    elif 'RemoTeCIS_Flux' in dataset:
        sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/RemoTeC+IS/2009010100/2019070100/output/"
    elif 'RemoTeC_Flux' in dataset:
        if res == 'glb6x4':
            sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/RemoTeC/2009010100/2019070100/output/"    
        if res == 'glb3x2':
            sourcepath = datapath + "/TM5Inversion/glb3x2_20220413/new_vpp_files/RemoTeC_2.4.0/2009010100/2019070100/output/"    
    elif 'RemoTeC238ISloc_Flux' in dataset and res == 'glb6x4':
        sourcepath = datapath + "/TM5InversionRT238/RemoTeC_2.3.8+IS-land_ocean_bc/2009010100/2019070100/output/"        
    elif 'RemoTeC238IS_Flux' in dataset and res == 'glb6x4':
        sourcepath = datapath + "/TM5InversionRT238/RemoTeC_2.3.8+IS/2009010100/2019070100/output/"
    elif 'RemoTeC238_Flux' in dataset and res == 'glb6x4':
        sourcepath = datapath + "/TM5InversionRT238/RemoTeC_2.3.8/2009010100/2019070100/output/"    
    elif 'RemoTeC238ISloc_Flux' in dataset and res == 'glb3x2':
        sourcepath = datapath + "/TM5Inversion/glb3x2_20220413/vpp_files/RemoTeC_2.3.8+IS-land_ocean_bc/2009010100/2019070100/output/"        
    elif 'RemoTeC238_Flux' in dataset and res == 'glb3x2':
        sourcepath = datapath + "/TM5Inversion/glb3x2_20220413/vpp_files/RemoTeC_2.3.8/2009010100/2019070100/output/"    
    elif 'ACOSIS_Flux' in dataset and res == 'glb6x4':
        sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/ACOS+IS/2009010100/2019070100/output/"
    elif 'ACOS_Flux' in dataset and res == 'glb6x4':
        sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/ACOS/2009010100/2019070100/output/"
    elif 'ACOSIS_Flux' in dataset and res == 'glb3x2':
        sourcepath = datapath + "/TM5Inversion/glb3x2_20220413/vpp_files/ACOS+IS/2009010100/2019070100/output/"
    elif 'ACOS_Flux' in dataset and res == 'glb3x2':
        sourcepath = datapath + "/TM5Inversion/glb3x2_20220413/vpp_files/ACOS/2009010100/2019070100/output/"
    elif 'IS_Flux' in dataset:
        sourcepath = datapath + "/TM5Inversion/"+res+"/ml60/tropo25/IS/2009010100/2019070100/output/"
    
    else:
        print('wrong dataset given')
        
    print(sourcepath)
    #loops through the posterior and apriori flux files (VPP_APOS and VPP_APRI)
    for num, filepath in enumerate(glob.glob(sourcepath+"VPP_A*.nc")):
        DS = xr.open_mfdataset(filepath,group = 'CO2',combine='by_coords',concat_dim='None',decode_times=False)

        #create Dataframe
        df2= CreateDF(DS, filepath)

        if num == 0:
            df = df2.copy()
        else:           
            df = pd.merge(df,df2,on=['Region','Date','Month','Year'],how="inner")

    
    print("finished reading data") 
    df.to_pickle(savepath+"DF2_"+dataset.split('_')[0]+res+".pkl")

        

#main
for nam in ['RemoTeCISloc_Flux','ACOSIS_Flux','IS_Flux']:
    CreateDataFrameTM5Flux(nam,'glb3x2')
