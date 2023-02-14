#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of CT2019B CO2 flux files
#Author: E.-M. Schoemann
#Date 26.11.2020

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas
import xarray as xr
import time
import math
from skimage.measure import block_reduce

datapath = "."
savepath = "."

def CreateDF2(DS):
    Flux_ocn = []
    Flux_bio = []
    Flux_ff = []
    Flux_fire = []

    Flux_ocn_apri = []
    Flux_bio_apri = []

    region = []
    Month = []
    Year = []
    Date =[]
    SecMonth = []
    for i in range(len(DS.bio_flux_prior_tc23.values)):
        print(i)
        
        for j in range(23):
            Flux_ocn.append(DS.ocn_flux_opt_tc23.values[i][j])
            Flux_bio.append(DS.bio_flux_opt_tc23.values[i][j])
            Flux_ff.append(DS.ff_flux_imp_tc23.values[i][j])
            Flux_fire.append(DS.fire_flux_imp_tc23.values[i][j])
            Flux_ocn_apri.append(DS.ocn_flux_prior_tc23.values[i][j])
            Flux_bio_apri.append(DS.bio_flux_prior_tc23.values[i][j])
                          
            Date.append(datetime.date(int(math.floor(DS.decimal_date.values[i])),int(math.ceil(12*(DS.decimal_date.values[i]-math.floor(DS.decimal_date.values[i])))),15))
            region.append(j)
            Month.append(int(math.ceil(12*(DS.decimal_date.values[i]-math.floor(DS.decimal_date.values[i])))))
            Year.append(int(math.floor(DS.decimal_date.values[i])))
            SecMonth.append(DS.seconds_in_month.values[i])
            
    kind = ''            
    d = {'Flux_ocn' + kind: Flux_ocn, 
         'Flux_bio' + kind: Flux_bio, 
         'Flux_ff' + kind:Flux_ff,
         'Flux_fire' + kind:Flux_fire,
         'Flux_bio_apri' + kind:Flux_bio_apri,
         'Flux_ocn_apri' + kind:Flux_ocn_apri,
         'Region': region, 
         'Year': Year, 
         'Month': Month, 
         'SecMonth':SecMonth,
         'Date': Date}
    df = pd.DataFrame(data=d)

    return df 

def CreateDataFrameCTFlux2():
    filepath = datapath + "/CT2019/fluxMonthly/CT2019B.regionfluxes_monthly.nc"
    DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',decode_times=False)
    df= CreateDF2(DS)
    print("finished reading data") 

    #create date variable
    print("create timestamp")
    

    df.insert(loc=1, column = 'Biotot', value = df.Flux_bio*df.SecMonth*44.01)
    df.insert(loc=1, column = 'Ocntot', value = df.Flux_ocn*df.SecMonth*44.01)
    df.insert(loc=1, column = 'fftot', value = df.Flux_ff*df.SecMonth*44.01)
    df.insert(loc=1, column = 'firetot', value = df.Flux_fire*df.SecMonth*44.01)
    df.insert(loc=1, column = 'Biotot_apri', value = df.Flux_bio_apri*df.SecMonth*44.01)
    df.insert(loc=1, column = 'Ocntot_apri', value = df.Flux_ocn_apri*df.SecMonth*44.01)


    df.to_pickle(datapath + "/CT2019/dataframes/GDF2MonthlyRegionFLUX.pkl")
    #create GeoDataFrame



def CreateDF(DS2,DS,Datatype):
     
    dNBE = pd.DataFrame(data=DS2.variables['area'].values,index = [-90] + list(range(-88,90,4)),
                                  columns =  list(np.array(range(-1775,1825,50))/10), dtype='float' )
    dNBE = pd.melt(dNBE.reset_index(), id_vars='index',
                                 value_name ='area',var_name = 'Long')
    dNBE["Lat"] = dNBE['index']
    dNBE2 = pd.DataFrame(data=DS2.variables['LandMask'].values,index = [-90] + list(range(-88,90,4)),
                                  columns =  list(np.array(range(-1775,1825,50))/10), dtype='float' )
    dNBE2 = pd.melt(dNBE2.reset_index(), id_vars='index',
                                 value_name ='LandMask',var_name = 'Long')
    dNBE2["Lat"] = dNBE2['index']
    dNBE = pd.merge(dNBE,dNBE2,on=['Lat','Long'],how='left')
    dNBE = dNBE[(dNBE.Lat > -90)]
    #print(dfy2)
        
    for d in range(len(DS['bio_flux_opt'])):
        for nam in ['bio_flux_opt','fire_flux_imp','fossil_flux_imp','ocn_flux_opt']:        
            #reduce gridresolution to 4x5°
            dfyCO2 = pd.DataFrame(data=block_reduce(DS.variables[nam][d].values,block_size=(4,5),func=np.mean),index = list(range(-88,92,4)), columns =list(np.array(range(-1775,1825,50))/10) , dtype='float' )
            dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                      value_name =nam+'0',var_name = 'Long')
            dfyCO2["Lat"] = dfyCO2['index']
            dfyCO2 = pd.merge(dfyCO2, dNBE, on=['Lat','Long'], how = 'left')
            dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
            dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
            dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
            dfyCO2.insert(loc=1,column='Hour',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[11:13]))
            dfyCO2.insert(loc=1,column=nam,value=dfyCO2[nam+'0']*dfyCO2['area']*1000000*dfyCO2.LandMask) #area is in km**2 and fluxes in m**2

            if nam in 'bio_flux_opt':
                dfy = dfyCO2.copy()
            else:
                dfy = pd.merge(dfyCO2, dfy, on=['Lat','Long'])        
                dfy = dfy.drop(columns=['Hour_y', 'Month_y','Year_y','Day_y'])
                dfy = dfy.rename(columns={"Hour_x": "Hour", "Day_x":"Day","Month_x": "Month","Year_x":"Year"})
            del dfyCO2
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
                     
    return df

def CreateDataFrameCTFlux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    DS2 = xr.open_mfdataset(datapath + "/CMS-Flux NBE 2020/CMS-Flux-NBE-2020.monthly.grid.nc",combine='by_coords',concat_dim='None',decode_times=False)
    for y in range(2009,2020): 
        print("Start reading data for : " +str(y))
        dp = datapath + "/CT2019/flux3hour/CT2019B.flux1x1."+str(y)+"*.nc"
        for num, filepath in enumerate(glob.glob(dp)):
            #print(filepath)
            DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_components')#decode_times=False)
       
            #create Dataframe
            df3= CreateDF(DS2,DS,Datatype)
            del DS
            df2 = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)]

            if num == 0:
                df = df2.copy()
            else:           
                df = df.append(df2, ignore_index=True)
        
        #print(df)
        
        #print(df.keys())
        print(df.keys)
       
        print("finished reading data") 

        #create date variable
        print("create timestamp")
        date = []
        for i,row in df.iterrows():
            date.append(datetime.date(int(row.Year),int(row.Month),int(row.Day)))    
        df.insert(loc=1,column='Date',value=date)


        biot = []
        fft = []
        ocnt = []
        firet = []
        biots = []
        ffts = []
        ocnts = []
        firets = []
        for i,row in df.iterrows():
            #       sec/min * min/h * 3h * g/mol #* m²/grid 
            biot.append(row.bio_flux_opt*60*60*3*44.01)#*111319*111000*math.cos(math.radians(row.Lat))*44.01)
            fft.append(row.fossil_flux_imp*60*60*3*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            ocnt.append(row.ocn_flux_opt*60*60*3*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            firet.append(row.fire_flux_imp*60*60*3*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            #        m²/grid * g/mol
            biots.append(row.bio_flux_opt*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            ffts.append(row.fossil_flux_imp*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            ocnts.append(row.ocn_flux_opt*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)
            firets.append(row.fire_flux_imp*44.01)#111319*111000*math.cos(math.radians(row.Lat))*44.01)

        df.insert(loc=1, column = 'Biotot', value = biot)
        df.insert(loc=1, column = 'Icntot', value = ocnt)
        df.insert(loc=1, column = 'fftot', value = fft)
        df.insert(loc=1, column = 'firetot', value = firet)
        df.insert(loc=1, column = 'Biotot_s', value = biots)
        df.insert(loc=1, column = 'Icntot_s', value = ocnts)
        df.insert(loc=1, column = 'fftot_s', value = ffts)
        df.insert(loc=1, column = 'firetot_s', value = firets)


        df.to_pickle(datapath + "/CT2019/dataframes/DF2FLUX"+str(y)+"_"+RegionName+str(Num)+".pkl")
        #create GeoDataFrame


        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
        gdf.crs = {'init' :'epsg:4326'}
     
        if Num >= 900:
            if y == 2009:
                Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
            igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
            gdf = gdf.loc[igdf]
       
        gdf.to_pickle(datapath + "/CT2019/dataframes/GDF2FLUX"+str(y)+"_"+RegionName+str(Num)+".pkl")

        del gdf,df
        
def CombineGDF(RegionName, Num):
    for y in range(2009,2020):
        df2 = pd.read_pickle(datapath + "/CT2019/dataframes/GDF2FLUX"+str(y)+"_"+RegionName+str(Num)+".pkl")
        if y == 2009:
            df = df2.copy()
        else:           
            df = df.append(df2, ignore_index=True)
            
    df.to_pickle(datapath + "/CT2019/dataframes/GDF2FLUX_"+RegionName+str(Num)+".pkl")

