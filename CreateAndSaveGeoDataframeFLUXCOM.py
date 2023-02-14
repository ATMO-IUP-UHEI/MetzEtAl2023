#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of FLUXCOM files
#Author: E.-M. Schoemann
#Date 26.11.2020

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame, getReferenceDate, getReferencesDateDay 
import glob, os
import geopandas
import xarray as xr
import time
from shapely.geometry import Polygon
import math
from skimage.measure import block_reduce
from RegionParam import getRegion
import warnings
warnings.filterwarnings('ignore', '.*Mean of empty slice.*', )

datapath = "."
savepath = "."

def CreateDF(DS,Datatype):
    for d in range(len(DS[Datatype])):
        for nam in [Datatype,Datatype+'_mad']:#,'GPP_mad','GPP_n']:        
            dfyCO2 = pd.DataFrame(data=block_reduce(DS.variables[nam][d].values,block_size=(12,12),func=np.nanmean),index = list(np.array(range(895,-905,-10))/10), columns = list(np.array(range(-1795,1805,10))/10), dtype='float' )
            dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                      value_name =nam,var_name = 'Long')
            dfyCO2["Lat"] = dfyCO2['index']
            dfyCO2.insert(loc=1,column='Year1',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
            dfyCO2.insert(loc=1,column='Month1',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
            dfyCO2.insert(loc=1,column='Day1',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
            dfyCO2.insert(loc=1,column='Year2',value= np.ones(len(dfyCO2["Lat"]))*((datetime.date(int(str(DS.time.values[d])[0:4]),int(str(DS.time.values[d])[5:7]),int(str(DS.time.values[d])[8:10]))+datetime.timedelta(7,0)).year))
            dfyCO2.insert(loc=1,column='Month2',value= np.ones(len(dfyCO2["Lat"]))*((datetime.date(int(str(DS.time.values[d])[0:4]),int(str(DS.time.values[d])[5:7]),int(str(DS.time.values[d])[8:10]))+datetime.timedelta(7,0)).month))
            dfyCO2.insert(loc=1,column='Day2',value= np.ones(len(dfyCO2["Lat"]))*((datetime.date(int(str(DS.time.values[d])[0:4]),int(str(DS.time.values[d])[5:7]),int(str(DS.time.values[d])[8:10]))+datetime.timedelta(7,0)).day))
            
    

            if nam in Datatype:
                dfy = dfyCO2.copy()
            else:
                dfy = pd.merge(dfyCO2, dfy, on=['Lat','Long','Day2', 'Month2', 'Year2', 'Day1', 'Month1','Year1'])

        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
                     
    return df

def CreateDataFrameSIFFC(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    
    a_corr2 = []

    for i in range(895,-905,-10):
        geom = [Polygon(zip([100,100,101,101],[i/10-0.5,i/10+0.5,i/10+0.5,i/10-0.5]))]
        GData = geopandas.GeoDataFrame({'data':[1]}, geometry=geom)
        GData.crs = 'epsg:4326'
        GData = GData.to_crs({'proj':'cea'})
        a_corr2.append(GData.geometry.area[0])
    
    if not os.path.isfile(datapath + "/SIF/FLUXCOM/DF2_"+Datatype+RegionName+str(Num)+".pkl"):
        if not os.path.isfile(datapath + "/SIF/FLUXCOM/DFD02_"+Datatype+RegionName+str(Num)+".pkl"):
            print("Start reading data:")
            if Datatype == 'NEE':
                print('get files for NEE')
                dp = glob.glob(datapath + "/FLUXCOM/NEE.RS_V006.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.8daily.201*.nc")
                dp.extend(glob.glob(datapath + "/FLUXCOM/NEE.RS_V*2009.nc"))
                
            elif Datatype == 'TER':
                dp = glob.glob(datapath + "/FLUXCOM/"+Datatype+".RS_V006.FP-ALL.MLM-ALL.METEO-NONE.4320_2160.8daily.201*.nc")
                dp.extend(glob.glob(datapath + "/FLUXCOM/"+Datatype+".RS_V*2009.nc"))
            else:
                dp = glob.glob("/mnt/data/users/lschrein/FLUXCOM/"+Datatype+".RS_V006.FP-ALL.MLM-ALL.METEO-NONE.4320_2160.8daily.201*.nc")
                dp.extend(glob.glob("/mnt/data/users/lschrein/FLUXCOM/"+Datatype+".RS_V*2009.nc"))
            for num, filepath in enumerate(dp):
                print(filepath)
                DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_bnds')#decode_times=False)
       
                #create Dataframe
                df3= CreateDF(DS,Datatype)
                
                df2 = df3[(df3.Long >= Long_min)
                 & (df3.Long <= Long_max)
                 & (df3.Lat >= Lat_min)
                 & (df3.Lat <= Lat_max)]

                if num == 0:
                    df = df2.copy()
                else:           
                    df = df.append(df2, ignore_index=True)
        

        
            print("finished reading data") 
            df.to_pickle(datapath + "/SIF/FLUXCOM/DFD02_"+Datatype+RegionName+str(Num)+".pkl")
        else:
            df = pd.read_pickle(datapath + "/SIF/FLUXCOM/DFD02_"+Datatype+RegionName+str(Num)+".pkl")
        #create date variable
        print("create timestamp")
        
        dfd = getReferencesDateDay(2009,2019,1,12,1,31)
       
        gpp = []
        gpp_mad = []
        #gpp_n = []
        gppt = []
        gpp_madt = []
        #gpp_nt = []
        lat = []
        lon = []
        date = []
        ml = []
        dl = []
        yl = []
        # the reference date set has an entry for ever day, each day gets a NEE value (as FLUXCOM has 8 day means, 8 days will get the same value)
        for i,row in dfd.iterrows():
            for ind in df[(row.Year <= df.Year2)&(row.Year >= df.Year1)&(row.Month <= df.Month2)&(row.Month >= df.Month1)&(row.Day <= df.Day2+31*(df.Month2-row.Month))&(row.Day >= df.Day1-31*(row.Month-df.Month1))].index.values:

                gpp.append(df.iloc[ind][Datatype])
                gpp_mad.append(df.iloc[ind][Datatype+'_mad'])
                gppt.append((df.iloc[ind][Datatype]*a_corr2[int(90-(df.iloc[ind].Lat+0.5))]))
                gpp_madt.append((df.iloc[ind][Datatype+'_mad']*a_corr2[int(90-(df.iloc[ind].Lat+0.5))]))
                lat.append(df.iloc[ind].Lat)
                lon.append(df.iloc[ind].Long)            
                ml.append(int(row.Month))
                yl.append(int(row.Year))
                dl.append(int(row.Day))
                date.append(datetime.date(int(row.Year),int(row.Month),int(row.Day)))    
        del dfd
        datad = {'Date':date,Datatype:gpp,Datatype+'tot':gppt,'Lat':lat,'Long':lon,'Year':yl,'Month':ml,'Day':dl}
        dfd = pd.DataFrame(data=datad)
        dfd.insert(loc=1,column=Datatype+'_mad',value=gpp_mad)
        dfd.insert(loc=1,column=Datatype+'_madtot',value=gpp_madt)
        



        dfd.to_pickle(datapath + "/SIF/FLUXCOM/DF2_"+Datatype+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    else:
        dfd = pd.read_pickle(datapath + "/SIF/FLUXCOM/DF2_"+Datatype+RegionName+str(Num)+".pkl") 


    gdf = geopandas.GeoDataFrame(
        dfd, geometry=geopandas.points_from_xy(dfd.Long, dfd.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(datapath + "/SIF/FLUXCOM/GDF2_"+Datatype+RegionName+str(Num)+".pkl")

'''main'''
'''
Numm  = 949
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

CreateDataFrameSIFFC(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'NEE')
'''