#!/usr/bin/env python
# Create a GeoDataFrame out of GFAS files
#Author: E.-M. Schoemann
#Date 26.11.2020

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas
import xarray as xr
import time
import math

datapath = "."
savepath = "."

def CreateDF(DS,Datatype):
    year = []
    month = []
    day = []
    date = []
    nam = 'co2fire'
    if 'Month' in Datatype:
        nam = 'fire'

    for d in range(len(DS[nam])):
        try:
            dfyCO2 = pd.DataFrame(data=DS.variables['co2fire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
        except:
            dfyCO2 = pd.DataFrame(data=DS.variables['fire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
        dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                  value_name ='CO2fire',var_name = 'Long')
        dfyCO2["Lat"] = dfyCO2['index']
        dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
        dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
        dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
        if 'Month' not in Datatype:
            dfyCO = pd.DataFrame(data=DS.variables['cofire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyCO = pd.melt(dfyCO.reset_index(), id_vars='index',
                                  value_name ='COfire',var_name = 'Long')
            dfyCO["Lat"] = dfyCO['index']

            dfyC = pd.DataFrame(data=DS.variables['cfire'][d].values,index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyC = pd.melt(dfyC.reset_index(), id_vars='index',
                                  value_name ='Cfire',var_name = 'Long')
            dfyC["Lat"] = dfyC['index']

            dfy = pd.merge(dfyCO2, dfyCO, on=['Lat','Long'])
            dfy = pd.merge(dfy, dfyC, on=['Lat','Long'])
        else:
            dfy = dfyCO2

        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
        
    return df

def CreateDataFrameGFAS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    if True: #not os.path.isfile(datapath + "/GFAS/dataframes/DF1_"+RegionName+str(Num)+".pkl"):
        print("Start reading data:")
        if 'Month' in Datatype:
            dp = datapath + "/GFAS/cams_gfas_CO2_2009_2020_daily_1x1_monthMeans.nc"
        else:
            dp = datapath + "/GFAS/Daily_1x1_regridded/*.nc"
        for num, filepath in enumerate(glob.glob(dp)):
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
       
            #create Dataframe
            df3= CreateDF(DS,Datatype)
            if Long_min >= 0:
                df2 = df3[(df3.Long >= Long_min)
                 & (df3.Long <= Long_max)
                 & (df3.Lat >= Lat_min)
                 & (df3.Lat <= Lat_max)]
            elif Long_max < 0:
                df2 = df3[(df3.Long >=(360 + Long_min))
                    & (df3.Long <= (360 + Long_max))
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]
            else:
                df2 = df3[(((df3.Long >=(360 + Long_min)) & (df3.Long <= 360))|
                      ((df3.Long >=(0)) & (df3.Long <= (Long_max))))
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]

            if num == 0:
                df = df2.copy()
            else:           
                df = df.append(df2, ignore_index=True)
        
       
        print("finished reading data") 

        #create date variable
        print("create timestamp")
 
        date = []
        for i,row in df.iterrows():
            date.append(datetime.date(int(row.Year),int(row.Month),int(row.Day)))    
        df.insert(loc=1,column='Date',value=date)

        long2 = []
        Cf = []
        COf = []
        CO2f = []
        for i,row in df.iterrows():
            #        (km*m**-2*s**-1*s/month*m**2/cell*g/kg)
            if 'Month' not in Datatype:
                Cf.append(row.Cfire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
                COf.append(row.COfire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
            CO2f.append(row.CO2fire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
            if row.Long > 180:
                long2.append(row.Long - 360)
            else:
                long2.append(row.Long)

        df = df.drop(columns = 'Long')
        df.insert(loc=1, column='Long',value = long2)
        if 'Month' not in Datatype:
            df.insert(loc=1, column = 'CfireE', value = Cf)
            df.insert(loc=1, column = 'COfireE', value = COf)
        df.insert(loc=1, column = 'CO2fireE', value = CO2f)


        df.to_pickle(datapath + "/GFAS/dataframes/DF1_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    else:
        df = pd.read_pickle(datapath + "/GFAS/dataframes/DF1_"+RegionName+str(Num)+".pkl") 

        Cf = []
        COf = []
        CO2f = []
        for i,row in df.iterrows():
            #        (km*m**-2*s**-1*s/month*m**2/cell*g/kg)
            if 'Month' not in Datatype:
                Cf.append(row.Cfire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
                COf.append(row.COfire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)
            CO2f.append(row.CO2fire*86400*111319*111000*math.cos(math.radians(row.Lat))*1000)

        if 'Month' not in Datatype:
            df.insert(loc=1, column = 'CfireE', value = Cf)
            df.insert(loc=1, column = 'COfireE', value = COf)
        df.insert(loc=1, column = 'CO2fireE', value = CO2f)

    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(datapath + "/GFAS/dataframes/GDF1_"+RegionName+str(Num)+".pkl")


