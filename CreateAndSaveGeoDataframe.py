#!/usr/bin/env python
# Create a GeoDataFrame out of  GOSAT .out files
#Author: E.-M. Schoemann
#Date 10.06.2020

import datetime
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas

datapath = "."

os.chdir(datapath + "/GOSAT_Markus/")
print("Start reading data:")
for num, file in enumerate(glob.glob("*.out")):
    print("read: "+file)
    data = read_remotec_out.read_output_file(datapath + "/GOSAT_Markus/" + file)
    if num == 0:        
        #create Dataframe
        df= createPandasDateFrame(data)
        
    else:      
        #create Dataframe
        df2= createPandasDateFrame(data)         
        df = df.append(df2, ignore_index=True)
        
print("finished reading data") 

#create date variable
print("create timestamp")
date = []
for i in range(len(df.Year)):
    date.append(datetime.date(df.Year[i],df.Month[i],df.Day[i]))    
df["Date"] = date

#create calender week number variable
print("create calender weeknumber")
week = []
daysMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
daysMonth_2 = [31,29,31,30,31,30,31,31,30,31,30,31]
for j in range(len(df.Year)):
    if df.Year[j] == 2012 or df.Year[j] == 2016:
        week.append((sum(daysMonth_2[0:df.Month[j]-1])+df.Day[j])/7)    
    else:
        week.append((sum(daysMonth[0:df.Month[j]-1])+df.Day[j])/7)  
week = np.floor(np.array(week))
df["Week"] = week


df.to_pickle(datapath + "/GOSAT_Markus/dataframes/DF09_19_8.pkl")
df.to_csv(datapath + "/GOSAT_Markus/dataframes/DF09_19_8.csv")
#create GeoDataFrame
gdf = geopandas.GeoDataFrame(
    df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
gdf.crs = {'init' :'epsg:4326'}
#Load world Data
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
countries = world[['geometry', 'name','continent']]

gdf2 = geopandas.sjoin(gdf, countries, how="left", op='intersects')
gdf2["transcom"] = np.nan


# add transcom region !NEED A CHANGE WITH BETTER SHAPEFILES
gdf2.loc[(gdf2.continent == 'North America') & (gdf2.Lat >= 50), "transcom"] = 1
gdf2.loc[(gdf2.continent == 'North America') & (gdf2.Lat < 50), "transcom"] = 2
gdf2.loc[(gdf2.continent == 'South America') & (gdf2.Lat >= -10), "transcom"] = 3
gdf2.loc[(gdf2.continent == 'South America') & (gdf2.Lat < -10), "transcom"] = 4
gdf2.loc[(gdf2.continent == 'Africa') & (gdf2.Lat >= 0), "transcom"] = 5
gdf2.loc[(gdf2.continent == 'Africa') & (gdf2.Lat < 0), "transcom"] = 6
gdf2.loc[(gdf2.continent == 'Europe') & (gdf2.Long > 60), "transcom"] = 7
gdf2.loc[(gdf2.continent == 'Asia')&~((gdf2.name == 'Thailand') |(gdf2.name == 'Myanmar') |(gdf2.name == 'Vietnam')|(gdf2.name ==  'Malaysia')|(gdf2.name == 'Indonesia')|(gdf2.name == 'Papua New Guinea')|(gdf2.name == 'Laos')|(gdf2.name == 'Cambodia')|(gdf2.name == 'Philippines')|(gdf2.name == 'Taiwan')|(gdf2.name == 'Brunei')), "transcom"] = 8
gdf2.loc[(gdf2.continent == 'Asia')&((gdf2.name == 'Thailand') |(gdf2.name == 'Myanmar') |(gdf2.name == 'Vietnam')|(gdf2.name ==  'Malaysia')|(gdf2.name == 'Indonesia')|(gdf2.name == 'Papua New Guinea')|(gdf2.name == 'Laos')|(gdf2.name == 'Cambodia')|(gdf2.name == 'Philippines')|(gdf2.name == 'Taiwan')|(gdf2.name == 'Brunei')), "transcom"] = 9
gdf2.loc[(gdf2.continent == 'Oceania'), "transcom"] = 10
gdf2.loc[(gdf2.continent == 'Europe') & (gdf2.Long < 60)& (gdf2.Long > -40), "transcom"] = 11

Dat = []
for i in range(len(gdf2.Year)):
    Dat.append(gdf2.Year[i]*100 +gdf2.Week[i])
gdf2["Week_date"] = Dat

gdf2.to_pickle(datapath + "/GOSAT_Markus/dataframes/GDF10_19_8.pkl")
gdf2.to_csv(datapath + "/GOSAT_Markus/dataframes/GDF10_19_8.csv")


