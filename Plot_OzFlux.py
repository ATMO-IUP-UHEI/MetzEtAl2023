#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Created on Fri Okt 22

@author: eschoema
Script to investigate OzFlux Carbon Fluxes

"""

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import getDayMeans, getDaySum, getReferenceDateDay, getMonthMeans, getReferenceDate, getMonthSum
import geopandas
import xarray as xr
import argparse
import os.path
from RegionParam import getRegion
#from CreateAndSaveGeoDataframeCTFlux import CreateDataFrameCTFlux, CombineGDF,CreateDataFrameCTFlux2
#from CreateAndSaveGeoDataframeFLUXCOM import CreateDataFrameSIFFC
#from CreateAndSaveGeoDataframeCAMSFlux import CreateDataFrameCAMSflux, CombineGDFCAMS
#from CreateAndSaveGeoDataframeCABLE import CreateDataFrameCABLE
from sklearn.linear_model import LinearRegression
import matplotlib as mpl

datapath = "."
savepath = "."

font = {'family' : 'Arial'}
mpl.rc('font', **font)

Numm = 949

startdate = '2009-04-15'#'2015-01-01'#'2009-04-15'#'2014-08-15'#
enddate = '2018-12-31'

year_min = int(startdate[0:4])#2009
month_min = int(startdate[5:7])#4
day_min = int(startdate[8:10])#1
year_max = int(enddate[0:4])#2019
month_max = int(enddate[5:7])#12
day_max = int(enddate[8:10])#31

RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
DateRef = getReferenceDateDay(year_min,year_max,month_min,month_max)
DateRefMonth = getReferenceDate(year_min,year_max,month_min,month_max)


Stations = ['ADdry','ADwet','ASM','BeaconFarm_IFR',
                'BeaconFarm_UUW','Boyagin1','Calperum','CapeTribulation',
                'CowBay','CumberlandPlain','CumberlandMelaleuca',
                'DalyUncleared','Dargo','Digby','DryRiver','DalyPasture',
                'Emerald','GatumPasture','Gingin','GWW','HowardSprings',
                'HowardUnderstory','Kopuatai','Litchfield','Longreach','Nimmo',
                'Otway','RDMF','Ridgefield','Riggs','Robson','Samford',
                'SturtPlains','TTE','Tumbarumba','Wallaby','Warra',
                'Whroo','WombatStateForest','Yanco']

dfswscharac = pd.DataFrame(data = {'Station':['test'],
                            'Max':0,
                            'Min':0,
                            'limit':0})
        
first = True

for stat in range(len(Stations)):
    onestation = False
    station_id = Stations[stat]
    print(station_id)
    dfStation = pd.read_pickle(datapath + "/OzFlux/DataFrames/DF5_"+station_id+".pkl")
    
    VegVar = 'Fc'
    if dfStation.Precipitation.min() == -9999:
        dfStation.loc[dfStation[(dfStation.Precipitation == -9999)].index,
                                'Precipitation'] = np.nan
    if dfStation[VegVar].min() == -9999:
        dfStation.loc[dfStation[(dfStation[VegVar] == -9999)].index,
                                VegVar] = np.nan
    if dfStation.Sws.min() == -9999:
        dfStation.loc[dfStation[(dfStation.Sws == -9999)].index,
                                'Sws'] = np.nan
    try:
        if dfStation.CO2.min() == -9999:
            dfStation.loc[dfStation[(dfStation.CO2 == -9999)].index,
                                    'CO2'] = np.nan
    except:
        pass
    try:
        if dfStation.ER_night.min() == -9999:
            dfStation.loc[dfStation[(dfStation.ER_night == -9999)].index,
                                    'ER_night'] = np.nan
    except:
        pass
    if dfStation.Ts.min() < -100:
            dfStation.loc[dfStation[(dfStation.Ts < -100)].index,
                                    'Ts'] = np.nan
   
    try:
        if dfStation.Ta.min() < -100:
            dfStation.loc[dfStation[(dfStation.Ta < -100)].index,
                                    'Ta'] = np.nan
    except:
        pass
    
    if not dfStation.Lat.any():
        dfStation.drop(columns = ['Lat','Long'],inplace = True)
        dfStation.insert(loc = 1, column = 'Lat', value = -30*np.ones(len(dfStation.Fc)))
        dfStation.insert(loc = 1, column = 'Long', value = 130*np.ones(len(dfStation.Fc)))
    
    MonthMeansNEE = getMonthMeans(dfStation,
                     year_min,month_min,day_min,
                     year_max,month_max,day_max,
                     Long_min,Long_max,
                     Lat_min,Lat_max,
                     VegVar)
    MonthSumPrecip = getMonthSum(dfStation,
                     year_min,month_min,day_min,
                     year_max,month_max,day_max,
                     Long_min,Long_max,
                     Lat_min,Lat_max,
                     'Precipitation')
    MonthdateNEE = MonthMeansNEE.apply(lambda x: datetime.date(int(x.Year),
                                                      int(x.Month),
                                                      15),axis=1)
    MonthMeansNEE.insert(loc = 1, value = MonthdateNEE, column = 'MonthDate')
    MonthdatePrecip = MonthSumPrecip.apply(lambda x: datetime.date(int(x.Year),
                                                      int(x.Month),
                                                      15),axis=1)
    MonthSumPrecip.insert(loc = 1, value = MonthdatePrecip, column = 'MonthDate')
    DayMeansNEE = getDayMeans(dfStation,
                     year_min,month_min,day_min,
                     year_max,month_max,day_max,
                     Long_min,Long_max,
                     Lat_min,Lat_max,
                     VegVar)
    DayMeansSws = getDayMeans(dfStation,
                     year_min,month_min,day_min,
                     year_max,month_max,day_max,
                     Long_min,Long_max,
                     Lat_min,Lat_max,
                     'Sws')
    DayMeansPrecip = getDaySum(dfStation,
                     year_min,month_min,day_min,
                     year_max,month_max,day_max,
                     Long_min,Long_max,
                     Lat_min,Lat_max,
                     'Precipitation')
    
    #set sws lmiit
    Swsthr = DayMeansSws.Sws.min() + 0.1*(DayMeansSws.Sws.max()-DayMeansSws.Sws.min())
    if Swsthr > 0.07:
        Swsthr = 0.07
    
    
    if not onestation:
        dfstat = pd.DataFrame(data = {'Station':[station_id],
                            'Max':DayMeansSws.Sws.max(),
                            'Min':DayMeansSws.Sws.min(),
                            'Mean':DayMeansSws.Sws.mean(),
                            'limit':Swsthr})
        dfswscharac = dfswscharac.append(dfstat,ignore_index=True)
        
    try:
        DayMeansCO2 = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'CO2')
        PDayMeansCO2 = pd.merge(DateRef, DayMeansCO2, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    try:
        DayMeansER_night = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'ER_night')
        PDayMeansER_night = pd.merge(DateRef, DayMeansER_night, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    try:
        DayMeansNEE2 = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'NEE')
        PDayMeansNEE2 = pd.merge(DateRef, DayMeansNEE2, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    try:
        DayMeansNEP = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'NEP')
        PDayMeansNEP = pd.merge(DateRef, DayMeansNEP, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    try:
        DayMeansTs = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'Ts')
        PDayMeansTs = pd.merge(DateRef, DayMeansTs, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    try:
        DayMeansTa = getDayMeans(dfStation,
                         year_min,month_min,day_min,
                         year_max,month_max,day_max,
                         Long_min,Long_max,
                         Lat_min,Lat_max,
                         'Ta')
        PDayMeansTa = pd.merge(DateRef, DayMeansTa, on=['Year','Month','Day','Date'], how = 'left')
    except:
        pass
    PMonthMeansNEE = pd.merge(DateRefMonth, MonthMeansNEE, on=['MonthDate'], how = 'left')
    PMonthSumPrecip = pd.merge(DateRefMonth, MonthSumPrecip, on=['MonthDate'], how = 'left')
    PDayMeansNEE = pd.merge(DateRef, DayMeansNEE, on=['Year','Month','Day','Date'], how = 'left')
    PDayMeansSws = pd.merge(DateRef, DayMeansSws, on=['Year','Month','Day','Date'], how = 'left')
    PDayMeansPrecip = pd.merge(DateRef, DayMeansPrecip, on=['Year','Month','Day','Date'], how = 'left')
    
     
    colorl = ['yellow','orange','red','pink',
              'violet','purple','darkblue','lightblue',
              'teal','green','darkgreen','brown']
    
    
    more4mm = PDayMeansPrecip.Precipitation >= 5
    SubPrec = PDayMeansPrecip[more4mm]
    PDayMeansNEE.insert(loc = 1,column = 'd'+VegVar,value =PDayMeansNEE[VegVar]-([0]+list(PDayMeansNEE[VegVar].values[:-1])))
    d2daysNEE = [0,0]
    for num in range(2,len(PDayMeansNEE)):
        d2daysNEE.append(PDayMeansNEE[VegVar].values[num:num+2].mean()-PDayMeansNEE[VegVar].values[num-2:num].mean())
    PDayMeansNEE.insert(loc = 1,column = 'd2'+VegVar,value =d2daysNEE)
    
    SubNEE = PDayMeansNEE[more4mm]
    
    #get all days with soilsmoisture beeing lower than threshold within the last 5 or or 1 days
    drybefore = []
    periodlength = 1
    swsthreshold = Swsthr 
    for i in range(len(PDayMeansSws.Sws)):
        if i >= periodlength:
            start = i-periodlength
        else:
            start = 0
        if PDayMeansSws.Sws[start:i].min() <= swsthreshold:
            drybefore.append(True)
        else:
            drybefore.append(False)
    
    
           
    PDayMeansNEE.insert(loc=1, column=VegVar+'nextday',value = np.concatenate((PDayMeansNEE[VegVar].values[1:],np.array([0]))))
    PDayMeansNEE.insert(loc=1, column=VegVar+'7daybefore',value = np.concatenate((np.array([0,0,0,0,0,0,0]),PDayMeansNEE[VegVar].values[:-7])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'6daybefore',value = np.concatenate((np.array([0,0,0,0,0,0]),PDayMeansNEE[VegVar].values[:-6])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'5daybefore',value = np.concatenate((np.array([0,0,0,0,0]),PDayMeansNEE[VegVar].values[:-5])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'4daybefore',value = np.concatenate((np.array([0,0,0,0]),PDayMeansNEE[VegVar].values[:-4])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'3daybefore',value = np.concatenate((np.array([0,0,0]),PDayMeansNEE[VegVar].values[:-3])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'2daybefore',value = np.concatenate((np.array([0,0]),PDayMeansNEE[VegVar].values[:-2])))
    PDayMeansNEE.insert(loc=1, column=VegVar+'1daybefore',value = np.concatenate((np.array([0]),PDayMeansNEE[VegVar].values[:-1])))
    
    PDayMeansPrecipdry = PDayMeansPrecip.loc[drybefore]
    PDayMeansNEEdry = PDayMeansNEE.loc[drybefore]
    PDayMeansdry = pd.merge(PDayMeansPrecipdry,PDayMeansNEEdry, on = ['Year','Month','Day','Date'],how = 'inner')
    
    if not any(drybefore):    
        PDayMeansPrecipBirch = pd.DataFrame(data ={'Year':[], 'FcMinBefore':[], 'FcMax':[], 'Month':[], 'Day':[], 'Date':[], 'Precipitation':[],
           'Fc1daybefore':[], 'Fc2daybefore':[], 'Fcnextday':[], 'd2Fc':[], 'dFc':[], 'Fc':[]})
    else:
        PDayMeansdry.insert(loc = 1, column = VegVar+'Max', value = PDayMeansdry.apply(lambda x: np.nanmax(np.array([x[VegVar], x[VegVar+'nextday']])),axis = 1))
        PDayMeansdry.insert(loc = 1, column = VegVar+'MinBefore', value = PDayMeansdry.apply(lambda x: np.nanmin(np.array([x[VegVar+'1daybefore'], x[VegVar+'2daybefore']])),axis = 1))
        PDayMeansdry.insert(loc = 1, column = VegVar+'Max7Before', value = PDayMeansdry.apply(lambda x: np.nanmax(np.array([x[VegVar+'1daybefore'], x[VegVar+'2daybefore'],x[VegVar+'3daybefore'], x[VegVar+'4daybefore'],x[VegVar+'5daybefore'], x[VegVar+'6daybefore'],x[VegVar+'7daybefore']])),axis = 1))
        try:
            PDayMeansPrecipBirch = PDayMeansdry[(PDayMeansdry.Precipitation > 5)
                                    & (PDayMeansdry[VegVar+'Max'] > 0)
                                    & (PDayMeansdry[VegVar+'Max7Before'] < PDayMeansdry[VegVar+'Max'])]
        except:
            PDayMeansPrecipBirch = PDayMeansdry[(PDayMeansdry.Precipitation <-100)] #empty dataframe
        
    
    NumBirch = PDayMeansPrecipBirch.groupby(['Year'])['Month'].count().reset_index()
    NumBirch.rename(columns={'Month':station_id},inplace=True)
    NumBirchMonth = PDayMeansPrecipBirch.groupby(['Year','Month'])['Day'].count().reset_index()
    NumBirchMonth.rename(columns={'Day':station_id},inplace=True)
    DateRefWeeks = DateRef.copy()
    DateRefWeeks.insert(loc = 1, column = 'QuaterMonth',value = DateRefWeeks.apply(lambda x: 1 + x.Day//8 + 4*(x.Month-1) + 48*(x.Year - 2009),axis=1))
    PDayMeansPrecipBirchWeek = pd.merge(PDayMeansPrecipBirch,DateRefWeeks,on=['Year','Month','Day','Date'])
    NumBirchWeek = PDayMeansPrecipBirchWeek.groupby('QuaterMonth')['Day'].count().reset_index()
    NumBirchWeek.rename(columns={'Day':station_id},inplace=True)
    
    PDayMeansNEE.insert(loc = 1, column = 'notnull',value=PDayMeansNEE[VegVar].notnull()&PDayMeansPrecip['Precipitation'].notnull()&PDayMeansSws['Sws'].notnull() )
    PDayMeansNEE.insert(loc = 1, column = 'isnan',value=np.invert(PDayMeansNEE['notnull']))
    NumNan = PDayMeansNEE.groupby(['Year'])['isnan'].sum().reset_index()
    NumNanMonth = PDayMeansNEE.groupby(['Year','Month'])['isnan'].sum().reset_index()
    PDayMeansNEEweek = pd.merge(PDayMeansNEE,DateRefWeeks,on=['Year','Month','Day','Date'])
    NumNanWeek = PDayMeansNEEweek.groupby(['QuaterMonth'])['isnan'].sum().reset_index()
    
    DateRefWeeksDate = DateRefWeeks.groupby('QuaterMonth')['Date'].min().reset_index()
    DateRefWeeksDate.rename(columns={'Date':'WeekDate'},inplace=True)
    DateRefWeeks  = pd.merge(DateRefWeeks,DateRefWeeksDate,on=['QuaterMonth'],how='left')
    
    BirchRelease = PDayMeansPrecipBirch.groupby(['Year'])['FcMax'].sum().reset_index()
    BirchRelease.rename(columns={'FcMax':station_id},inplace=True)
    if stat == 0 or onestation:
        print('Yes')
        Years = pd.DataFrame(data={'Year':range(2009,2019)})
        dfBirchYears = pd.merge(Years,NumBirch, on = ['Year'], how = 'outer')
        dfBirchMonth = pd.merge(DateRefMonth,NumBirchMonth, on = ['Year','Month'], how = 'outer')
        dfBirchWeek = pd.merge(DateRefWeeksDate,NumBirchWeek, on = ['QuaterMonth'], how = 'outer')
        dfNansYears = pd.merge(Years,NumNan, on = ['Year'], how = 'outer')
        dfNansMonth = pd.merge(DateRefMonth,NumNanMonth, on = ['Year','Month'], how = 'outer')
        dfNansWeek = pd.merge(DateRefWeeksDate,NumNanWeek, on = ['QuaterMonth'], how = 'outer')
        
        dfBirchRelease = pd.merge(Years,BirchRelease, on = ['Year'], how = 'outer')
        
    else:
        dfBirchYears = pd.merge(dfBirchYears,NumBirch, on = ['Year'], how = 'outer')
        dfBirchMonth = pd.merge(dfBirchMonth,NumBirchMonth, on = ['Year','Month'], how = 'outer')
        dfBirchWeek = pd.merge(dfBirchWeek,NumBirchWeek, on = ['QuaterMonth'], how = 'outer')
        dfNansYears = pd.merge(dfNansYears,NumNan, on = ['Year'], how = 'outer')
        dfNansMonth = pd.merge(dfNansMonth,NumNanMonth, on = ['Year','Month'], how = 'outer')
        dfNansWeek = pd.merge(dfNansWeek,NumNanWeek, on = ['QuaterMonth'], how = 'outer')
        
        dfBirchRelease = pd.merge(dfBirchRelease,BirchRelease, on = ['Year'], how = 'outer')
    
    #Figure S 11
    if station_id == 'DryRiver' or station_id == 'ASM' or station_id == 'DalyUncleared':
        
        if station_id == 'DryRiver':
            Letter = 'A'
            rangeyear = 2014
            axisnum = 0
        elif station_id == 'ASM':
            Letter = 'C'
            rangeyear = 2015
            axisnum = 2
            
        elif station_id == 'DalyUncleared':
            Letter = 'B'
            rangeyear = 2018
            axisnum = 1
        if first:
            cm = 1/2.54  # centimeters in inches
            lw = 0.8
            cs = 1.3
            mpl.rcParams['hatch.linewidth'] = 0.5
            fig2, ax2 = plt.subplots(3, 1, 
                            figsize = (14.4*cm,22*cm),
                            gridspec_kw={'width_ratios': [1],'height_ratios':[1,1,1]})
    
            fig2.subplots_adjust(right=0.75)
            plt.xticks(fontsize=6)
            plt.yticks(fontsize=6)
            first = False
            
        for y in range(rangeyear,rangeyear+1):
            print(station_id)
            PDayMeansNEE = pd.merge(DateRef[(DateRef.Year == y)], DayMeansNEE, on=['Date'], how = 'left')
            PDayMeansSws = pd.merge(DateRef[(DateRef.Year == y)], DayMeansSws, on=['Date'], how = 'left')
            PDayMeansPrecip = pd.merge(DateRef[(DateRef.Year == y)], DayMeansPrecip, on=['Date'], how = 'left')
            PDayMeansCO2 = pd.merge(DateRef[(DateRef.Year == y)], DayMeansCO2, on=['Date'], how = 'left')
            PDayMeansNEE2 = pd.merge(DateRef[(DateRef.Year == y)], DayMeansNEE2, on=['Date'], how = 'left')
            PDayMeansNEP = pd.merge(DateRef[(DateRef.Year == y)], DayMeansNEP, on=['Date'], how = 'left')
            PDayMeansER_night = pd.merge(DateRef[(DateRef.Year == y)], DayMeansER_night, on=['Date'], how = 'left')
            
            PDayMeansPrecipBirchy = PDayMeansPrecipBirch[(PDayMeansPrecipBirch.Year == y)]
            PDayMeansPrecipdryy = PDayMeansPrecipdry[(PDayMeansPrecipdry.Year == y)]
            
            
        
            secAx = ax2[axisnum].twinx()
            secAx2 = ax2[axisnum].twinx()
            secAx2.spines["right"].set_position(("axes", 1.12))
            secAx2.set_frame_on(True)
            secAx2.patch.set_visible(False)
            for sp in secAx2.spines.values():
                sp.set_visible(False)
            secAx2.spines["right"].set_visible(True)
            
            ax2[axisnum].plot([PDayMeansNEE.Date[0],PDayMeansNEE.Date[len(PDayMeansNEE.Date)-1]],[0,0],linewidth = lw,ls='-',marker= '',color= 'grey')
            secAx.plot([PDayMeansNEE.Date[0],PDayMeansNEE.Date[len(PDayMeansNEE.Date)-1]],[swsthreshold,swsthreshold],linewidth = lw,ls='-',marker= '',color= 'darkred')
            ax2[axisnum].plot(PDayMeansNEE.Date,np.array(PDayMeansNEE[VegVar]),linewidth = lw,ls='-',marker= '',color= 'green',label=r'$\rm CO_2~flux$')
            secAx.plot(PDayMeansSws.Date,PDayMeansSws.Sws,linewidth = lw,ls=':',marker= '',color= 'brown',label='Sws')
            secAx2.bar(PDayMeansPrecip.Date,PDayMeansPrecip.Precipitation,width=datetime.timedelta(days=1.1),color= 'blue',label='Precipitation')
            secAx2.plot(PDayMeansPrecipBirchy.Date,PDayMeansPrecipBirchy .Precipitation, marker='.',markersize = 2, ls='', color = 'red')
            
            ax2[axisnum].legend(fontsize=6,loc=(0.4,0.89))#,loc=2)
            secAx.legend(fontsize=6,loc=(0.6,0.9))#,loc=9)
            secAx2.legend(fontsize=6,loc=(0.75,0.9))#,loc=1)
            ax2[axisnum].set_title(station_id +", "+str(y),fontsize=7)
            ax2[axisnum].grid(True,which = 'both', axis='x')
            ax2[axisnum].set_xlabel(r'Date',fontsize=7)
            ax2[axisnum].set_ylabel(r'$\rm CO_2~flux~[umol/(m^2s)]$',fontsize=6)
            ax2[axisnum].yaxis.label.set_color('green')
            ax2[axisnum].yaxis.label.set_fontsize(6)
            ax2[axisnum].tick_params(axis='y', colors='green',labelsize = 6)
            secAx.set_ylabel(r'Soil Water Fraction (Sws)',fontsize=6)
            secAx.yaxis.label.set_color('brown')
            secAx.yaxis.label.set_fontsize(7)
            secAx.tick_params(axis='y', colors='brown',labelsize = 6)
            secAx2.set_ylabel(r'Precipitation [mm/d]',fontsize=6)
            secAx2.yaxis.label.set_color('blue')
            secAx2.tick_params(axis='y', colors='blue',labelsize = 6)
            
            ax2[axisnum].text(0.03, 0.92, Letter, horizontalalignment='center', verticalalignment='center', transform=ax2[axisnum].transAxes,fontsize=10,weight='bold')
            ax2[axisnum].tick_params(axis = 'x',labelsize = 6)
plt.subplots_adjust(wspace=0,  
                    hspace=0.4) 
plt.savefig(savepath + "/Results/Plots/Paper2021/final/OzFLuxTimeSeries.png", dpi=300, bbox_inches = "tight")
plt.savefig(savepath + "/Results/Plots/Paper2021/final/OzFLuxTimeSeries.pdf", dpi=300, bbox_inches = "tight",format = 'pdf')
    
    
            
#Number of Birch Events
dfBirchWeek.insert(loc = 1, column = 'NumStationBirch', value = len(Stations) - dfBirchWeek.isnull().sum(axis=1))
dfBirchMonth.insert(loc = 1, column = 'NumStationBirch', value = len(Stations) - dfBirchMonth.isnull().sum(axis=1))
dfBirchYears.insert(loc = 1, column = 'NumStationBirch', value = len(Stations) - dfBirchYears.isnull().sum(axis=1))

dfBirchWeek.insert(loc = 1, column = 'NumBirch', value = dfBirchWeek[Stations].sum(axis=1))
dfBirchMonth.insert(loc = 1, column = 'NumBirch', value = dfBirchMonth[Stations].sum(axis=1))
dfBirchYears.insert(loc = 1, column = 'NumBirch', value = dfBirchYears[Stations].sum(axis=1))


dfNansWeek0 = dfNansWeek[['isnan_x','isnan_y']] < 4
dfNansMonth0 = dfNansMonth[['isnan_x','isnan_y']] < 15
dfNansYear0 = dfNansYears[['isnan_x','isnan_y']] < 150

dfNansWeek1 = dfNansWeek[['isnan_x','isnan_y']] == 0
dfNansMonth1 = dfNansMonth[['isnan_x','isnan_y']] == 0
dfNansYear1 = dfNansYears[['isnan_x','isnan_y']] == 0


dfNansWeek.insert(loc = 1, column = 'NumStationMeas', value = dfNansWeek0.sum(axis=1))
dfNansMonth.insert(loc = 1, column = 'NumStationMeas', value = dfNansMonth0.sum(axis=1))
dfNansYears.insert(loc = 1, column = 'NumStationMeas', value = dfNansYear0.sum(axis=1))

dfNansWeek.insert(loc = 1, column = 'NumStationAllMeas', value = dfNansWeek1.sum(axis=1))
dfNansMonth.insert(loc = 1, column = 'NumStationAllMeas', value = dfNansMonth1.sum(axis=1))
dfNansYears.insert(loc = 1, column = 'NumStationAllMeas', value = dfNansYear1.sum(axis=1))

#Birch Release
dfBirchRelease2 = dfBirchRelease.transpose()
dfBirchRelease2.rename(columns = dfBirchRelease.Year,inplace = True)
dfBirchRelease2.drop(dfBirchRelease2.index[[0]],inplace=True) 
meanrelease = dfBirchRelease2.mean(axis=1)
sumrelease = dfBirchRelease2.sum(axis=1)
dfBirchRelease2.insert(loc=0, column = 'Mean',value = meanrelease) 
dfBirchRelease2.insert(loc=0, column = 'Sum',value = sumrelease) 
dfBirchRelease2.to_pickle(savepath + "/Results/Plots/OzFlux/BirchMarker/20230120_1dsws/dfBirchRelease2.pkl")


# Figure S12
cm = 1/2.54  # centimeters in inches
lw = 0.8
cs = 1.3
fig4, ax4 = plt.subplots(figsize = (12.0*cm,6*cm))

ax4.bar(dfNansMonth.MonthDate,dfNansMonth.NumStationMeas,width=31, color = 'black',label = 'Measuring')
ax4.bar(dfNansMonth.MonthDate,dfNansMonth.NumStationAllMeas,width=31, color = 'grey',label = 'Gapless measuring')
ax4.bar(dfBirchMonth.MonthDate,dfBirchMonth.NumStationBirch,width=31,label = 'At least one pulse')

ax4.set_xlabel(r'Date',fontsize=7)
ax4.set_ylabel(r'# OzFlux stations',fontsize=7)
ax4.legend(fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.savefig(savepath + "/Results/Plots/Paper2021/final//BarPlotsBirchEventsMonthsWithAll_Partlymeas.png", dpi=400, bbox_inches = "tight")
plt.savefig(savepath + "/Results/Plots/Paper2021/final//BarPlotsBirchEventsMonthsWithAll_Partlymeas.pdf", dpi=400, bbox_inches = "tight",format = 'pdf')

     