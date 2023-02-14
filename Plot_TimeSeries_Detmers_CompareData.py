#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Script to investigate CO2 Concentrations
# June 2020
# author: E.-M. Schoemann

import datetime
#import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import DetrendMonthdate, getMonthMeansTCCON, getMonthMeans, getTenDayMeans,getTenDaySum, createDataFrameNBE, createDataFrameGFED, getMonthSum, addMonthDate, IndividualDetrend, getReferenceDate
import geopandas
import glob
import xarray as xr
import argparse
import os.path
import math
from CreateAndSaveGeoDataframeGFAS import CreateDataFrameGFAS
from CreateAndSaveGeoDataframeFINN import CreateDataFrameFINN
from CreateAndSaveGeoDataframeGOSATold import CreateDataFrameRTold
from CreateAndSaveGeoDataframeOCO2 import CreateDataFrameOCO
from CreateAndSaveGeoDataframeACOS import CreateDataFrameACOS
from CreateAndSaveGeoDataframeNIES import CreateDataFrameNIES
from CreateAndSaveGeoDataframeTM54DVAR import CreateDataFrameTK5
from CreateAndSaveGeoDataframeMOPITT import CreateDataFrameMOPITT
from CreateAndSaveGeoDataframeOCO2SIF import CreateDataFrameOCOSIF
from CreateAndSaveGeoDataframeCTnoncs import CreateDataFrameCTnoncs
from CreateAndSaveGeoDataframeTCCON import CreateDataFrameTCCON
from CreateAndSaveGeoDataframeTM5noncs import CreateDataFrameTM5noncs
from CreateAndSaveGeoDataframeCAMSnoncs import CreateDataFrameCAMSnoncs
import h5py
from RegionParam import getRegion
import matplotlib as mpl

datapath = "."
savepath = "."

font = {'family' : 'Arial'}
mpl.rc('font', **font)

#SETTINGS
WithOCO2 = False # 
# which Data should be treated
CTm = False
RemotecCO2 = False#False
ACOS = False#False
OCO2 = False#False

TM5_Inversion = False
TM5nocs = True
CAMSnocs = False

CO2cs = False#True
CTcs_comp = False#True
CTcsRT238 = False#True
NIES = False
TM5 = False#TTrue
TM5_2 = False#True

GFED_CO2 = False
GFAS = False#False#True
FINN = False#False#True
TCCONnocs = False; TCCONStationIDList = ['db','wg']; TCCONStationNameList = ['Darwin','Wollongong']

TM5list =['RemoTeCISloc3x2_Inversion','IS3x2_Inversion','ACOSIS3x2_Inversion']
dataRTG = ''
 
# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to Plot a TimeSeries of CO2 concentrations ")
    parser.add_argument("Number", type=int, help="Enter value to choose a region: 949, Australia, 950: semi-arid Australia, 951: not semi-arid Australia")
    parser.add_argument("startdate", type=str, help="enter startdate as string in form 'yyyy-mm-dd'")
    parser.add_argument("enddate", type=str, help="enter enddate as string in form 'yyyy-mm-dd'")    
    parser.add_argument("plot_type", type=int, help="type of plot which should be created: 0 omly avail.")
    args = parser.parse_args()
    return args

# SETTINGS:
args = parse_arguments()

year_min = int(args.startdate[0:4])#2009
month_min = int(args.startdate[5:7])#4
day_min = int(args.startdate[8:10])#1
year_max = int(args.enddate[0:4])#2019
month_max = int(args.enddate[5:7])#12
day_max = int(args.enddate[8:10])#31


Numm = args.Number


RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

# -------------------------------------------------------------------------
# prepare reference datelist to plot: Add Nan values for month without data
# --------------------------------------------------------------------------

DateRef = getReferenceDate(year_min,year_max,month_min,month_max)


dataTK = ''
datap = ''
datapA = ''


# -----------------------------------------------------
# get and prepare GOSAT data
# -----------------------------------------------------

# ACOS GOSAT CO2 DATA
# -----------------------------------------------------
if ACOS:
    if not os.path.isfile(datapath + "/ACOS/DataFrames/GDF6_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameACOS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)

    datap = ''
    gdfACOS = pd.read_pickle(datapath + "/ACOS/DataFrames/GDF6_"+RegionName+str(Numm)+".pkl")
    gdfACOS = gdfACOS[(gdfACOS.quality == 0)] #only measurements with good quality
    
    if False:#os.path.isfile(datapath + "/ACOS/DataFrames/monthFrames/MonthMeans_ACOSCO2_"+datap+RegionName+ str(Numm)+".pkl"):
        Month_means_ACOS = pd.read_pickle(datapath + "/ACOS/DataFrames/monthFrames/MonthMeans_ACOSCO2_"+datap+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_ACOS_1 = getMonthMeans(gdfACOS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2','CO2_error','CO2_uncorr')

        Month_means_ACOS = DetrendMonthdate(Month_means_ACOS_1,2009,5,0,0,0,'CO2','ACOS')
        Month_means_ACOS_uncorr = DetrendMonthdate(Month_means_ACOS_1,2009,5,0,0,0,'CO2_uncorr','ACOS')

        Month_means_ACOS.insert(loc=1,column='Detrend_uncorr',value=Month_means_ACOS_uncorr.Detrend)
        # move detrended data on y-axis to zero
        # get mean values
        meanA = Month_means_ACOS.Detrend.mean()
        mean_uncorrA = Month_means_ACOS.Detrend_uncorr.mean()

        Month_means_ACOS.insert(loc=1,column='Detrend0',value=(Month_means_ACOS.Detrend- meanA))
        Month_means_ACOS.insert(loc=1,column='Detrend_uncorr0',value=(Month_means_ACOS.Detrend_uncorr- mean_uncorrA))

        # calculate IAV
        YearaverageACOS = Month_means_ACOS.groupby(['Month'])['Detrend'].mean().reset_index()
        IAVlist = []
        for index, row in Month_means_ACOS.iterrows():
            IAVlist.append(row["Detrend"]-YearaverageACOS[(YearaverageACOS.Month == row["Month"])].Detrend)
        Month_means_ACOS["IAV"] = IAVlist
        Month_means_ACOS.to_pickle(datapath + "/ACOS/DataFrames/monthFrames/MonthMeans_ACOSCO2_"+datap+RegionName+ str(Numm)+".pkl")

    PMonth_means_ACOS = pd.merge(DateRef, Month_means_ACOS, on=['MonthDate'], how = 'left')


#RemoTeC
# ---------------------------------------------------

if RemotecCO2:
    gdf = pd.read_pickle(datapath + "/GOSAT_Markus/dataframes/GDF10_19_8.pkl")
    gdf = gdf[(gdf.meas_geom == '0')]
    gdf = gdf.drop_duplicates(keep = 'first')
    
    if Numm >= 900:
        if True:
            Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
            gdf = gdf[(gdf.Long >= Long_min)&(gdf.Long <= Long_max)&(gdf.Lat >= Lat_min)&(gdf.Lat <= Lat_max)]
            igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
            gdf = gdf.loc[igdf]
            gdf.to_pickle(datapath + "/GOSAT_Markus/dataframes/GDF_"+RegionName+str(Numm)+".pkl")
        
        gdf = pd.read_pickle(datapath + "/GOSAT_Markus/dataframes/GDF_"+RegionName+str(Numm)+".pkl")
        gdf = gdf[(gdf.meas_geom == '0')]
        gdf = gdf.drop_duplicates(keep = 'first')

    Ggain = 'no'#'MM' #'no'

    if 'no' not in Ggain:
        gdf = gdf[(gdf.gain == Ggain)]
    if False:#os.path.isfile(datapath + "/GOSAT_Markus/dataframes/monthFrames/MonthMeans_GOSAT_"+datapA+RegionName+ str(Numm)+".pkl") and 'no' in Ggain:
        Month_means_CO2 = pd.read_pickle(datapath + "/GOSAT_Markus/dataframes/monthFrames/MonthMeans_GOSAT_"+datapA+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_CO2_1 = getMonthMeans(gdf,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2','CO2_error','CO2_uncorr')          
        
        Month_means_CO2 = DetrendMonthdate(Month_means_CO2_1,2009,5,0,0,0,'CO2','RTnew')#'CH4')
        
        Month_means_CO2_uncorr = DetrendMonthdate(Month_means_CO2_1,2009,5,0,0,0,'CO2_uncorr','RTnew')
        Month_means_CO2.insert(loc=1,column='Detrend_uncorr',value=Month_means_CO2_uncorr.Detrend)

        Month_means_CO2_linDetr = DetrendMonthdate(Month_means_CO2_1,2009,5,0,3,0,'CO2','RTnew')#'CH4')
        Month_means_CO2.insert(loc=1,column='Detrend_lin',value=Month_means_CO2_linDetr.Detrend)

        try:
            Month_means_CO2 = Month_means_CO2.drop(index = list(Month_means_CO2.MonthDate).index(datetime.date(2015,1,15)))
        except:
            pass

        # move detrended data on y-axis to zero
        # get mean values
        meanR = Month_means_CO2.Detrend.mean()
        mean_uncorrR = Month_means_CO2.Detrend_uncorr.mean()
    
        Month_means_CO2.insert(loc=1,column='Detrend0',value=(Month_means_CO2.Detrend- meanR))
        Month_means_CO2.insert(loc=1,column='Detrend_uncorr0',value=(Month_means_CO2.Detrend_uncorr- mean_uncorrR))

        # calculate IAV
        YearaverageRT = Month_means_CO2.groupby(['Month'])['Detrend'].mean().reset_index()
    
        IAVlist = []
        for index, row in Month_means_CO2.iterrows():
            IAVlist.append(row["Detrend"]-YearaverageRT[(YearaverageRT.Month == row["Month"])].Detrend)
        Month_means_CO2["IAV"] = IAVlist
        if 'no' in Ggain:
            Month_means_CO2.to_pickle(datapath + "/GOSAT_Markus/dataframes/monthFrames/MonthMeans_GOSAT_"+datapA+RegionName+ str(Numm)+".pkl")
        
    PMonth_means_CO2 = pd.merge(DateRef, Month_means_CO2, on=['MonthDate'], how = 'left')


if TM5nocs:
    PMonth_means_TM5nocs = []
    TMDtypel = ['IS']#,'RemoTeC+IS-land_ocean_bc']
    for TMDtype in TMDtypel:
        print(TMDtype)
        if not os.path.isfile(datapath + "/TK5_4DVAR/dataframes/Unsampled/GDF_MonthlynoncsMean_"+TMDtype+"_"+RegionName+ str(Numm)+"_merged.pkl"):
            CreateDataFrameTM5noncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, year_min,year_max, month_min,month_max, TMDtype)
        dat = pd.read_pickle(datapath + "/TK5_4DVAR/dataframes/Unsampled/GDF_MonthlynoncsMean_"+TMDtype+"_"+RegionName+ str(Numm)+"_merged.pkl")
    
        if False:
            Month_means_TM5nocs = pd.read_pickle(datapath + "/TK5_4DVAR/dataframes/Unsampled/monthFrames/MonthMeans_TM5_"+TMDtype+"_"+RegionName+ str(Numm)+".pkl")
        else:
            Month_means_TM5nocs = DetrendMonthdate(dat,2009,5,0,0,0,'CO2','TM5nocs')
    
            # move detrended data on y-axis to zero
            meanTM5nocs = Month_means_TM5nocs.Detrend.mean()
            Month_means_TM5nocs.insert(loc=1,column='Detrend0',value=(Month_means_TM5nocs.Detrend- meanTM5nocs))
            Month_means_TM5nocs.to_pickle(datapath + "/TK5_4DVAR/dataframes/Unsampled/monthFrames/MonthMeans_TM5_newSuppPaper_"+TMDtype+"_"+RegionName+ str(Numm)+".pkl")
        
        PMonth_means_TM5nocs.append(pd.merge(DateRef, Month_means_TM5nocs, on=['MonthDate'], how = 'left'))
        
#GFAS
if GFAS:
    if not os.path.isfile(datapath + "/GFAS/dataframes/GDF1_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameGFAS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Month')

    gdfGFAS = pd.read_pickle(datapath + "/GFAS/dataframes/GDF1_"+RegionName+str(Numm)+".pkl")

    if os.path.isfile(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFAS_"+RegionName+ str(Numm)+".pkl") and os.path.isfile(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFASten_"+RegionName+ str(Numm)+".pkl"):
        Month_means_GFAS = pd.read_pickle(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFAS_"+RegionName+ str(Numm)+".pkl")
        Month_means_GFASt = pd.read_pickle(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFASten_"+RegionName+ str(Numm)+".pkl")

    else:
        Month_means_GFASCO21 = getMonthSum(gdfGFAS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2fireE')
        
        Month_means_GFASCO21 = addMonthDate(Month_means_GFASCO21)

        try:
            Month_means_GFASCO1 = getMonthSum(gdfGFAS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'COfireE')

            Month_means_GFASC1 = getMonthSum(gdfGFAS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CfireE')

            Month_means_GFASCO1 = addMonthDate(Month_means_GFASCO1)
            Month_means_GFASC1 = addMonthDate(Month_means_GFASC1)


            Month_means_GFAS = pd.merge(Month_means_GFASCO21, Month_means_GFASCO1[['MonthDate','COfireE']], on='MonthDate', how = 'left')
            Month_means_GFAS = pd.merge(Month_means_GFAS, Month_means_GFASC1[['MonthDate','CfireE']], on='MonthDate', how = 'left')

        except:
            Month_means_GFAS = Month_means_GFASCO21.copy()

        Month_means_GFASt = getTenDaySum(gdfGFAS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2fireE')

        Month_means_GFAS.to_pickle(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFAS_"+RegionName+ str(Numm)+".pkl")
        Month_means_GFASt.to_pickle(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFASten_"+RegionName+ str(Numm)+".pkl")


    PMonth_means_GFAS = pd.merge(DateRef, Month_means_GFAS, on=['MonthDate'], how = 'left')


#FINN
if FINN:
    if not os.path.isfile(datapath + "/FINN/dataframes/GDF1_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameFINN(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Month')

    gdfFINN = pd.read_pickle(datapath + "/FINN/dataframes/GDF1_"+RegionName+str(Numm)+".pkl")

    if os.path.isfile(datapath + "/FINN/dataframes/monthFrames/MonthMeans_FINN_"+RegionName+ str(Numm)+".pkl"):
        Month_means_FINN = pd.read_pickle(datapath + "/FINN/dataframes/monthFrames/MonthMeans_FINN_"+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_FINNCO21 = getMonthSum(gdfFINN,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2fireE')

        Month_means_FINN = addMonthDate(Month_means_FINNCO21)
        

        Month_means_FINN.to_pickle(datapath + "/FINN/dataframes/monthFrames/MonthMeans_FINN_"+dataTK+RegionName+ str(Numm)+".pkl")

    PMonth_means_FINN = pd.merge(DateRef, Month_means_FINN, on=['MonthDate'], how = 'left')


#CAMS non sampled
if CAMSnocs:
    print('CAMSnocs')
    #only done for surface CAMS
    PMonth_means_CAMSnocs = []
    if not os.path.isfile(datapath + "/CAMS/dataframes/GDF_Monthlynoncs3_"+RegionName+ str(Numm)+".pkl"):
        CreateDataFrameCAMSnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
    dat = pd.read_pickle(datapath + "/CAMS/dataframes/GDF_Monthlynoncs3_"+RegionName+ str(Numm)+".pkl")
    dat.insert(loc=1,column='XCO2',value=dat.CO2*1000000)
    if False:#os.path.isfile(datapath + "/CAMS/dataframes/monthFrames/MonthMeans_CAMSnoncs_"+RegionName+ str(Numm)+".pkl"): 
        Month_means_CAMSnocs = pd.read_pickle(datapath + "/CAMS/dataframes/monthFrames/MonthMeans_CAMSnoncs_"+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_CAMSnocs1 = getMonthMeans(dat,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'XCO2')
        
        Month_means_CAMSnocs = DetrendMonthdate(Month_means_CAMSnocs1,2009,5,0,0,0,'XCO2','CAMSnocs')

        # move detrended data on y-axis to zero
        meanCAMSnocs = Month_means_CAMSnocs.Detrend.mean()
        Month_means_CAMSnocs.insert(loc=1,column='Detrend0',value=(Month_means_CAMSnocs.Detrend- meanCAMSnocs))
        Month_means_CAMSnocs.to_pickle(datapath + "/CAMS/dataframes/monthFrames/MonthMeans_CAMSnoncs_"+RegionName+ str(Numm)+".pkl")
    PMonth_means_CAMSnocs = pd.merge(DateRef, Month_means_CAMSnocs, on=['MonthDate'], how = 'left')


# CT CO2 
# ----------------------------------------------------
# not cosampled Monthly val
if CTm:
    print('CTm')
    if not os.path.isfile(datapath + "/CT2019/dataframes/GDF_Monthlynoncs2_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameCTnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
    
    gdfCTm = pd.read_pickle(datapath + "/CT2019/dataframes/GDF_Monthlynoncs2_"+RegionName+str(Numm)+".pkl")


    if False:#os.path.isfile(datapath + "/CT2019/dataframes/monthFrames/MonthMeans_CTm_"+RegionName+ str(Numm)+".pkl"):
        Month_means_CTm = pd.read_pickle(datapath + "/CT2019/dataframes/monthFrames/MonthMeans_CTm_"+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_CTm_1 = getMonthMeans(gdfCTm,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2')
        
        Month_means_CTm = DetrendMonthdate(Month_means_CTm_1,2009,5,0,0,0,'CO2','CTm')

        # move detrended data on y-axis to zero
        meanCTm = Month_means_CTm.Detrend.mean()
        Month_means_CTm.insert(loc=1,column='Detrend0',value=(Month_means_CTm.Detrend- meanCTm))
        
        # calculate IAV
        YearaverageCTm = Month_means_CTm.groupby(['Month'])['Detrend'].mean().reset_index()
        IAVlist = []
        for index, row in Month_means_CTm.iterrows():
            IAVlist.append(row["Detrend"]-YearaverageCTm[(YearaverageCTm.Month == row["Month"])].Detrend)
        Month_means_CTm["IAV"] = IAVlist
        Month_means_CTm.to_pickle(datapath + "/CT2019/dataframes/monthFrames/MonthMeans_CTm_"+RegionName+ str(Numm)+".pkl")

    PMonth_means_CTm = pd.merge(DateRef, Month_means_CTm, on=['MonthDate'], how = 'left')

# -----------------------------------------------------
# get and prepare OCO-2 data
# -----------------------------------------------------
if OCO2:
    if not os.path.isfile(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl"):
        CreateDataFrameOCO(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)

    gdfOCO = pd.read_pickle(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl")
    
    if Numm >= 900:
        if not os.path.isfile(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl"):
            Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
            gdfOCO = gdfOCO[(gdfOCO.Long >= Long_min)&(gdfOCO.Long <= Long_max)&(gdfOCO.Lat >= Lat_min)&(gdfOCO.Lat <= Lat_max)]
            igdfOCO = gdfOCO.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
            gdfOCO = gdfTK5.loc[igdfOCO]
            gdfOCO.to_pickle(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl")

        gdfOCO = pd.read_pickle(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Numm)+".pkl")

    #only good flag
    gdfOCO = gdfOCO[gdfOCO.quality == 0]

    if False:#os.path.isfile(datapath + "/OCO-2/DataFrames/monthFrames/MonthMeans_OCO2CO2_"+RegionName+ str(Numm)+".pkl"):
        Month_means_OCO = pd.read_pickle(datapath + "/OCO-2/DataFrames/monthFrames/MonthMeans_OCOCO2_"+RegionName+ str(Numm)+".pkl")
    else:
        Month_means_OCO_1 = getMonthMeans(gdfOCO,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2','CO2_error')

        Month_means_OCO = DetrendMonthdate(Month_means_OCO_1,2009,5,0,0,0,'CO2')#,'OCO')
        Month_means_OCO_C = DetrendMonthdate(Month_means_OCO_1,2009,5,0,0,0,'CO2','OCO')
        Month_means_OCO_C = Month_means_OCO_C.rename(columns={"Detrend": "DetrendCAMS"},errors='raise')
        Month_means_OCO = pd.merge(Month_means_OCO, Month_means_OCO_C[['MonthDate','DetrendCAMS']], on='MonthDate', how = 'left')

        # move detrended data on y-axis to zero
        meanOCO2 = Month_means_OCO.Detrend.mean()
        Month_means_OCO.insert(loc=1,column='Detrend0',value=(Month_means_OCO.Detrend- meanOCO2))
        

        # calculate IAV
        YearaverageOCO = Month_means_OCO.groupby(['Month'])['Detrend'].mean().reset_index()
        IAVlist = []
        for index, row in Month_means_OCO.iterrows():
            IAVlist.append(row["Detrend"]-YearaverageOCO[(YearaverageOCO.Month == row["Month"])].Detrend)
        Month_means_OCO["IAV"] = IAVlist
        
    PMonth_means_OCO = pd.merge(DateRef, Month_means_OCO, on=['MonthDate'], how = 'left')


#get and prepare GFED
# --------------------------------------------------------

if GFED_CO2:
    if os.path.isfile(datapath + "/GFED4/MonthMeansCO2_"+RegionName+str(Numm)+".pkl"):
        GFED_CO2e = pd.read_pickle(datapath + "/GFED4/MonthMeansCO2_"+RegionName+str(Numm)+".pkl")
        #GFED_CO2et = pd.read_pickle(datapath + "/GFED4/MonthMeansCO2ten_"+RegionName+str(Numm)+".pkl")
        print('got GFED df')
    else:
        GFED = pd.read_pickle(datapath + "/GFED4/GFED_CO2.pkl")
        GFED = GFED[(GFED.Long >= Long_min) & (GFED.Long <= Long_max) & (GFED.Lat >= Lat_min) & (GFED.Lat <= Lat_max)]
        print('got GFED')

        if Numm >= 900:
            if not os.path.isfile(datapath + "/GFED4/GDF_CO2_"+RegionName+str(Numm)+".pkl"):
                GFED = GFED[(GFED.Long >= Long_min)&(GFED.Long <= Long_max)&(GFED.Lat >= Lat_min)&(GFED.Lat <= Lat_max)]
                gGFED = geopandas.GeoDataFrame(
                      GFED, geometry=geopandas.points_from_xy(GFED.Long, GFED.Lat))
                gGFED.crs = {'init' :'epsg:4326'}

                Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
                igGFED = gGFED.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
                gGFED = gGFED.loc[igGFED]
                gGFED.to_pickle(datapath + "/GFED4/GDF_CO2_"+RegionName+str(Numm)+".pkl")

            GFED = pd.read_pickle(datapath + "/GFED4/GDF_CO2_"+RegionName+str(Numm)+".pkl")
        if 'Grid_area' not in GFED.keys():
            GFED.insert(loc=1,column='Grid_area',value= (np.cos(GFED.Lat*math.pi/180)*773200902))
        if 'total_emission' not in GFED.keys():
            GFED.insert(loc=1, column='total_emission', value = GFED.Grid_area*GFED.emission/(1000000000000))

        GFED_CO2e = getMonthSum(GFED,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'total_emission')

        md = []
        for num,row2 in GFED_CO2e.iterrows():
            md.append(datetime.date(int(row2.Year),int(row2.Month),15))
        GFED_CO2e.insert(loc=1,column='MonthDate',value= md)
        print('calculated monthly means')
        GFED_CO2e.to_pickle(datapath + "/GFED4/MonthMeansCO2_"+RegionName+str(Numm)+".pkl")


    PGFED_CO2e = pd.merge(DateRef, GFED_CO2e, on=['MonthDate'], how = 'left')


if TCCONnocs:
    
    PMonth_means_TCCONnocs = []
    for TCCONStationID in TCCONStationIDList:
        if not os.path.isfile(datapath + "/TCCON/DataFrames/DF2_"+TCCONStationID+".pkl"):
            print('Create TCCON DF')
            CreateDataFrameTCCON(TCCONStationID)
        dfTCCON = pd.read_pickle(datapath + "/TCCON/DataFrames/DF2_"+TCCONStationID+".pkl")
        
        if False:#os.path.isfile(datapath + "/TCCON/DataFrames/monthFrames/Month_means_noCosampling_"+TCCONStationID+".pkl"): 
            Month_means_TCCONnocs = pd.read_pickle(datapath + "/TCCON/DataFrames/monthFrames/Month_means_noCosampling_"+TCCONStationID+".pkl")
        else:
            print('Get TCCON monthly means')
            #dfTCCON = dfTCCON[(dfTCCON.Hour >= 1)&(dfTCCON.Hour <=4)] #11:30 - 14:30 in Darwin -> ca. time of overflight of GoSAT
            Month_means_TCCONnocs1 = getMonthMeansTCCON(dfTCCON,year_min,month_min,day_min,year_max,month_max,day_max,'CO2',10,'CO2_error')
            
            Month_means_TCCONnocs = DetrendMonthdate(Month_means_TCCONnocs1,2009,5,0,0,0,'CO2','TCCON'+TCCONStationID)
    
            # move detrended data on y-axis to zero
            meanTCCONnocs = Month_means_TCCONnocs.Detrend.mean()
            Month_means_TCCONnocs.insert(loc=1,column='Detrend0',value=(Month_means_TCCONnocs.Detrend- meanTCCONnocs))
            #Month_means_TCCONnocs.to_pickle(datapath + "/TCCON/DataFrames/monthFrames/Month_means_noCosampling_"+TCCONStationID+".pkl")
        PMonth_means_TCCONnocs.append(pd.merge(DateRef, Month_means_TCCONnocs, on=['MonthDate'], how = 'left'))

if CTm and TM5nocs and CAMSnocs:#mean of the nocosampled inverse models
    detrendkind = 'Detrend0'
    allModel = PMonth_means_TM5nocs[0][['MonthDate',detrendkind]].copy()
    MonthlistM = allModel.apply(lambda x: x.MonthDate.month,axis=1)
    allModel.insert(loc = 1, column = 'Month', value = MonthlistM)
    print(allModel.MonthDate.max())
    allModel = allModel.rename(columns = {detrendkind:'TM5IS'})
    allModel = pd.merge(allModel,PMonth_means_CAMSnocs[['MonthDate',detrendkind]])
    allModel = allModel.rename(columns = {detrendkind:'CAMS'})
    allModel = pd.merge(allModel,PMonth_means_CTm[['MonthDate',detrendkind]])
    allModel = allModel.rename(columns = {detrendkind:'CT'})
    allModel.insert(loc = 1, column = 'mean_CO2', value = allModel[['TM5IS','CAMS','CT']].mean(axis = 1))
    allModel.insert(loc = 1, column = 'std_CO2', value = allModel[['TM5IS','CAMS','CT']].std(axis = 1,ddof = 0))
    allModel.insert(loc = 1, column = 'min_CO2', value = allModel[['TM5IS','CAMS','CT']].min(axis = 1))
    allModel.insert(loc = 1, column = 'max_CO2', value = allModel[['TM5IS','CAMS','CT']].max(axis = 1))
    allModelYmean = allModel.groupby('Month')['mean_CO2'].mean().reset_index()
    allModelYstd = allModel.groupby('Month')['mean_CO2'].std(ddof = 0).reset_index()
    print(allModel.MonthDate.max())
    print(allModel[(allModel.MonthDate == datetime.date(2020,6,15))].mean_CO2.max())

if ACOS and RemotecCO2:#mean of the GOSAT data
    detrendkind = 'Detrend0'
    allGOSAT = PMonth_means_ACOS[['MonthDate',detrendkind]].copy()
    MonthlistG = allGOSAT.apply(lambda x: x.MonthDate.month,axis=1)
    allGOSAT.insert(loc = 1, column = 'Month', value = MonthlistG)
    allGOSAT = allGOSAT.rename(columns = {detrendkind:'ACOS'})
    allGOSAT = pd.merge(allGOSAT,PMonth_means_CO2[['MonthDate',detrendkind]])
    allGOSAT = allGOSAT.rename(columns = {detrendkind:'RT'})
    varlist = ['ACOS','RT']
    if OCO2 and WithOCO2:
        allGOSAT = pd.merge(allGOSAT,PMonth_means_OCO[['MonthDate',detrendkind]])
        allGOSAT = allGOSAT.rename(columns = {detrendkind:'OCO2'})
        varlist = ['ACOS','RT','OCO2']
    allGOSAT.insert(loc = 1, column = 'mean_CO2', value = allGOSAT[varlist].mean(axis = 1))
    allGOSAT.insert(loc = 1, column = 'std_CO2', value = allGOSAT[varlist].std(axis = 1,ddof = 0))
    allGOSAT.insert(loc = 1, column = 'min_CO2', value = allGOSAT[varlist].min(axis = 1))
    allGOSAT.insert(loc = 1, column = 'max_CO2', value = allGOSAT[varlist].max(axis = 1))
    allGOSATYmean = allGOSAT.groupby('Month')['mean_CO2'].mean().reset_index()
    allGOSATYstd = allGOSAT.groupby('Month')['mean_CO2'].std(ddof = 0).reset_index()
   
    
# create Plot

if False: # Plot 1 Paper
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(1, 2, 
                            figsize = (18.0*cm,6.98*cm),
                            gridspec_kw={'width_ratios': [20, 4],'height_ratios':[1]})
    #Panel a)
    #minmax
    ax2[0].fill_between(allGOSAT.MonthDate, allGOSAT.min_CO2, allGOSAT.max_CO2,color = 'red',zorder = 1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[0].fill_between(allModel.MonthDate, allModel.min_CO2, allModel.max_CO2,color = 'blue',zorder = 1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    #mean
    ax2[0].plot(allGOSAT.MonthDate, allGOSAT.mean_CO2,color = 'darkred',ls = '-',linewidth = lw,marker = '.',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    ax2[0].plot(allModel.MonthDate, allModel.mean_CO2,color = 'blue',ls = '-',linewidth = lw,marker = '.',markersize = 2,zorder = 2, label = r'$\rm Inverse~model_{in-situ}$')# $\rm CO_2$ ')
       
    # Panel b)
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]
    ax2[1].fill_between(range(1,13), allGOSATYmean.mean_CO2.iloc[indexrange]-allGOSATYstd.mean_CO2.iloc[indexrange], allGOSATYmean.mean_CO2.iloc[indexrange]+allGOSATYstd.mean_CO2.iloc[indexrange],color = 'red',zorder = 1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[1].fill_between(range(1,13), allModelYmean.mean_CO2.iloc[indexrange]-allModelYstd.mean_CO2.iloc[indexrange], allModelYmean.mean_CO2.iloc[indexrange]+allModelYstd.mean_CO2.iloc[indexrange],color = 'blue',zorder = 1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2[1].plot(range(1,13), allGOSATYmean.mean_CO2.iloc[indexrange],color = 'darkred',ls = '-',linewidth = lw,marker = '.',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    ax2[1].plot(range(1,13), allModelYmean.mean_CO2.iloc[indexrange],color = 'blue',ls = '-',linewidth = lw,marker = '.',markersize = 2,zorder = 2, label = r'Inverse insitu Models')# $\rm CO_2$ ')
    
    #SETTINGS
    # Panel a)
    ax2[0].text(0.03, 0.95, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10)#, weight='bold')
    ax2[0].legend(fontsize=6,loc = 1)
    ##ax2[0].set_ylim(-1.65,1.8)
    ax2[0].set_ylabel(r'$\rm Detrended~CO_2~(ppm)$',fontsize=7)
    ax2[0].tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2[0].tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2[0].set_xlabel(r'Date',fontsize=7)
    ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    #ax2[0].set_xlim(datetime.date(2008,12,15),datetime.date(2022,12,31))
    ax2[0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0].set_axisbelow(True)

    # Panel b) 
    ax2[1].text(0.15, 0.95, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10)#, weight='bold')      
    secaxD = ax2[1].twinx()   
    secaxD.set_ylim(ax2[0].get_ylim())
    ##secaxD.set_ylim(-1.65,1.8)
    secaxD.set_yticklabels([])
    secaxD.tick_params(axis = 'y',direction = 'in')
    ax2[1].set_xticks(ticks = [3,6,9,12])
    ax2[1].set_xticklabels(['9','12','3','6'])
    ax2[1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1,1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1].set_xlabel(r'Month',fontsize=7)
    ax2[1].set_ylim(ax2[0].get_ylim())
    ##ax2[1].set_ylim(-1.65,1.8)
    ax2[1].set_yticklabels([])
    ax2[1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1].tick_params(axis = 'x',labelsize = 6, zorder = 0)
    ax2[1].set_xlim(0.5,12.5)
    ax2[1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig1complete2_3x2res_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')#,format = 'pdf'
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig1complete2_3x2res_"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')#,format = 'pdf'


#Paper figure Supp with OCO-2 S1   
if False:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6.98*cm))
    ax2.fill_between(allGOSAT.MonthDate, allGOSAT.min_CO2, allGOSAT.max_CO2,color = 'red',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2.fill_between(allModel.MonthDate, allModel.min_CO2, allModel.max_CO2,color = 'blue',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    ax2.plot(allGOSAT.MonthDate, allGOSAT.mean_CO2,color = 'darkred',ls = '-',linewidth= lw,marker = '.',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    ax2.plot(allModel.MonthDate, allModel.mean_CO2,color = 'blue',ls = '-',linewidth= lw,marker = '.',markersize = 2,zorder = 2, label = r'$\rm Inverse~model_{in-situ}$')# $\rm CO_2$ ')
    ax2.plot(PMonth_means_OCO.MonthDate,PMonth_means_OCO.Detrend0,color = 'black',linewidth= lw,marker = '.',markersize = 2,zorder = 2,ls = '-',label = 'OCO-2')# individual detrend')#CO2 NOAA + CAMS detrended')
    savestr = 'Supp2_OCO2_4_2_3x2'
    
    ax2.set_ylabel(r'$\rm Detrended~CO_2~(ppm)$',fontsize=7)
    ax2.tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2.tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.legend(fontsize=6,loc = 1)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_axisbelow(True)

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".png", dpi=400, bbox_inches = "tight")
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".pdf", dpi=400, bbox_inches = "tight",format = 'pdf')

#Paper figure Supp with TCCON, S2    
if False:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6.98*cm))
    
    ax2.fill_between(allGOSAT.MonthDate, allGOSAT.min_CO2, allGOSAT.max_CO2,color = 'red',zorder =1,alpha = 0.3)#, label = r'GOSAT mean')# $\rm CO_2$ ')
    
    ax2.plot(allGOSAT.MonthDate, allGOSAT.mean_CO2,color = 'darkred',ls = '-',linewidth= lw,marker = '.',markersize = 2,zorder = 2, label = r'GOSAT')# $\rm CO_2$ ')
    
    ax2.plot(PMonth_means_TCCONnocs[0].MonthDate,PMonth_means_TCCONnocs[0].Detrend0,color = 'grey',linewidth= lw,marker = '.',markersize = 2,zorder = 2,ls = '--',label = 'TCCON '+TCCONStationNameList[0])# individual detrend')#CO2 NOAA + CAMS detrended')
    ax2.plot(PMonth_means_TCCONnocs[1].MonthDate,PMonth_means_TCCONnocs[1].Detrend0,color = 'black',linewidth= lw,marker = '.',markersize = 2,zorder = 2,ls = ':',label = 'TCCON '+TCCONStationNameList[1])# individual detrend')#CO2 NOAA + CAMS detrended')
    
    savestr = 'Supp2_TCCON_4_3x2'
    
    ax2.set_ylabel(r'$\rm Detrended~CO_2~(ppm)$',fontsize=7)
    ax2.tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2.tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.legend(fontsize=6,loc = 1)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_axisbelow(True)

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".png", dpi=400, bbox_inches = "tight")
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".pdf", dpi=400, bbox_inches = "tight",format = 'pdf')

#Fire plot S8
if True:
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6.98*cm))
    
    #amplified!!!!!
    print('plotting')
    ax2.plot(PMonth_means_FINN.MonthDate,PMonth_means_FINN.CO2fireE*10**(-12)*10*12/44,color = 'purple',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '--',label = 'FINN * 10')
    ax2.plot(PMonth_means_GFAS.MonthDate,PMonth_means_GFAS.CO2fireE/1000000000000*12/44,color = 'red',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'GFAS')
    ax2.plot(PGFED_CO2e.MonthDate,PGFED_CO2e.total_emission*12/44,color = 'orange',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'GFED')#r'$\rm GFED~fire~CO_2$ emissions')
    ax2.plot(PMonth_means_FINN.MonthDate,PMonth_means_FINN.CO2fireE*10**(-12)*12/44,color = 'purple',linewidth= lw,marker = '.',markersize = 2,zorder = 2, ls = '-',label = 'FINN')
    savestr = 'Supp_GFASFINNGFED_ampliTgC_2'
    ax2.set_ylabel(r'$\rm Fire~CO_2~emission~(TgC/month)$',fontsize=7)
    ax2.tick_params(axis = 'x',length = 0.01, zorder = 0,labelsize = 6)
    ax2.tick_params(axis = 'y',labelsize = 6,zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.legend(fontsize=6,loc = 1)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_axisbelow(True)

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".png", dpi=400, bbox_inches = "tight")
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/"+savestr+"_"+RegionName+str(Numm)+str(year_min)+"_"+str(year_max)+".pdf", dpi=400, bbox_inches = "tight",format = 'pdf')
