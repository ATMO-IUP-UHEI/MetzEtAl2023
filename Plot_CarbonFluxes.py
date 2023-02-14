#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:38:45 2021

@author: eschoema
Script to investigate Carbon Fluxes 

"""
from pathlib import Path
import datetime
from tarfile import REGULAR_TYPES
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import DetrendMonthdate, getMeanAmplitude, getMonthMeans, getTenDayMeans,getTenDaySum, createDataFrameNBE, getMonthSum, addMonthDate, IndividualDetrend, getReferenceDate,  getReferencesDateDay,getMonthSumFCday,getNumDayOfMonth
import geopandas
import glob
import xarray as xr
import argparse
import os.path
import math
from RegionParam import getRegion
from CreateAndSaveGeoDataframeCTFlux import CreateDataFrameCTFlux, CombineGDF,CreateDataFrameCTFlux2
#from CreateAndSaveGeoDataframeFLUXCOM import CreateDataFrameSIFFC
from CreateAndSaveGeoDataframeCAMSFlux1x1 import CreateDataFrameCAMSflux, CombineGDFCAMS
from CreateAndSaveGeoDataframeCABLE import CreateDataFrameCABLE
from CreateAndSaveGeoDataframeMIP import CreateDataFrameMIPflux
from sklearn.linear_model import LinearRegression
import matplotlib.font_manager as fm

datapath = "."
savepath = "."

font = {'family' : 'Arial'}
mpl.rc('font', **font)

WithOCO2 = False

#choose which datasets are getting loaded
CAMS_b = False###
CT_flux = False###
CT_type = 'MonthlyRegion'#''
SIF_FC = False#
TM5Flux = False###
TM5Flux1x1 = False###
TRENDY = True###False#True; 
ERA5precip = False###
Fires = False###
MIPdata = False###

#further specifications, which TRENDY models are taken into the analyses
TrendyModels = ['DLEM', 'IBIS', 'ISAM', 'JSBACH', 
                'LPX-Bern', 'OCN','ORCHIDEE', 'ORCHIDEE-CNP', 'ORCHIDEEv3', 
                'CLM5.0','ISBA-CTRIP','LPJ','SDGVM',
                'VISIT','YIBs','CABLE-POP', 'CLASSIC','JULES-ES-1p0'] #'JULES-ES-1p0'

# define the subselections of TRENDY models
def GoodModelsNames(Datatype):
    return ['JSBACH'+Datatype,
             'LPJ'+Datatype,
             'CLASSIC'+Datatype,
             'YIBs'+Datatype,
             'OCN'+Datatype]            
def EarlyGPPModelNames(Datatype):
    return ['ORCHIDEE-CNP'+Datatype,
             'CABLE-POP'+Datatype,
             'ORCHIDEE'+Datatype,
             'CLM5.0'+Datatype,#
             'VISIT'+Datatype,#
             'DLEM'+Datatype,#
             'ISAM'+Datatype,
             'IBIS'+Datatype]
def LateRespNames(Datatype): 
    return ['LPX-Bern'+Datatype,
             'ISBA-CTRIP'+Datatype,
             #'VISIT'+Datatype,
             #'DLEM'+Datatype,
             'SDGVM'+Datatype,
             'ORCHIDEEv3'+Datatype,
             #'CLM5.0'+Datatype,
             'JULES-ES-1p0'+Datatype]
def BadModelNames(Datatype):
    return EarlyGPPModelNames(Datatype) + LateRespNames(Datatype)

#select the Inversion TM5-4DVar datasets
# for 1x1 resolution
TM51x1dataset = ["flux_1x1_RemoTeC_2.3.8+IS","flux_1x1_RemoTeC_2.4.0+IS","flux_1x1_IS","flux_1x1_ACOS+IS","flux_1x1_LNLGIS","flux_1x1_LNLG"]#,
#for regional fluxes
TM5list = ['RemoTeCISlocglb3x2','ISglb3x2','ACOSISglb3x2']

#input (could be moved to arg parser)
# which region and which period
Numm = 949
RegionArea = 7.965*10**12 ##DryERA5 (950):5422341009727 #949: 7.965*10**12 #wetERA5 (951): 2525018121358 #dryLC(952): 4902147095848 #wetLC (953): 3045212035237
startdate = '2009-01-01'#'2009-01-01'#'2015-01-01'#'2009-04-15'#'2014-08-15'#
enddate = '2018-12-31'
minyearmean = 2010 #needs to be changed for annual flux calculation later starting than 2009
                   # exception: the june-july calculations are starting mid of 2009
    
#main settings
year_min = int(startdate[0:4])#2009
month_min = int(startdate[5:7])#4
day_min = int(startdate[8:10])#1
year_max = int(enddate[0:4])#2019
month_max = int(enddate[5:7])#12
day_max = int(enddate[8:10])#31

RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

DateRef = getReferenceDate(year_min,year_max,month_min,month_max)

#IAV and mean flux values dataframe
dfIAV = pd.DataFrame(data = {'Dataset':[],'AnnualMean':[],'StDev':[],'AnnualMin':[],'AnnualMax':[]})


# get the data sets
# --------------------------------------------------------
if Fires:
    #read GFED
    GFEDregion = pd.read_pickle(datapath + "/GFED4/MonthMeansCO2_"+RegionName+str(Numm)+".pkl")
    PGFEDregion = pd.merge(DateRef, GFEDregion, on=['MonthDate'], how = 'left')
    PGFEDregion.rename(columns={'Month_x':'Month','Year_x':'Year'},inplace=True)
    # get annual mean and Stdev
    GFEDMonthlyAnnualMean = PGFEDregion.groupby('Month')['total_emission'].mean().reset_index()
    GFEDMonthlyAnnualStd = PGFEDregion.groupby('Month')['total_emission'].std(ddof=0).reset_index()
    #get mean Amplitude
    GFEDmeanAmplitude, GFEDmeanAmplitudeStDev = getMeanAmplitude(PGFEDregion,'total_emission','Year',2009,2018)[0:2]
    GFEDmeanAmplitude_no09, GFEDmeanAmplitudeStDev_no09 = getMeanAmplitude(PGFEDregion,'total_emission','Year',2010,2018)[0:2]
    GFEDmeanAmplitudejj, GFEDmeanAmplitudeStDevjj,GFEDmeanAmplitudeMin,GFEDmeanAmplitudeMax = getMeanAmplitude(PGFEDregion,'total_emission','JJYear',2009,2017)
    
    #read GFAS
    FireRegion = pd.read_pickle(datapath + "/GFAS/dataframes/monthFrames/MonthMeans_GFAS_AU902.pkl")
    FireRegion.insert(loc=1, value= FireRegion.CO2fireE/10**12,column = 'CO2fireET')
    PFireRegion = pd.merge(DateRef, FireRegion, on=['MonthDate'], how = 'left')
    MonthlistG = PFireRegion.apply(lambda x: x.MonthDate.month,axis=1)
    PFireRegion.insert(loc = 1, column = 'Month', value = MonthlistG)
    AMeanGFAS = PFireRegion.groupby(['Month'])['CO2fireET'].mean().reset_index()

if MIPdata:
    if False:#os.path.isfile(datapath + "/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl"):
        MIPmean = pd.read_pickle(datapath + "/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl")
    else:
        if os.path.isfile(datapath + "/OCO-2_v10_MIP/dataframes/GDFMIP_"+RegionName+str(Numm)+".pkl"):
            MIP = pd.read_pickle(datapath + "/OCO-2_v10_MIP/dataframes/GDFMIP_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,version = 'regular')
            MIP = pd.read_pickle(datapath + "/OCO-2_v10_MIP/dataframes/GDFMIP_"+RegionName+str(Numm)+".pkl")
    
        MIPmean = getMonthSum(MIP,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot')
        #gaussian errorpropagation for Std: reg Std = sqrt(sum(grid std^2))
        MIPStd = getMonthMeans(MIP,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Landtot','LandStdtot')
        # * number of grid cell as getMonthMeans divides Std by that
        MIPStd.insert(loc = 1, column = 'LandStdtot2',value= MIPStd.LandStdtot*MIPStd.number)
        MIPmonth = pd.merge(MIPmean,MIPStd[['LandStdtot2','Month','Year']], on=['Month','Year'], how = 'left')
        # insert MontDate in datframe
        md = MIPmonth.apply(lambda x: datetime.date(int(x.Year), int(x.Month), 15),axis=1)
        MIPmonth.insert(loc=1,column='MonthDate',value=md)
        MIPmonth.to_pickle(datapath + "/OCO-2_v10_MIP/dataframes/monthFrames/MonthMeansMIP_"+RegionName+str(Numm)+".pkl")
        del MIPStd,MIPmean,MIP,md

    PMIP = pd.merge(DateRef, MIPmonth, on=['MonthDate'], how = 'left')
    del MIPmonth

if CAMS_b:       
    if False:#os.path.isfile(datapath + "/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl"):
        CAMS_region = pd.read_pickle(datapath + "/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl")
    else:
        CAMSGDFversion = '1x1'#'1x1'#'2'
        
        if os.path.isfile(datapath + "/CAMS/dataframes/GDF"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl"):
            CAMS = pd.read_pickle(datapath + "/CAMS/dataframes/GDF"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl")
        else:
            CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Datatype='Sur')
            CombineGDFCAMS(RegionName, Numm,'Sur')
            CAMS = pd.read_pickle(datapath + "/CAMS/dataframes/GDF"+CAMSGDFversion+"FLUXSur_"+RegionName+str(Numm)+".pkl")
    
        #CAMS_region = CAMS[(CAMS.Long <= Long_max) & (CAMS.Long > Long_min) & (CAMS.Lat <= Lat_max) & (CAMS.Lat > Lat_min)]
        CAMS_region_b = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_bio_tot')
        CAMS_region_o = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apos_ocean_tot')
        CAMS_region_ba = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_bio_tot')
        CAMS_region_oa = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_apri_ocean_tot')
        CAMS_region_ff = getMonthSum(CAMS,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'flux_foss_tot')

        CAMS_region = pd.merge(CAMS_region_b,CAMS_region_o, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ba, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_oa, on=['Month','Year'], how = 'left')
        CAMS_region = pd.merge(CAMS_region,CAMS_region_ff, on=['Month','Year'], how = 'left')
        
        md = []
        for num,row2 in CAMS_region.iterrows():
            md.append(datetime.date(int(row2.Year),int(row2.Month),15))
        CAMS_region.insert(loc=1,column='MonthDate',value= md)
        CAMS_region.to_pickle(datapath + "/CAMS/dataframes/monthFrames/MonthMeansFluxSur_"+RegionName+str(Numm)+".pkl")

    PCAMS_region_Sur = pd.merge(DateRef, CAMS_region, on=['MonthDate'], how = 'left')

if TM5Flux1x1:
    TM51x1_regionl = []
    for TM51x1name in TM51x1dataset:
        try:
            TM51x1 = pd.read_pickle(datapath + "/TM5Inversion/glb3x2_20220413/new_gridded_flux/dataframes/GDF2FLUXmonthly1x1_"+TM51x1name+RegionName+str(Numm)+"V2.pkl")
        except:
            TM51x1 = pd.read_pickle(datapath + "/TM5Inversion/glb3x2_20220413/gridded_flux/dataframes/GDF2FLUXmonthly1x1_"+TM51x1name+RegionName+str(Numm)+".pkl")
        
        try:
            TM51x1_region_nee = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_nee_total')
            TM51x1_region_fire = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_fire_total')
            TM51x1_region = pd.merge(TM51x1_region_nee,TM51x1_region_fire, on =['Year','Month'])      
        except:
            TM51x1_region = getMonthSum(TM51x1,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'CO2_flux_nbe_total')
            
        date = TM51x1_region.apply(lambda x: datetime.date(int(x.Year),
                                                        int(x.Month),
                                                        15),axis=1)
        TM51x1_region.insert(loc=1,column='MonthDate',value=date)
        
        PTM51x1_region = pd.merge(DateRef, TM51x1_region, on=['MonthDate'], how = 'left') 
        TM51x1_regionl.append(PTM51x1_region)


if CT_flux:   
    CTF = pd.read_pickle(datapath + "/CT2019/dataframes/GDF1FLUXmonthly1x1_"+RegionName+str(Numm)+".pkl")
    CTF_region_b = getMonthSum(CTF,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Biotot')
    CTF_region_o = getMonthSum(CTF,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'Ocntot')
    CTF_region_ff = getMonthSum(CTF,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'fftot')
    CTF_region_f = getMonthSum(CTF,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'firetot')

    CTF_region = pd.merge(CTF_region_b,CTF_region_o, on=['Month','Year'], how = 'left')
    CTF_region = pd.merge(CTF_region,CTF_region_ff, on=['Month','Year'], how = 'left')
    CTF_region = pd.merge(CTF_region,CTF_region_f, on=['Month','Year'], how = 'left')
    date = CTF_region.apply(lambda x: datetime.date(int(x.Year),
                                                 int(x.Month),
                                                 15),axis=1)
    CTF_region.insert(loc=1,column='MonthDate',value=date)

    PCTF_region = pd.merge(DateRef, CTF_region, on=['MonthDate'], how = 'left')   

if SIF_FC:
    PSIFFC_region = []
    PSIFFC_region2 = []
    for Datatype in ['NEE']:#,'GPP']:#['GPP','NEE']:#,'TER']:
        print(datapath + "/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
        if os.path.isfile(datapath + "/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl"):
            SIFFC_region = pd.read_pickle(datapath + "/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
        else:
            print('get dataframe')
            if os.path.isfile(datapath + "/SIF/FLUXCOM/GDF2_"+Datatype+RegionName+str(Numm)+".pkl"):
                print('yes')
                SIFFC = pd.read_pickle(datapath + "/SIF/FLUXCOM/GDF2_"+Datatype+RegionName+str(Numm)+".pkl")
            else:
                CreateDataFrameSIFFC(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,Datatype)
                SIFFC = pd.read_pickle(datapath + "/SIF/FLUXCOM/GDF2_"+Datatype+RegionName+str(Numm)+".pkl")
        
            SIFFC_region = getMonthSum(SIFFC,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,Datatype+'tot',Datatype+'_madtot')
            SIFFC_region2 = getMonthSumFCday(SIFFC,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,Datatype)
            
            md = []
            for num,row2 in SIFFC_region.iterrows():
                md.append(datetime.date(int(row2.Year),int(row2.Month),15))
            SIFFC_region.insert(loc=1,column='MonthDate',value= md)
            SIFFC_region2.insert(loc=1,column='MonthDate',value= md)
            SIFFC_region.to_pickle(datapath + "/SIF/FLUXCOM/MonthMeans_"+Datatype+RegionName+str(Numm)+".pkl")
            SIFFC_region2.to_pickle(datapath + "/SIF/FLUXCOM/MonthMeans2_"+Datatype+RegionName+str(Numm)+".pkl")        
        print('got FC' + Datatype)
        PSIFFC_region.append(pd.merge(DateRef, SIFFC_region, on=['MonthDate'], how = 'left'))   
          
    FCGFED = pd.merge(PSIFFC_region[0],PGFEDregion,on='MonthDate')
    FCGFED.insert(loc=1, column = 'nbp',value=FCGFED.NEEtot/10**12*(44/12)+FCGFED.total_emission)
    FCGFED = pd.merge(DateRef, FCGFED, on=['MonthDate','Year','Month'], how = 'left')
    FCGFED.to_pickle(datapath + "/SIF/FLUXCOM/MonthMeans_FCGFED_"+RegionName+str(Numm)+"_newForSuppPaper.pkl")
    #get mean Amplitude
    FCGFEDmeanAmplitude, FCGFEDmeanAmplitudeStDev = getMeanAmplitude(FCGFED,'nbp','Year',2009,2018)[0:2]
    FCGFEDmeanAmplitude_no09, FCGFEDmeanAmplitudeStDev_no09 = getMeanAmplitude(FCGFED,'nbp','Year',2010,2018)[0:2]
    FCGFEDmeanAmplitudejj, FCGFEDmeanAmplitudeStDevjj,FCGFEDmeanAmplitudeMin,FCGFEDmeanAmplitudeMax = getMeanAmplitude(FCGFED,'nbp','JJYear',2009,2017)
        
    FCannualMean = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Year'])['nbp'].sum().mean()
    FCannualStDev = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Year'])['nbp'].sum().std(ddof = 0)
    FCMonthlyAnnualMean = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Month'])['nbp'].mean().reset_index()
    FCMonthlyAnnualStDev = FCGFED[(FCGFED.Year >= minyearmean)].groupby(['Month'])['nbp'].std(ddof = 0).reset_index()
    dfIAVFC = pd.DataFrame(data = {'Dataset':['FC+GFED'],
                            'AnnualMean':[FCannualMean],
                            'StDev':[FCannualStDev],
                            'AnnualMin':[np.nan],
                            'AnnualMax':[np.nan]})
    dfIAV = dfIAV.append(dfIAVFC,ignore_index=True)
    
if ERA5precip:
    ERA5 = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/GDF1_AU949.pkl")
    MonthSumERA5 = getMonthMeans(ERA5,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'monthlytp')
    date = MonthSumERA5.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthSumERA5.insert(loc=1,column='MonthDate',value=date)
    MeanMonthERA5 = MonthSumERA5.groupby(['Month'])['monthlytp'].mean().reset_index()

    #Mask region with ERA5 Mask
    masklimitmm = '0.02'
    masklimitmonth = '4'
    MaskDry = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskDry"+masklimitmm+"_"+masklimitmonth+".pkl")
    ERA5dry = geopandas.overlay(ERA5,MaskDry,how='intersection')
    ERA5dry.rename(columns= {'Lat_1':'Lat','Long_1':'Long','monthlytp_1':'monthlytp','Month_1':'Month'},inplace = True)
    ERA5dry.drop(columns=['Lat_2','Long_2'] )
    MonthSumERA5dry = getMonthMeans(ERA5dry,year_min,month_min,day_min,year_max,month_max,day_max,Long_min,Long_max,Lat_min,Lat_max,'monthlytp')
    date = MonthSumERA5dry.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
    MonthSumERA5dry.insert(loc=1,column='MonthDate',value=date)
    
    
if TRENDY:
    PMonthSumTrendy = []
    PMonthMeansTrendy = []
    AllMMTrendyNBP = DateRef.copy()
    MeanCycle = pd.DataFrame(data = {'Month':np.array(range(1,13))})
    firstmodel = True
    for Model in TrendyModels:
        DataList = ['nbp','gpp','ra','rh','fFire']# ['nbp','gpp','ra','rh']#['gpp','nbp','npp','ra','rh','fFire','fLuc']
        for Datatype in DataList:
            if os.path.isfile(datapath + "/TRENDY/dataframes/MonthFrames/MonthMeans_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl"):
                MonthMeansTrendy = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/MonthMeans_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")
                MonthSumTrendy = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/MonthSum_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")
            else:
                print('get dataframe')
                if Model in ['CABLE-POP','DLEM'] and Datatype == 'nbp':
                    print(Model +' does not provide nbp, nbp will be calculated as npp - rh')
                    gdfTrendyNPP = pd.read_pickle(datapath + "/TRENDY/dataframes/GDF1"+Model+"npp_"+RegionName+str(Numm)+".pkl")
                    gdfTrendyRH = pd.read_pickle(datapath + "/TRENDY/dataframes/GDF1"+Model+"rh_"+RegionName+str(Numm)+".pkl")
                  
                    MonthSumTrendyNPP = getMonthSum(gdfTrendyNPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'npptot')
                    MonthMeansTrendyNPP = getMonthMeans(gdfTrendyNPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'npp')
                    MonthSumTrendyRH = getMonthSum(gdfTrendyRH,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'rhtot')
                    MonthMeansTrendyRH = getMonthMeans(gdfTrendyRH,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'rh')
                    
                    MonthMeansTrendy = pd.merge(MonthMeansTrendyNPP,MonthMeansTrendyRH,on=['Month','Year'])
                    #convert units from mugCo2/m2s to TgCO2 in Australia
                    SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                    MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                    # value * SecPerMonth * g/kg * Tg/g * A_Australia * gCO2/gC
                    rhTg = MonthMeansTrendy.rh*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * RegionArea *44/12
                    nppTg = MonthMeansTrendy.npp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * RegionArea * 44/12
                    MonthMeansTrendy.insert(loc = 1, column = 'rhTg', value=rhTg)
                    MonthMeansTrendy.insert(loc = 1, column = 'nppTg', value=nppTg)
                    MonthMeansTrendy.drop(['npp','rh'],axis = 1, inplace = True)
                    MonthMeansTrendy.rename(columns={'nppTg':'npp','rhTg':'rh'}, inplace = True)
                    MonthSumTrendy = pd.merge(MonthSumTrendyNPP,MonthSumTrendyRH,on=['Month','Year'])
                    MonthMeansTrendy.insert(loc = 1, value = MonthMeansTrendy.npp-MonthMeansTrendy.rh,column = 'nbp')
                    MonthSumTrendy.insert(loc = 1, value = MonthSumTrendy.npptot-MonthSumTrendy.rhtot,column = 'nbptot')
                elif Model in ['JULES-ES-1p0'] and Datatype == 'ra':
                    print(Model +' does not provide ra, calculated as gpp - npp as npp=gpp-ra')
                    gdfTrendyNPP = pd.read_pickle(datapath + "/TRENDY/dataframes/GDF1"+Model+"npp_"+RegionName+str(Numm)+".pkl")
                    gdfTrendyGPP = pd.read_pickle(datapath + "/TRENDY/dataframes/GDF1"+Model+"gpp_"+RegionName+str(Numm)+".pkl")
                  
                    MonthSumTrendyNPP = getMonthSum(gdfTrendyNPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'npptot')
                    MonthMeansTrendyNPP = getMonthMeans(gdfTrendyNPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'npp')
                    MonthSumTrendyGPP = getMonthSum(gdfTrendyGPP,year_min,month_min,day_min,
                                                   year_max,month_max,day_max,
                                                   Long_min,Long_max,Lat_min,Lat_max,
                                                   'gpptot')
                    MonthMeansTrendyGPP = getMonthMeans(gdfTrendyGPP,year_min,month_min,day_min,
                                                      year_max,month_max,day_max,
                                                      Long_min,Long_max,Lat_min,Lat_max,
                                                      'gpp')
                    MonthMeansTrendy = pd.merge(MonthMeansTrendyNPP,MonthMeansTrendyGPP,on=['Month','Year'])
                    #convert units from kgCo2/m2s to TgCO2 in Australia
                    SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                    MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                    # value * SecPerMonth * g/kg * Tg/g * A_Australia * gCo2/gC 
                    gppTg = MonthMeansTrendy.gpp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * RegionArea * 44/12
                    nppTg = MonthMeansTrendy.npp*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * RegionArea * 44/12
                    MonthMeansTrendy.insert(loc = 1, column = 'gppTg', value=gppTg)
                    MonthMeansTrendy.insert(loc = 1, column = 'nppTg', value=nppTg)
                    MonthMeansTrendy.drop(['npp','gpp'],axis = 1, inplace = True)
                    MonthMeansTrendy.rename(columns={'nppTg':'npp','gppTg':'gpp'}, inplace = True)
                    MonthSumTrendy = pd.merge(MonthSumTrendyNPP,MonthSumTrendyGPP,on=['Month','Year'])
                    MonthMeansTrendy.insert(loc = 1, value = MonthMeansTrendy.gpp-MonthMeansTrendy.npp,column = 'ra')
                    MonthSumTrendy.insert(loc = 1, value = MonthSumTrendy.gpptot-MonthSumTrendy.npptot,column = 'ratot')
                else:
                    if Model in ['CABLE-POP','DLEM','IBIS','ISAM','JULES-ES-1p0','OCN','ORCHIDEEv3','YIBs'] and Datatype == 'fFire':#models without fire emissions
                        nanarray = np.empty(len(MonthSumTrendy.MonthDate))
                        nanarray = np.nan
                        MonthSumTrendy.insert(loc = 1, column = 'fFire', value=nanarray)
                        MonthMeansTrendy.insert(loc = 1, column = 'fFire', value=nanarray)
                        MonthSumTrendy = MonthSumTrendy[['fFire','Year','Month']]
                        MonthMeansTrendy = MonthMeansTrendy[['fFire','Year','Month']]
                    else:
                        gdfTrendy = pd.read_pickle(datapath + "/TRENDY/dataframes/GDF1"+Model+Datatype+"_"+RegionName+str(Numm)+".pkl")
                    
                        MonthSumTrendy = getMonthSum(gdfTrendy,year_min,month_min,day_min,
                                                    year_max,month_max,day_max,
                                                    Long_min,Long_max,Lat_min,Lat_max,
                                                    Datatype+'tot')
                        MonthMeansTrendy = getMonthMeans(gdfTrendy,year_min,month_min,day_min,
                                                        year_max,month_max,day_max,
                                                        Long_min,Long_max,Lat_min,Lat_max,
                                                        Datatype)
                    if Datatype is not 'mrso':
                        print('convert units from mugCo2/m2s to TgCO2 in Australia')
                        SpM = MonthMeansTrendy.apply(lambda x: 60*60*24*getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1)
                        MonthMeansTrendy.insert(loc = 1, column = 'SecInMonth', value=SpM)
                        # value * SecPerMonth * g/kg * Tg/g * A_Australia *gCO2/gC
                        varTg = MonthMeansTrendy[Datatype]*MonthMeansTrendy.SecInMonth * 10**(3) * 10**(-12) * RegionArea * 44/12
                        MonthMeansTrendy.insert(loc = 1, column = 'varTg', value=varTg)
                        MonthMeansTrendy.drop([Datatype],axis = 1, inplace = True)
                        MonthMeansTrendy.rename(columns={'varTg':Datatype}, inplace = True)
                Monthdate = MonthSumTrendy.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
                MonthSumTrendy.insert(loc=1,column='MonthDate',value=Monthdate)
                Monthdate = MonthMeansTrendy.apply(lambda x: datetime.date(int(x.Year),int(x.Month),15),axis=1)
                MonthMeansTrendy.insert(loc=1,column='MonthDate',value=Monthdate)
                MonthMeansTrendy.to_pickle(datapath + "/TRENDY/dataframes/MonthFrames/MonthMeans_"+Model+"_"+Datatype+RegionName+str(Numm)+"_newForSuppPaper.pkl")
                MonthSumTrendy.to_pickle(datapath + "/TRENDY/dataframes/MonthFrames/MonthSum_"+Model+"_"+Datatype+RegionName+str(Numm)+".pkl")                
                print('got Trendy ' + Datatype)
            MeanCycle = MeanCycle.merge(MonthMeansTrendy.groupby(['Month'])[Datatype].mean(), on = 'Month')
            MeanCycle = MeanCycle.rename(columns = {Datatype: Model+Datatype})
            PMonthSumTrendy.append(pd.merge(DateRef, MonthSumTrendy, on=['MonthDate'], how = 'left'))   
            PMonthMeansTrendy.append(pd.merge(DateRef, MonthMeansTrendy, on=['MonthDate'], how = 'left'))   
            AllMMTrendyNBP = pd.merge(AllMMTrendyNBP,MonthMeansTrendy[[Datatype,'MonthDate']],on='MonthDate', how = 'left')
            AllMMTrendyNBP = AllMMTrendyNBP.rename(columns={Datatype:Model+Datatype})
        if True: #set false if using precipitation
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'ra']+MeanCycle[Model+'rh'],column = Model+'rarh')
            maxFlux = max([MeanCycle[Model+'rarh'].max(),-MeanCycle[Model+'gpp'].max()])
            MeanCycle.insert(loc=1,
                            value= MeanCycle[Model+'rarh']/maxFlux,column = Model+'rarhNorm')
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'gpp']/maxFlux,column = Model+'gppNorm')
            MeanCycle.insert(loc=1,value= MeanCycle[Model+'nbp']/maxFlux,column = Model+'nbpNorm')

    if len(TrendyModels) >= 16:        
            #get mean and stdev for all models
        meanModel = AllMMTrendyNBP[BadModelNames('nbp')+GoodModelsNames('nbp')].mean(axis=1)
        stDevModel = AllMMTrendyNBP[BadModelNames('nbp')+GoodModelsNames('nbp')].std(axis=1,ddof = 0)
        minModel = AllMMTrendyNBP[BadModelNames('nbp')+GoodModelsNames('nbp')].min(axis=1)
        maxModel = AllMMTrendyNBP[BadModelNames('nbp')+GoodModelsNames('nbp')].max(axis=1)
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBP', value = minModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBP', value = maxModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBP', value = meanModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBP', value = stDevModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp')].min(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp')].max(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBPbadModel', value = AllMMTrendyNBP[BadModelNames('nbp')].std(axis=1,ddof = 0))
        AllMMTrendyNBP.insert(loc = 1, column = 'minNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp')].min(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'maxNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp')].max(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'stdNBPgoodModel', value = AllMMTrendyNBP[GoodModelsNames('nbp')].std(axis=1,ddof = 0))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('nbp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPlateModel', value = AllMMTrendyNBP[LateRespNames('nbp')].mean(axis=1))
        
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppbadModel', value = AllMMTrendyNBP[BadModelNames('gpp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppgoodModel', value = AllMMTrendyNBP[GoodModelsNames('gpp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangppearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('gpp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meangpplateModel', value = AllMMTrendyNBP[LateRespNames('gpp')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhrabadModel', value = AllMMTrendyNBP[BadModelNames('ra')].mean(axis=1)+AllMMTrendyNBP[BadModelNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhragoodModel', value = AllMMTrendyNBP[GoodModelsNames('ra')].mean(axis=1)+AllMMTrendyNBP[GoodModelsNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhraearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('ra')].mean(axis=1)+AllMMTrendyNBP[EarlyGPPModelNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhralateModel', value = AllMMTrendyNBP[LateRespNames('ra')].mean(axis=1)+AllMMTrendyNBP[LateRespNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhbadModel', value = AllMMTrendyNBP[BadModelNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhgoodModel', value = AllMMTrendyNBP[GoodModelsNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrhlateModel', value = AllMMTrendyNBP[LateRespNames('rh')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanrabadModel', value = AllMMTrendyNBP[BadModelNames('ra')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanragoodModel', value = AllMMTrendyNBP[GoodModelsNames('ra')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanraearlyModel', value = AllMMTrendyNBP[EarlyGPPModelNames('ra')].mean(axis=1))
        AllMMTrendyNBP.insert(loc = 1, column = 'meanralateModel', value = AllMMTrendyNBP[LateRespNames('ra')].mean(axis=1))
        
        #get mean amplitude
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPgoodModelSwitch',value = -1*AllMMTrendyNBP.meanNBPgoodModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPbadModelSwitch',value = -1*AllMMTrendyNBP.meanNBPbadModel)
        AllMMTrendyNBP.insert(loc = 1, column = 'meanNBPSwitch',value = -1*AllMMTrendyNBP.meanNBP)
        AllMMTrendyNBP.to_pickle(datapath + "/TRENDY/dataframes/MonthFrames/AllMMTrendyNBP"+RegionName+str(Numm)+"_newforSuppPaper.pkl")   
        goodTRENDYmeanAmplitude, goodTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','Year',2009,2018)[0:2]
        goodTRENDYmeanAmplitude_no09, goodTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','Year',2010,2018)[0:2]
        goodTRENDYmeanAmplitudejj, goodTRENDYmeanAmplitudeStDevjj,goodTRENDYmeanAmplitudeMin,goodTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPgoodModelSwitch','JJYear',2009,2017)
        badTRENDYmeanAmplitude, badTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','Year',2009,2018)[0:2]
        badTRENDYmeanAmplitude_no09, badTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','Year',2010,2018)[0:2]
        badTRENDYmeanAmplitudejj, badTRENDYmeanAmplitudeStDevjj,badTRENDYmeanAmplitudeMin,badTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPbadModelSwitch','JJYear',2009,2017)
        allTRENDYmeanAmplitude, allTRENDYmeanAmplitudeStDev = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','Year',2009,2018)[0:2]
        allTRENDYmeanAmplitude_no09, allTRENDYmeanAmplitudeStDev_no09 = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','Year',2010,2018)[0:2]
        allTRENDYmeanAmplitudejj, allTRENDYmeanAmplitudeStDevjj,allTRENDYmeanAmplitudeMin,allTRENDYmeanAmplitudeMax = getMeanAmplitude(AllMMTrendyNBP,'meanNBPSwitch','JJYear',2009,2017)
        
        TRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBP'].sum().mean()
        TRENDYannualStDev = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBP'].sum().std(ddof = 0)
        TRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp')+BadModelNames('nbp')].sum().mean(axis=0).min()
        TRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp')+BadModelNames('nbp')].sum().mean(axis=0).max()
        
        goodTRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPgoodModel'].sum().mean()
        badTRENDYannualMean = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPbadModel'].sum().mean()
        goodTRENDYannualStd = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPgoodModel'].sum().std(ddof = 0)
        badTRENDYannualStd = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])['meanNBPbadModel'].sum().std(ddof = 0)
        goodTRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp')].sum().mean(axis=0).min()
        
        badTRENDYannualMin = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[BadModelNames('nbp')].sum().mean(axis=0).min()
        goodTRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[GoodModelsNames('nbp')].sum().mean(axis=0).max()
        badTRENDYannualMax = AllMMTrendyNBP[(AllMMTrendyNBP.Year >= minyearmean)].groupby(['Year'])[BadModelNames('nbp')].sum().mean(axis=0).max()
        
        dfIAVTrendy = pd.DataFrame(data = {'Dataset':['allTrendy','goodTrendy','badTrendy'],
                                'AnnualMean':[-TRENDYannualMean,-goodTRENDYannualMean,-badTRENDYannualMean],
                                'StDev':[TRENDYannualStDev,goodTRENDYannualStd,badTRENDYannualStd],
                                'AnnualMin':[-TRENDYannualMax,-goodTRENDYannualMax,-badTRENDYannualMax],
                                'AnnualMax':[-TRENDYannualMin,-goodTRENDYannualMin,-badTRENDYannualMin]})
                                
        dfIAV = dfIAV.append(dfIAVTrendy,ignore_index=True)



        MeanGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanAll = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        MeanLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StErrorLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevGood = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevBad = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevAll = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevEarlyGPP = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        StDevLateResp = pd.DataFrame(data = {'Month':np.array(range(1,13))})
        
    if len(TrendyModels) >= 16:
        for Datatype in DataList+['rarh','rarhNorm','gppNorm','nbpNorm']:
            GoodModel = MeanCycle[GoodModelsNames(Datatype)].copy()
            AllModel = MeanCycle[GoodModelsNames(Datatype)+BadModelNames(Datatype)].copy()
            EarlyGPPModel = MeanCycle[EarlyGPPModelNames(Datatype)].copy()
            LateResp = MeanCycle[LateRespNames(Datatype)].copy()
            BadModel = MeanCycle[BadModelNames(Datatype)].copy()
            MeanGood.insert(loc=1,
                            value=GoodModel.mean(axis=1), column = Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.mean(axis=1), column = Datatype)
            MeanAll.insert(loc=1,
                            value=AllModel.mean(axis=1), column = Datatype)
            MeanGood.insert(loc=1,
                            value=GoodModel.min(axis=1), column = 'min_'+Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.min(axis=1), column = 'min_'+Datatype)
            MeanGood.insert(loc=1,
                            value=GoodModel.max(axis=1), column = 'max_'+Datatype)
            MeanBad.insert(loc=1,
                            value=BadModel.max(axis=1), column = 'max_'+Datatype)
            MeanEarlyGPP.insert(loc=1, value=EarlyGPPModel.mean(axis=1),
                                column = Datatype)
            MeanLateResp.insert(loc=1, value=LateResp.mean(axis=1),
                                column = Datatype)
            StErrorGood.insert(loc=1,value=GoodModel.sem(axis=1,ddof = 0), 
                               column = Datatype)
            StErrorBad.insert(loc=1,value=BadModel.sem(axis=1,ddof = 0), 
                               column = Datatype)
            StErrorEarlyGPP.insert(loc=1, value=EarlyGPPModel.sem(axis=1,ddof = 0),
                                column = Datatype)
            StErrorLateResp.insert(loc=1, value=LateResp.sem(axis=1,ddof = 0),
                                column = Datatype)
            StDevGood.insert(loc=1,value=GoodModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevBad.insert(loc=1,value=BadModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevAll.insert(loc=1,value=AllModel.std(axis=1,ddof = 0), 
                               column = Datatype)
            StDevEarlyGPP.insert(loc=1, value=EarlyGPPModel.std(axis=1,ddof = 0),
                                column = Datatype)
            StDevLateResp.insert(loc=1, value=LateResp.std(axis=1,ddof = 0),
                                column = Datatype)
            

if TM5Flux:
    TM5F_region = []
    #rnum = [4,18,20,16,3] 
    rnum = [4,18,21,16,3]
    
    MeanCycleTM5 = pd.DataFrame(data = {'Month':np.array(range(1,13))})
    for nam in TM5list:
        if 'wAllU' in nam:
            data = pd.read_pickle(datapath + "/TM5Inversion/glb3x2_20220413/vpp_files/dataframes_flux/DF2"+nam+".pkl")        
        elif '238' in nam:
            data = pd.read_pickle(datapath + "/TM5InversionRT238/dataframes_flux/DF1wU_"+nam+".pkl")        
        else:
            data = pd.read_pickle(datapath + "/TM5Inversion/dataframes_flux/DF2_"+nam+".pkl")        
        
        if Numm == 949 or Numm == 950 or Numm == 951:
            datadf = data[(data.Region == rnum[902 - 900])]
        else:
            datadf = data[(data.Region == rnum[Numm - 900])]
        datadf = datadf.rename(columns={'Date':'MonthDate'})
        SpM = datadf.apply(lambda x: 60*60*24*getNumDayOfMonth(x.Year,x.Month),axis=1)
        datadf.insert(loc = 1, column = 'SecInMonth', value=SpM)
        datadf.to_pickle(datapath + "/TM5Inversion/dataframes_flux/MonthFrames/MonthMeans_"+nam+"_regional_"+RegionName+str(Numm)+".pkl")
        TM5F_region.append(pd.merge(DateRef, datadf, on=['MonthDate'], how = 'left'))
        datadf.insert(loc = 1, column = 'Flux_fire_mug_m2s', value=datadf.Flux_fire*10**18/(RegionArea * datadf.SecInMonth))
        datadf.insert(loc = 1, column = 'Flux_bio_mug_m2s', value=datadf.Flux_bio*10**18/(RegionArea * datadf.SecInMonth))
        MeanCycleTM5 = MeanCycleTM5.merge(datadf.groupby(['Month'])[['Flux_fire_mug_m2s','Flux_bio_mug_m2s']].mean(), on = 'Month')
        MeanCycleTM5 = MeanCycleTM5.rename(columns = {'Flux_fire_mug_m2s': nam+'Flux_fire_mug_m2s','Flux_bio_mug_m2s': nam+'Flux_bio_mug_m2s'})
            

#calculate mean seasonal cycle and key numbers
if TM5Flux:#1x1:#mean of the satellite based TM5
    #1x1 flux
    '''
    allSatModel = TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_ACOS+IS')[0][0]][['MonthDate','JJYear','CO2_flux_nee_total','CO2_flux_fire_total']]
    allSatModel.insert(loc = 1, column = 'TM5ACOSIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
    allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
    allSatModel = pd.merge(allSatModel,TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_RemoTeC_2.4.0+IS')[0][0]][['MonthDate','CO2_flux_nee_total','CO2_flux_fire_total']])
    allSatModel.insert(loc = 1, column = 'TM5RTIS_nbp', value = (allSatModel.CO2_flux_nee_total+allSatModel.CO2_flux_fire_total)*10**(-12))
    allSatModel = allSatModel.drop(columns = ['CO2_flux_nee_total','CO2_flux_fire_total'])
    #regional fluxes
    '''
    
    allSatModel = TM5F_region[0][['MonthDate','JJYear','Flux_bio','Flux_fire']].copy()
    allSatModel.insert(loc = 1, column = 'TM5RTIS_nbp', value = allSatModel.Flux_bio+allSatModel.Flux_fire)
    allSatModel = allSatModel.drop(columns = ['Flux_bio','Flux_fire'])
    allSatModel = pd.merge(allSatModel,TM5F_region[2][['MonthDate','Flux_bio','Flux_fire']])
    allSatModel.insert(loc = 1, column = 'TM5ACOSIS_nbp', value = allSatModel.Flux_bio+allSatModel.Flux_fire)
    allSatModel = allSatModel.drop(columns = ['Flux_bio','Flux_fire'])
    
    ParamList = ['TM5RTIS_nbp','TM5ACOSIS_nbp']
    
    if TM5Flux1x1 and WithOCO2:
        print(np.where(np.array(TM51x1dataset) == 'flux_1x1_LNLGIS')[0][0])
        print(len(TM51x1dataset))
        allSatModel = pd.merge(allSatModel,TM51x1_regionl[np.where(np.array(TM51x1dataset) == 'flux_1x1_LNLGIS')[0][0]][['MonthDate','CO2_flux_nbe_total']])
        allSatModel.insert(loc = 1, column = 'TM5OCO2IS_nbp', value = allSatModel.CO2_flux_nbe_total*10**(-12))
        print(allSatModel.TM5OCO2IS_nbp.max())
        print(allSatModel.TM5RTIS_nbp.max())
        print(allSatModel.TM5ACOSIS_nbp.max())
        allSatModel = allSatModel.drop(columns = ['CO2_flux_nbe_total'])
        ParamList = ['TM5RTIS_nbp','TM5ACOSIS_nbp','TM5OCO2IS_nbp']
    allSatModel.insert(loc = 1, column = 'mean_nbp', value = allSatModel[ParamList].mean(axis = 1))
    allSatModel.insert(loc = 1, column = 'std_nbp', value = allSatModel[ParamList].std(axis = 1,ddof = 0))
    allSatModel.insert(loc = 1, column = 'min_nbp', value = allSatModel[ParamList].min(axis = 1))
    allSatModel.insert(loc = 1, column = 'max_nbp', value = allSatModel[ParamList].max(axis = 1))
    Monthlist = allSatModel.apply(lambda x: x.MonthDate.month,axis=1)
    allSatModel.insert(loc = 1, column = 'Month', value = Monthlist)
    Yearlist = allSatModel.apply(lambda x: x.MonthDate.year,axis=1)
    allSatModel.insert(loc = 1, column = 'Year', value = Yearlist) 
    allSatModel.to_pickle(datapath + "/TM5Inversion/dataframes_flux/MonthFrames/MonthMeansAllSatModels_newForSuppPaper_"+RegionName+str(Numm)+".pkl")          
    AMeanSatModel = allSatModel.groupby(['Month'])['mean_nbp'].mean().reset_index()
    
    #get mean amplitude
    allSatModelYmeanAmplitude, allSatModelmeanAmplitudeStDev = getMeanAmplitude(allSatModel,'mean_nbp','Year',2009,2018)[0:2]
    allSatModelYmeanAmplitude_no09, allSatModelmeanAmplitudeStDev_no09 = getMeanAmplitude(allSatModel,'mean_nbp','Year',2010,2018)[0:2]
    allSatModelYmeanAmplitudejj, allSatModelmeanAmplitudeStDevjj,allSatModelmeanAmplitudeMin,allSatModelmeanAmplitudeMax = getMeanAmplitude(allSatModel,'mean_nbp','JJYear',2009,2017)
    
    allSatModelannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().mean()
    allSatModelannualStDev = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().std(ddof = 0)
    allSatModelMeanSais = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].mean().reset_index()
    allSatModelMeanSaisStDev = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].std(ddof = 0).reset_index()
    
    RTannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5RTIS_nbp'].sum().mean()
    ACOSannualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5ACOSIS_nbp'].sum().mean()
    if WithOCO2:
        OCO2annualMean = allSatModel[(allSatModel.Year >= minyearmean)].groupby(['Year'])['TM5OCO2IS_nbp'].sum().mean()
    
    if WithOCO2:
        satmeans = [ACOSannualMean,RTannualMean,OCO2annualMean]
    else:
        satmeans = [ACOSannualMean,RTannualMean]
    dfIAVGOSAT = pd.DataFrame(data = {'Dataset':['GOSAT'],
                            'AnnualMean':[allSatModelannualMean],
                            'StDev':[allSatModelannualStDev],
                            'AnnualMin':[min(satmeans)],
                            'AnnualMax':[max(satmeans)]})
    dfIAV = dfIAV.append(dfIAVGOSAT,ignore_index=True)

if CAMS_b and CT_flux and TM5Flux:#mean of the inverse models
    if np.any(np.where('ISglb6x4' == np.array(TM5list))):
        ISnum = np.where('ISglb6x4' == np.array(TM5list))
    else:
        ISnum = np.where('ISglb3x2' == np.array(TM5list))
    allModel = TM5F_region[ISnum[0][0]][['MonthDate','JJYear','Flux_bio','Flux_fire']].copy()
    allModel.insert(loc = 1, column = 'TM5_nbp', value = allModel.Flux_bio+allModel.Flux_fire)
    allModel = allModel.drop(columns = ['Flux_bio','Flux_fire'])
    allModel = pd.merge(allModel,PCTF_region[['MonthDate','Biotot','firetot']])
    allModel.insert(loc = 1, column = 'CT_nbp', value = allModel.Biotot*10**-12+allModel.firetot*10**-12)
    allModel = allModel.drop(columns = ['Biotot','firetot'])
    allModel = pd.merge(allModel,PCAMS_region_Sur[['MonthDate','flux_apos_bio_tot']])
    allModel.insert(loc = 1, column = 'CAMSsur_nbp', value = allModel.flux_apos_bio_tot*10**-12)
    allModel = allModel.drop(columns = ['flux_apos_bio_tot'])
    allModel.insert(loc = 1, column = 'mean_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].mean(axis = 1))
    allModel.insert(loc = 1, column = 'std_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].std(axis = 1,ddof = 0))
    allModel.insert(loc = 1, column = 'min_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].min(axis = 1))
    allModel.insert(loc = 1, column = 'max_nbp', value = allModel[['CAMSsur_nbp','CT_nbp','TM5_nbp']].max(axis = 1))
    Monthlist = allModel.apply(lambda x: x.MonthDate.month,axis=1)
    allModel.insert(loc = 1, column = 'Month', value = Monthlist)
    Yearlist = allModel.apply(lambda x: x.MonthDate.year,axis=1)
    allModel.insert(loc = 1, column = 'Year', value = Yearlist)     
    allModel.to_pickle(datapath + "/CT2019/dataframes/monthFrames/MonthMeansAllISModels_"+RegionName+str(Numm)+"_newForSuppPaper_.pkl")
        
    #get mean Amplitude
    allModelYmeanAmplitude, allModelmeanAmplitudeStDev = getMeanAmplitude(allModel,'mean_nbp','Year',2009,2018)[0:2]
    allModelYmeanAmplitude_no09, allModelmeanAmplitudeStDev_no09 = getMeanAmplitude(allModel,'mean_nbp','Year',2010,2018)[0:2]
    allModelYmeanAmplitudejj, allModelmeanAmplitudeStDevjj,allModelmeanAmplitudeMin,allModelmeanAmplitudeMax = getMeanAmplitude(allModel,'mean_nbp','JJYear',2009,2017)
    
    allModelannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().mean()
    allModelannualStDev = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['mean_nbp'].sum().std(ddof = 0)
    allModelMeanSais = allModel[(allModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].mean().reset_index()
    allModelMeanSaisStDev = allModel[(allModel.Year >= minyearmean)].groupby(['Month'])['mean_nbp'].std(ddof = 0).reset_index()
    
    CTannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['CT_nbp'].sum().mean()
    CAMSannualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['CAMSsur_nbp'].sum().mean()
    TM5annualMean = allModel[(allModel.Year >= minyearmean)].groupby(['Year'])['TM5_nbp'].sum().mean()
    
    dfIAVInvM = pd.DataFrame(data = {'Dataset':['InvModels'],
                            'AnnualMean':[allModelannualMean],
                            'StDev':[allModelannualStDev],
                            'AnnualMin':[CAMSannualMean],
                            'AnnualMax':[TM5annualMean]})
    dfIAV = dfIAV.append(dfIAVInvM,ignore_index=True)
    #dfIAV.to_csv('/home/eschoema/Results/TableIAV_withoutOCO2_without2009_regsatflux.csv')



#Plot Data

#Figure 2 Paper
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(2, 2, 
                            figsize = (18.0*cm,12.2*cm),
                            gridspec_kw={'width_ratios': [4, 1],'height_ratios':[1,1]})

    #get TRENDY std
    StdCycnbp = AllMMTrendyNBP.groupby('Month')['meanNBP'].std(ddof=0).reset_index()
    StdCycgoodnbp = AllMMTrendyNBP.groupby('Month')['meanNBPgoodModel'].std(ddof=0).reset_index()
    StdCycbadnbp = AllMMTrendyNBP.groupby('Month')['meanNBPbadModel'].std(ddof=0).reset_index()
        

    # Panel a)
    ax2[0,0].text(0.03, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0,0].transAxes,fontsize=10, weight='bold')
    #Shading
    ## stdDev for TRENDY
    ax2[0,0].fill_between(AllMMTrendyNBP.MonthDate,
                        -(AllMMTrendyNBP.meanNBP-AllMMTrendyNBP.stdNBP)*factor,
                        -(AllMMTrendyNBP.meanNBP+AllMMTrendyNBP.stdNBP)*factor,
                        color= 'grey', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    ##Min Max for Sat and Models
    ax2[0,0].fill_between(allModel.MonthDate,allModel.max_nbp*factor,
                          allModel.min_nbp*factor
                          ,color= 'blue', alpha=0.3,zorder = 2)
    ax2[0,0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    #monthly means
    ax2[0,0].plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[0,0].plot(allModel.MonthDate,allModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm Inverse~model_{in-situ}$',zorder = 4)
    ax2[0,0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2[0,0].plot(PSIFFC_region[0].MonthDate,(PSIFFC_region[0].NEEtot/10**12*(44/12)+PGFEDregion.total_emission)*factor,ls='-',linewidth = lw,marker= '',color= 'gold',label='FLUXCOM+GFED',zorder = 4)
    ax2[0,0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBP*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label=r'$\rm TRENDY_{all}$',zorder = 4)
       
    # Panel c)
    ax2[1,0].text(0.03, 0.9, r'$\rm \bf{C}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1,0].transAxes,fontsize=10, weight='bold')
    # SHADING
    ## stDev for TRENDY
    ax2[1,0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPbadModel-AllMMTrendyNBP.stdNBPbadModel)*factor,-(AllMMTrendyNBP.meanNBPbadModel+AllMMTrendyNBP.stdNBPbadModel)*factor,color= 'grey', alpha=0.3,zorder = 2)#,label='bad TRENDY')#,label = 'StDev')
    ax2[1,0].fill_between(AllMMTrendyNBP.MonthDate,-(AllMMTrendyNBP.meanNBPgoodModel-AllMMTrendyNBP.stdNBPgoodModel)*factor,-(AllMMTrendyNBP.meanNBPgoodModel+AllMMTrendyNBP.stdNBPgoodModel)*factor,color= 'black', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    ## min max for satellite
    ax2[1,0].fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)     
    #monthly means
    ax2[1,0].plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1,0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBPgoodModel*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{selection}$',zorder = 4)
    ax2[1,0].plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBPbadModel*factor,ls='-',linewidth = lw,marker= '',color= 'dimgrey',label=r'$\rm TRENDY_{others}$',zorder = 4)
    ax2[1,0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2[1,0].plot(PSIFFC_region[0].MonthDate,(PGFEDregion.total_emission)*factor,ls='-',linewidth = lw,marker= '',color= 'orange',label=r'$\rm GFED~fire~CO_2$',zorder = 4)
    
    # Panel b)
    indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]

    ax2[0,1].text(0.15, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0,1].transAxes,fontsize=10, weight='bold')
    ax2[0,1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[0,1].plot(range(1,13),-MeanAll.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ax2[0,1].plot(range(1,13),(FCMonthlyAnnualMean.nbp.iloc[indexrange])*factor,ls='-',linewidth = lw,marker= '',color= 'gold',label='FLUXCOM+GFED',zorder = 4)
    ax2[0,1].plot(range(1,13),allModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm Inverse~model_{in-situ}$',zorder = 4)
    ax2[0,1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    # SHADING
    ax2[0,1].fill_between(range(1,13),(allModelMeanSais.mean_nbp.iloc[indexrange]-allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allModelMeanSais.mean_nbp.iloc[indexrange]+allModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'blue', alpha=0.3,zorder = 2)
    ax2[0,1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    #ax2[0,1].fill_between(range(1,13),-(MeanAll.nbp.iloc[indexrange]-StDevAll.nbp.iloc[indexrange])*factor,-(MeanAll.nbp.iloc[indexrange]+StDevAll.nbp.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[0,1].fill_between(range(1,13),-(MeanAll.nbp.iloc[indexrange]-StdCycnbp.meanNBP.iloc[indexrange])*factor,-(MeanAll.nbp.iloc[indexrange]+StdCycnbp.meanNBP.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[0,1].fill_between(range(1,13),(FCMonthlyAnnualMean.nbp.iloc[indexrange]-FCMonthlyAnnualStDev.nbp.iloc[indexrange])*factor,(FCMonthlyAnnualMean.nbp.iloc[indexrange]+FCMonthlyAnnualStDev.nbp.iloc[indexrange])*factor,color= 'gold', alpha=0.3,zorder = 2)
    # ERRORBAR
    ax2[0,1].errorbar([14.5],0,np.array([-allTRENDYmeanAmplitudeMin*factor,allTRENDYmeanAmplitudeMax*factor]).reshape(2,1)
                ,ls='-',linewidth = lw,marker= '',color= 'grey',ecolor = 'grey',capsize = cs,zorder = 2)
    ax2[0,1].errorbar([14],0,np.array([-FCGFEDmeanAmplitudeMin*factor,FCGFEDmeanAmplitudeMax*factor]).reshape(2,1)
                ,ls='-',linewidth = lw,marker= '',color= 'gold',ecolor = 'gold',capsize = cs,zorder = 2)
    ax2[0,1].errorbar([13.5],0,np.array([-allModelmeanAmplitudeMin*factor,allModelmeanAmplitudeMax*factor]).reshape(2,1)
                ,ls='-',linewidth = lw,marker= '',color= 'blue',ecolor = 'blue',capsize = cs,zorder = 2)
    ax2[0,1].errorbar([13],0,np.array([-allSatModelmeanAmplitudeMin*factor,allSatModelmeanAmplitudeMax*factor]).reshape(2,1)
                ,ls='-',linewidth = lw,marker= '',color= 'darkred',ecolor ='darkred',capsize = cs,zorder = 2)
    print(np.array([-allSatModelmeanAmplitudeMin*factor,allSatModelmeanAmplitudeMax*factor]).reshape(2,1))
    #Panel d)
    ax2[1,1].text(0.15, 0.9, r'$\rm \bf{D}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1,1].transAxes,fontsize=10, weight='bold')
    
    ax2[1,1].plot([0.5,15.5],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2[1,1].plot(range(1,13),-MeanGood.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ax2[1,1].plot(range(1,13),-MeanBad.nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ax2[1,1].plot(range(1,13),allSatModelMeanSais.mean_nbp.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2[1,1].plot(range(1,13),GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]*factor,ls='-',linewidth = lw,marker= '',color= 'orange',label=r'$\rm GFED~fire~CO_{2}$',zorder = 4)    
    # SHADING
    ax2[1,1].fill_between(range(1,13),(allSatModelMeanSais.mean_nbp.iloc[indexrange]-allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,(allSatModelMeanSais.mean_nbp.iloc[indexrange]+allSatModelMeanSaisStDev.mean_nbp.iloc[indexrange])*factor,color= 'red', alpha=0.3,zorder = 2)
    #ax2[1,1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StDevGood.nbp.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StDevGood.nbp.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    #ax2[1,1].fill_between(range(1,13),-(MeanBad.nbp.iloc[indexrange]-StDevBad.nbp.iloc[indexrange])*factor,-(MeanBad.nbp.iloc[indexrange]+StDevBad.nbp.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[1,1].fill_between(range(1,13),-(MeanGood.nbp.iloc[indexrange]-StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,-(MeanGood.nbp.iloc[indexrange]+StdCycgoodnbp.meanNBPgoodModel.iloc[indexrange])*factor,color= 'black', alpha=0.3,zorder = 2)
    ax2[1,1].fill_between(range(1,13),-(MeanBad.nbp.iloc[indexrange]-StdCycbadnbp.meanNBPbadModel.iloc[indexrange])*factor,-(MeanBad.nbp.iloc[indexrange]+StdCycbadnbp.meanNBPbadModel.iloc[indexrange])*factor,color= 'grey', alpha=0.3,zorder = 2)
    ax2[1,1].fill_between(range(1,13),(GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]-GFEDMonthlyAnnualStd.total_emission.iloc[indexrange])*factor,(GFEDMonthlyAnnualMean.total_emission.iloc[indexrange]+GFEDMonthlyAnnualStd.total_emission.iloc[indexrange])*factor,color= 'orange', alpha=0.3,zorder = 2)
    # ERRORBAR
    ax2[1,1].errorbar([13],0,np.array([-allSatModelmeanAmplitudeMin*factor,allSatModelmeanAmplitudeMax*factor]).reshape(2,1)
                ,ls='-',linewidth = lw,marker= '',color= 'darkred',ecolor ='darkred',capsize = cs,zorder = 2)
    ax2[1,1].errorbar([13.5],0,np.array([-goodTRENDYmeanAmplitudeMin*factor,goodTRENDYmeanAmplitudeMax*factor]).reshape(2,1)
                   ,ls='-',linewidth = lw,marker= '',color= 'black',ecolor= 'black',capsize = cs,zorder = 2)
    ax2[1,1].errorbar([14],0,np.array([-badTRENDYmeanAmplitudeMin*factor,badTRENDYmeanAmplitudeMax*factor]).reshape(2,1)
                   ,ls='-',linewidth = lw,marker= '',color= 'grey',ecolor= 'grey',capsize = cs,zorder = 2)
    ax2[1,1].errorbar([14.5],0,np.array([GFEDmeanAmplitudeMin*factor,GFEDmeanAmplitudeMax*factor]).reshape(2,1)
                   ,ls='-',linewidth = lw,marker= '',color= 'orange',ecolor= 'orange',capsize = cs,zorder = 2)
    
    # SETTINGS
    
    #Panel a)
    ax2[0,0].legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[0,0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[0,0].set_ylim(-280,280)
    ax2[0,0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[0,0].set_xticklabels([])
    #ax2[0,0].set_yticklabels(ax2[0,0].get_yticks(),fontsize=7)
    ax2[0,0].tick_params(axis = 'x',length = 0.01, zorder = 0)
    ax2[0,0].tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2[0,0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[0,0].set_axisbelow(True)

    #Panel c)
    ax2[1,0].legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2[1,0].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[1,0].set_xlabel(r'Date',fontsize=7)
    ax2[1,0].set_ylim(-280,280)
    #ax2[1,0].set_xticklabels(ax2[1,0].get_xticks(),fontsize=7)
    #ax2[1,0].set_yticklabels(ax2[1,0].get_yticks(),fontsize=7)
    ax2[1,0].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2[1,0].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[1,0].set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2[1,0].set_axisbelow(True)

    #Panel b)  
    secaxB = ax2[0,1].twinx()     
    secaxB.set_ylim(-280,280)
    secaxB.set_yticklabels([])
    secaxB.tick_params(axis = 'y',direction = 'in')          
    ax2[0,1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    ax2[0,1].set_ylim(-280,280)
    ax2[0,1].set_xticklabels([])
    ax2[0,1].set_yticklabels([])
    ax2[0,1].tick_params(axis = 'both',length = 0.01, zorder = 0)
    ax2[0,1].set_xlim(0.5,15.5)
    ax2[0,1].set_axisbelow(True)
    #Panel d)       
    secaxD = ax2[1,1].twinx()     
    secaxD.set_ylim(-280,280)
    secaxD.set_yticklabels([])
    secaxD.tick_params(axis = 'y',direction = 'in')
    ax2[1,1].set_xticks(ticks = [3,6,9,12])
    ax2[1,1].set_xticklabels(['9','12','3','6'])
    ax2[1,1].axvline(6.5, linestyle='-', color='darkgrey',linewidth = 0.75,zorder = 0) # vertical lines
    #ax2[1,1].grid(True,which = 'major', axis='x',zorder = 0) 
    ax2[1,1].set_xlabel(r'Month',fontsize=7)
    ax2[1,1].set_ylim(-280,280)
    ax2[1,1].set_yticklabels([])
    ax2[1,1].tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2[1,1].tick_params(axis = 'x',labelsize = 6)
    ax2[1,1].set_xlim(0.5,15.5)
    ax2[1,1].set_axisbelow(True)
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig2complete_3x2res_corrSTdTrendy_Region43TM5"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig2complete_3x2res_corrSTdTrendy_Region43TM5"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')
    
    
#Figure new S7 Paper OCO-2 MIP flux 
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (14.4*cm,6.1*cm))
    # Panel a)
    #Shading
    ## stdDev for TRENDY
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        -(AllMMTrendyNBP.meanNBP-AllMMTrendyNBP.stdNBP)*factor,
                        -(AllMMTrendyNBP.meanNBP+AllMMTrendyNBP.stdNBP)*factor,
                        color= 'grey', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    ## stdDev for MIP
    ax2.fill_between(PMIP.MonthDate,
                        (PMIP.Landtot+PMIP.LandStdtot2)*(44/12)*factor,
                        (PMIP.Landtot-PMIP.LandStdtot2)*(44/12)*factor,
                        color= 'black', alpha=0.3,zorder = 2)#,label='good TRENDY')#,label = 'StDev')
    
    ##Min Max for Sat and Models
    ax2.fill_between(allModel.MonthDate,allModel.max_nbp*factor,
                          allModel.min_nbp*factor
                          ,color= 'blue', alpha=0.3,zorder = 2)
    ax2.fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    #monthly means
    ax2.plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    ax2.plot(allModel.MonthDate,allModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm Inverse~model_{in-situ}$',zorder = 4)
    ax2.plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2.plot(PSIFFC_region[0].MonthDate,(PSIFFC_region[0].NEEtot/10**12*(44/12)+PGFEDregion.total_emission)*factor,ls='-',linewidth = lw,marker= '',color= 'gold',label='FLUXCOM+GFED',zorder = 4)
    ax2.plot(AllMMTrendyNBP.MonthDate,-AllMMTrendyNBP.meanNBP*factor,ls='-',linewidth = lw,marker= '',color= 'grey',label=r'$\rm TRENDY_{all}$',zorder = 4)
    ax2.plot(PMIP.MonthDate,PMIP.Landtot*(44/12)*factor,ls='-',linewidth = lw,marker= '',color= 'black',label=r'$\rm Inverse~model_{+OCO-2}$',zorder = 4)
    # SETTINGS
    #Panel a)
    ax2.legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_ylim(-280,280)
    ax2.set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    #ax2.set_xticklabels([])
    #ax2[0,0].set_yticklabels(ax2[0,0].get_yticks(),fontsize=7)
    #ax2[0,0].tick_params(axis = 'x',length = 0.01, zorder = 0)
    #ax2[0,0].tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2.tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.set_axisbelow(True)

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/FigS4_MIP_3_Region43TM5_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/FigS4_MIP_3_Region43TM5_"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')



#Figure S5 Apriori
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (14.4*cm,6.1*cm))
    
    ax2.plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',color= 'grey',zorder = 0)
    #Shading
    ax2.fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.3,zorder = 2)  
    #monthly means
    ax2.plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$',zorder = 4)
    ax2.plot(TM5F_region[1].MonthDate,(TM5F_region[1].Flux_bio+TM5F_region[1].Flux_fire)*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r'$\rm TM5-4DVar_{in-situ}$')
    ax2.plot(TM5F_region[2].MonthDate,(TM5F_region[2].Flux_bio_apri+TM5F_region[2].Flux_fire_apri)*factor,ls=':',linewidth = lw,marker= '',color= 'black',label=r'$\rm TM5-4DVar_{apri}$')
    
    # SETTINGS
    ax2.legend(fontsize=6,loc = 4,ncol = 2)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    ax2.set_ylim(-280,280)
    ax2.set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    #ax2.set_xticklabels([])
    #ax2[0,0].set_yticklabels(ax2[0,0].get_yticks(),fontsize=7)
    #ax2[0,0].tick_params(axis = 'x',length = 0.01, zorder = 0)
    #ax2[0,0].tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2.tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.set_axisbelow(True)

    plt.savefig(savepath + "/Results/Plots/Paper2021/final/FigS0_Apriori_Apost_3x2_Region43TM5"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/FigS0_Apriori_Apost_3x2_Region43TM5"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')
    

#Figure 3 Paper
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    mpl.rcParams['hatch.linewidth'] = 0.5
    fig2, ax2 = plt.subplots(4, 1, 
                            figsize = (18.0*cm,15.3*cm),
                            gridspec_kw={'width_ratios': [1],'height_ratios':[2,2,0.4,1]})
    # TRENDY, Panel a) and b)
    for numTR in [0,1]:
        clusterlist = ['good','bad']
        clusterkind = clusterlist[numTR]
        ax2[numTR].plot([AllMMTrendyNBP.MonthDate.values[0],AllMMTrendyNBP.MonthDate.values[-1]],[0,0],color = 'grey',ls='-', marker='')
        y_gpp = AllMMTrendyNBP['meangpp'+clusterkind+'Model']*12/44
        y_rarh = (AllMMTrendyNBP['meanrhra'+clusterkind+'Model'])*12/44
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp,linewidth = 0.7,
                            ls='-',facecolor = 'none',edgecolor= 'green',hatch='//',
                            label = 'GPP')
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_rarh,linewidth = 0.7,
                            ls='-',facecolor = 'none',edgecolor= 'indigo',hatch='\\\\',
                            label = 'TER')
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp, y_rarh,linewidth = 0.7,
                            ls='-',where =(y_gpp>y_rarh),color= 'green',
                            alpha = 0.4,interpolate=True)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp,y_rarh,linewidth = 0.7,
                            ls='-',where =(y_gpp<y_rarh),color= 'indigo',
                            alpha = 0.6,interpolate=True)
        yshading = np.array([850,800,700,600,500,425])*12/44
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[0],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[1],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[2],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[3],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[4],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[5],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        
        y1 = -AllMMTrendyNBP['meanNBP'+clusterkind+'Model']*12/44#*10**9*44/12
        yStd = AllMMTrendyNBP['stdNBP'+clusterkind+'Model']*12/44
        ax2[numTR].plot(AllMMTrendyNBP.MonthDate,y1,ls='-',marker= '',color= 'black',label='NBP',linewidth = 0.7)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1>0),color= 'indigo',#'cornflowerblue'
                            alpha = 0.6,interpolate=True)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1<0),color= 'green'
                            ,alpha = 0.4,interpolate=True)
                    
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            (y1-yStd),#*10**9*44/12,
                            (y1+yStd),#*10**9*44/12,
                            ls='-',color= 'grey',alpha = 0.4)
        ax2[numTR].grid(True,which = 'both', axis='x')
        ax2[numTR].set_ylim(-300,800)
        ax2[numTR].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    

    ax2[1].set_xlabel(r'Date',fontsize=7)
    ax2[0].plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='--',linewidth = 0.7,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT}$')       
    ax2[0].set_xticklabels([])
    ax2[0].tick_params(axis = 'x', length = 0.01, zorder = 0, labelsize = 6)
    ax2[0].tick_params(axis = 'y', labelsize = 6, zorder = 0)
    ax2[1].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[0].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[1].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[0].text(0.02, 0.93, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=9, fontweight='bold')
    ax2[0].text(0.095, 0.93, 'Selection', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=9)
    ax2[1].text(0.02, 0.93, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=9, fontweight='bold')
    ax2[1].text(0.075, 0.93, 'Other', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=9)
    ax2[0].set_axisbelow(True)
    ax2[1].set_axisbelow(True)
    #empty plot to vary subplot spacing
    ax2[2].set_visible(False)        

    #ERA 5
    ax2[3].bar(MonthSumERA5.MonthDate,MonthSumERA5.monthlytp*1000,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'White',label='Total area')
    ax2[3].bar(MonthSumERA5dry.MonthDate,MonthSumERA5dry.monthlytp*1000,width = 15,linewidth = 0.7,edgecolor = 'blue',color = 'none',label='Semi-arid area')
           
    ax2[3].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[3].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[3].set_xlabel(r'Date',fontsize=7)
    ax2[3].set_ylabel('Mean monthly \n precipitation \n (mm)',fontsize=7)
    ax2[3].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[3].text(0.021, 0.85, r'$\rm \bf{C}$', horizontalalignment='center', verticalalignment='center', transform=ax2[3].transAxes,fontsize=9, fontweight='bold')
    ax2[3].set_axisbelow(True) #to put grid behind the plots, zorder in grid does not work as grid is part of the axis
    
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    #plt.savefig(savepath + "/Results/Plots/Paper2021/Fig3completeTER_3x2_"+RegionName+str(Numm)+".pdf", dpi=250,bbox_inches='tight',format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/Fig3completeTER_3x2_newcolor"+RegionName+str(Numm)+".png", dpi=250,bbox_inches='tight')
    #plt.savefig(savepath + "/Results/Plots/Paper2021/Fig3_withra_rh_"+RegionName+str(Numm)+".png", dpi=250,bbox_inches='tight')#,format = 'pdf')

#Figure 3 without fires TRENDY NBP -> GPP - resp, GOSAT - GFED
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    mpl.rcParams['hatch.linewidth'] = 0.5
    fig2, ax2 = plt.subplots(4, 1, 
                            figsize = (18.0*cm,15.3*cm),
                            gridspec_kw={'width_ratios': [1],'height_ratios':[2,2,0.4,1]})
    # TRENDY, Panel a) and b)
    for numTR in [0,1]:
        clusterlist = ['good','bad']
        clusterkind = clusterlist[numTR]
        #zero line
        ax2[numTR].plot([AllMMTrendyNBP.MonthDate.values[0],AllMMTrendyNBP.MonthDate.values[-1]],[0,0],color = 'grey',ls='-', marker='')
        y_gpp = AllMMTrendyNBP['meangpp'+clusterkind+'Model']*12/44
        y_rarh = (AllMMTrendyNBP['meanrhra'+clusterkind+'Model'])*12/44
        #hatch
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp,linewidth = 0.7,
                            ls='-',facecolor = 'none',edgecolor= 'green',hatch='//',
                            label = 'GPP')
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_rarh,linewidth = 0.7,
                            ls='-',facecolor = 'none',edgecolor= 'indigo',hatch='\\\\',
                            label = 'TER')
        #diff filling gpp and resp
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp, y_rarh,linewidth = 0.7,
                            ls='-',where =(y_gpp>y_rarh),color= 'green',
                            alpha = 0.4,interpolate=True)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,
                            y_gpp,y_rarh,linewidth = 0.7,
                            ls='-',where =(y_gpp<y_rarh),color= 'indigo',
                            alpha = 0.6,interpolate=True)
        #fading                    
        yshading = np.array([850,800,700,600,500,425])*12/44
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[0],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[1],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[2],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[3],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[4],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate, yshading[5],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
        
        y1 = y_rarh-y_gpp
        yStd = AllMMTrendyNBP['stdNBP'+clusterkind+'Model']*12/44
        #nbp
        ax2[numTR].plot(AllMMTrendyNBP.MonthDate,y1,ls='-',marker= '',color= 'black',label='TER - GPP',linewidth = 0.7)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1>0),color= 'indigo',
                            alpha = 0.6,interpolate=True)
        ax2[numTR].fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1<0),color= 'green'
                            ,alpha = 0.4,interpolate=True)
                    
        ax2[numTR].grid(True,which = 'both', axis='x')
        ax2[numTR].set_ylim(-300,800)
        ax2[numTR].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    
    
    ax2[1].set_xlabel(r'Date',fontsize=7)
    ax2[0].plot(allSatModel.MonthDate,(allSatModel.mean_nbp-PGFEDregion.total_emission)*factor,ls='--',linewidth = 0.7,marker= '',color= 'darkred',label=r'$\rm Inverse~model_{+GOSAT} - GFED~fire~emissions$')       
    ax2[0].set_xticklabels([])
    ax2[0].tick_params(axis = 'x', length = 0.01, zorder = 0, labelsize = 6)
    ax2[0].tick_params(axis = 'y', labelsize = 6, zorder = 0)
    ax2[1].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[0].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[1].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[0].text(0.02, 0.93, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10, fontweight='bold')
    ax2[0].text(0.095, 0.93, 'Selection', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10)
    ax2[1].text(0.02, 0.93, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10, fontweight='bold')
    ax2[1].text(0.075, 0.93, 'Other', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10)
    ax2[0].set_axisbelow(True)
    ax2[1].set_axisbelow(True)
    #empty plot to vary subplot spacing
    ax2[2].set_visible(False)        

    #ERA 5
    ax2[3].bar(MonthSumERA5.MonthDate,MonthSumERA5.monthlytp*1000,width = 15,linewidth = 0.7,edgecolor = 'Black',color = 'White',label='Total area')
    ax2[3].bar(MonthSumERA5dry.MonthDate,MonthSumERA5dry.monthlytp*1000,width = 15,linewidth = 0.7,edgecolor = 'blue',color = 'none',label='Semi-arid area')
           
    ax2[3].legend(fontsize=6,loc=1,ncol = 3)#,loc=3)#, loc = 9)
    ax2[3].grid(True,which = 'both', axis='x',zorder = 0)
    ax2[3].set_xlabel(r'Date',fontsize=7)
    ax2[3].set_ylabel('Mean monthly \n precipitation \n (mm)',fontsize=7)
    ax2[3].tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2[3].text(0.021, 0.85, r'$\rm \bf{C}$', horizontalalignment='center', verticalalignment='center', transform=ax2[3].transAxes,fontsize=10, fontweight='bold')
    ax2[3].set_axisbelow(True) #to put grid behind the plots, zorder in grid does not work as grid is part of the axis
    
    plt.subplots_adjust(wspace=0,  
                    hspace=0) 
    
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig3_replaceNBP_3x2_newcolor_Region43TM5_"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Fig3_replaceNBP_3x2_newcolor_Region43TM5_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')
    

#Supp Figure 9 Semi-arid, not semi arid 
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    mpl.rcParams['hatch.linewidth'] = 0.5
    fig2, ax2 = plt.subplots(figsize = (18.0*cm,6*cm))

    clusterkind = 'good'

    ax2.plot([AllMMTrendyNBP.MonthDate.values[0],AllMMTrendyNBP.MonthDate.values[-1]],[0,0],color = 'grey',ls='-', marker='')
        
    y_gpp = AllMMTrendyNBP['meangpp'+clusterkind+'Model']*factor
    y_rarh = (AllMMTrendyNBP['meanrhra'+clusterkind+'Model'])*factor
    print(AllMMTrendyNBP.MonthDate)
    
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp,
                        ls='-',linewidth = lw,facecolor = 'none',edgecolor= 'green',hatch='//',
                        label = 'GPP')
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        y_rarh,
                        ls='-',linewidth = lw,facecolor = 'none',edgecolor= 'indigo',hatch='\\\\',
                        label = ' Total respiration')
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp, y_rarh,
                        ls='-',linewidth = lw,where =(y_gpp>y_rarh),color= 'green',
                        alpha = 0.4,interpolate=True)
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        y_gpp,y_rarh,
                        ls='-',linewidth = lw,where =(y_gpp<y_rarh),color= 'indigo',
                        alpha = 0.6,interpolate=True)
    #yshading = np.array([850,800,700,600,500,425])*12/44 #all
    #yshading = [100,90,80,70,60,50] #wet good ind scaling
    yshading = [80,75,70,65,60,50] #dry good ind scaling
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[0],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[1],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[2],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[3],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[4],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    ax2.fill_between(AllMMTrendyNBP.MonthDate, yshading[5],0, edgecolor = 'none',facecolor= 'white', alpha = 0.4)
    
    y1 = -AllMMTrendyNBP['meanNBP'+clusterkind+'Model']*factor#*10**9*44/12
    yStd = AllMMTrendyNBP['stdNBP'+clusterkind+'Model']*factor
    ax2.plot(AllMMTrendyNBP.MonthDate,y1,linewidth = lw,ls='-',marker= '',color= 'black',label='NBP')
    
    ax2.fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1>0),color= 'indigo',
                        alpha = 0.6,linewidth = lw,interpolate=True)
    ax2.fill_between(AllMMTrendyNBP.MonthDate,y1,0,where =(y1<0),color= 'green'
                        ,alpha = 0.4,linewidth = lw,interpolate=True)
                
    ax2.fill_between(AllMMTrendyNBP.MonthDate,
                        (y1-yStd),#*10**9*44/12,
                        (y1+yStd),#*10**9*44/12,
                        ls='-',linewidth = lw,color= 'grey',alpha = 0.4)
    #ax2.plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='--',marker= '',color= 'darkred',label=r'$\rm Inverse~Model_{+GOSAT}$')
    ax2.text(0.11, 0.95, r'$\rm \bf{B}$' +' Not semi-arid', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=9)
    #ax2.text(0.09, 0.95, r'$\rm \bf{A}$' +' Semi-arid', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=9)
    
    ax2.legend(fontsize=6,loc = 1,ncol = 3)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    ax2.set_ylim(-150,420)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    ax2.set_xlabel(r'Date',fontsize=7)
    #ax2.set_xticklabels(ax2.get_xticks(),fontsize=7)
    #ax2.set_yticklabels(ax2.get_yticks(),fontsize=7)
    ax2.set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2.tick_params(axis = 'both',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    ax2.set_axisbelow(True)
    
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp9_PanelA_goodModels_NewNewColor_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp9_PanelA_goodModels_NewNewColor_"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')

#S6 Annual Fluxes
if False:
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (14.4*cm,6.1*cm))
    secAx = ax2.twinx()
    secAx.plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',markersize=1,color= 'grey',zorder = 4)
    
    ACOSAf = allSatModel.groupby(['Year'])['TM5ACOSIS_nbp'].sum().reset_index()
    RTAf = allSatModel.groupby(['Year'])['TM5RTIS_nbp'].sum().reset_index()
    MeanAf = allSatModel.groupby(['Year'])['mean_nbp'].sum().reset_index()
    #Figure 8 MA  
    yeard = []
    RTAfd = []
    ACOSAfd = []
    MeanAfd = []
    MeanAfd2 = []
    yeard2 = [datetime.date(2009,4,15),datetime.date(2009,12,31)]
    for i in range(year_min,year_max+1):    
        yeard.append(datetime.date(i,6,30))
        if i >2009:
            yeard2.append(datetime.date(i,1,1))
            yeard2.append(datetime.date(i,12,31))
        ACOSAfd = ACOSAfd + [ACOSAf[(ACOSAf.Year == i)].TM5ACOSIS_nbp.values[0]]
        RTAfd = RTAfd + [RTAf[(RTAf.Year == i)].TM5RTIS_nbp.values[0]]
        MeanAfd = MeanAfd + [MeanAf[(MeanAf.Year == i)].mean_nbp.values[0],MeanAf[(MeanAf.Year == i)].mean_nbp.values[0]]
        MeanAfd2 = MeanAfd2 + [MeanAf[(MeanAf.Year == i)].mean_nbp.values[0]]
    secAx.fill_between(yeard2,np.array(MeanAfd)*factor,color= 'red', alpha=0.3, label = 'monthly fluxes',zorder = 3)
    secAx.errorbar(yeard,np.array(MeanAfd2)*factor,(np.array(RTAfd)-np.array(ACOSAfd))*factor/2,ls = '',color = 'grey',capsize = cs,linewidth = 0.8, zorder = 3, alpha=0.5)

    ax2.set_yticklabels([])
    ax2.tick_params(axis = 'y',length = 0.01, zorder = 0,labelsize = 6)
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    
    secAx.legend(fontsize=6,loc = 3,ncol = 1)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    secAx.set_ylim(-610,610)
    secAx.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    ax2.tick_params(axis = 'x', zorder = 0,labelsize = 6)
    #ax2.set_xticklabels([])
    #ax2[0].set_yticklabels(ax2[0].get_yticks(),fontsize=7)
    secAx.tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    #ax2.set_axisbelow(True)
    ax2.set_xlabel(r'Date',fontsize=7)
    #ax2.set_zorder(secAx.get_zorder()+1)
    #ax2.set_frame_on(False)
    
    
    plt.savefig(savepath + "/Results/Plots/AGU/AnnualfluxBar_GOSAT_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')

#S6 Paper Supp Annual Fluxes July - June
if False:
    withIS = True
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig2, ax2 = plt.subplots(figsize = (14.4*cm,6.1*cm))
    secAx = ax2.twinx()
    secAx.plot([TM5F_region[0].MonthDate.values[0],TM5F_region[0].MonthDate.values[-1]],[0,0],ls='-',linewidth = lw,marker= '',markersize=1,color= 'grey',zorder = 0)
    
    ACOSAf = allSatModel.groupby(['JJYear'])['TM5ACOSIS_nbp'].sum().reset_index()
    RTAf = allSatModel.groupby(['JJYear'])['TM5RTIS_nbp'].sum().reset_index()
    MeanAf = allSatModel.groupby(['JJYear'])['mean_nbp'].sum().reset_index()
    #Figure 8 MA  
    yeard = []
    RTAfd = []
    ACOSAfd = []
    MeanAfd = []
    MeanAfd2 = []
    yeard2 = [datetime.date(2009,7,1),datetime.date(2010,6,30)]
    for i in range(year_min,year_max):    
        yeard.append(datetime.date(i,12,31))
        if i >2009:
            yeard2.append(datetime.date(i,7,1))
            yeard2.append(datetime.date(i+1,6,30))
        ACOSAfd = ACOSAfd + [ACOSAf[(ACOSAf.JJYear == i)].TM5ACOSIS_nbp.values[0]]
        RTAfd = RTAfd + [RTAf[(RTAf.JJYear == i)].TM5RTIS_nbp.values[0]]
        MeanAfd = MeanAfd + [MeanAf[(MeanAf.JJYear == i)].mean_nbp.values[0],MeanAf[(MeanAf.JJYear == i)].mean_nbp.values[0]]
        MeanAfd2 = MeanAfd2 + [MeanAf[(MeanAf.JJYear == i)].mean_nbp.values[0]]
    
    if withIS:
        CAMSAf = allModel.groupby(['JJYear'])['CAMSsur_nbp'].sum().reset_index()
        CTAf = allModel.groupby(['JJYear'])['CT_nbp'].sum().reset_index()
        TM5Af = allModel.groupby(['JJYear'])['TM5_nbp'].sum().reset_index()
        MeanAfIS = allModel.groupby(['JJYear'])['mean_nbp'].sum().reset_index()
        #Figure 8 MA  
        yeard = []
        MaxAfdIS = []
        MinAfdIS = []
        MeanAfdIS = []
        MeanAfd2IS = []
        yeard2 = [datetime.date(2009,7,1),datetime.date(2010,6,30)]
        for i in range(year_min,year_max):    
            yeard.append(datetime.date(i,12,31))
            if i >2009:
                yeard2.append(datetime.date(i,7,1))
                yeard2.append(datetime.date(i+1,6,30))
            MinAfdIS = MinAfdIS + [min([CAMSAf[(CAMSAf.JJYear == i)].CAMSsur_nbp.values[0],CTAf[(CTAf.JJYear == i)].CT_nbp.values[0],TM5Af[(TM5Af.JJYear == i)].TM5_nbp.values[0]])]
            MaxAfdIS = MaxAfdIS + [max([CAMSAf[(CAMSAf.JJYear == i)].CAMSsur_nbp.values[0],CTAf[(CTAf.JJYear == i)].CT_nbp.values[0],TM5Af[(TM5Af.JJYear == i)].TM5_nbp.values[0]])]
            MeanAfdIS = MeanAfdIS + [MeanAfIS[(MeanAfIS.JJYear == i)].mean_nbp.values[0],MeanAfIS[(MeanAfIS.JJYear == i)].mean_nbp.values[0]]
            MeanAfd2IS = MeanAfd2IS + [MeanAfIS[(MeanAfIS.JJYear == i)].mean_nbp.values[0]]
        secAx.fill_between(yeard2,np.array(MeanAfdIS)*factor,color= 'blue', alpha=0.4,label = r"Inverse mode$\rm l_{in-situ}$",zorder = 6)
        secAx.errorbar(yeard,np.array(MeanAfd2IS)*factor,(np.array(MaxAfdIS)-np.array(MinAfdIS))*factor/2,ls = '',color = 'black',capsize = cs,linewidth = 0.8, zorder = 3, alpha=0.4)#, label = r'$\rm range~GOSAT_{RemoTeC}~GOSAT_{ACOS}$')
        
        ax2.fill_between(allModel.MonthDate,allModel.max_nbp*factor,allModel.min_nbp*factor,color= 'blue', alpha=0.1,zorder = 4)  
        ax2.plot(allModel.MonthDate,allModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'blue',label=r"Inverse mode$\rm l_{in-situ}$",zorder = 6)
        
    secAx.fill_between(yeard2,np.array(MeanAfd)*factor,color= 'red', alpha=0.3, label = r"Inverse mode$\rm l_{+GOSAT}$",zorder = 6)secAx.errorbar(yeard,np.array(MeanAfd2)*factor,(np.array(RTAfd)-np.array(ACOSAfd))*factor/2,ls = '',color = 'black',capsize = cs,linewidth = 0.8, zorder = 3, alpha=0.4)#, label = r'$\rm range~GOSAT_{RemoTeC}~GOSAT_{ACOS}$')
    
    ax2.fill_between(allSatModel.MonthDate,allSatModel.max_nbp*factor,allSatModel.min_nbp*factor,color= 'red', alpha=0.1,zorder = 4)  
    ax2.plot(allSatModel.MonthDate,allSatModel.mean_nbp*factor,ls='-',linewidth = lw,marker= '',color= 'darkred',label=r"Inverse mode$\rm l_{+GOSAT}$",zorder = 6)#label="monthly fluxes\n"+r"(TM5-4DVa$\rm r_{GOSAT+in-situ})$",zorder = 6)
    
    ax2.text(0.075, 0.21, r'$\rm Monthly~flux$', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=6, weight='bold')
    ax2.text(0.86, 0.21, r'$\rm Annual~flux~June-July$', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=6, weight='bold')
    

    ax2.set_ylim(-280,280)
    ax2.set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    ax2.tick_params(axis = 'y', zorder = 0,labelsize = 6)
    ax2.legend(fontsize=6,loc = 3,ncol = 1)#,loc=3)#, loc = 9)
    
    ax2.grid(True,which = 'both', axis='x',zorder = 0)
    
    secAx.legend(fontsize=6,loc = 4,ncol = 1)#,loc=3)#, loc = 9)
    #ax2.set_zorder(ax2.get_zorder()+1)
    #ax2.patch.set_visible(False)
    secAx.set_ylim(-760,760)
    secAx.set_ylabel(r'$\rm CO_2~flux~(TgC/year)$',fontsize=7)
    ax2.tick_params(axis = 'x', zorder = 0,labelsize = 6)
    #ax2.set_xticklabels([])
    #ax2[0].set_yticklabels(ax2[0].get_yticks(),fontsize=7)
    secAx.tick_params(axis = 'y',labelsize = 6, zorder = 0)
    ax2.set_xlim(datetime.date(2008,12,15),datetime.date(2019,3,15))
    #ax2.set_axisbelow(True)
    ax2.set_xlabel(r'Date',fontsize=7)
    #ax2.set_zorder(secAx.get_zorder()+1)
    #ax2.set_frame_on(False)
    
    
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/AnnualfluxBarJJ_MonthFlux_GOSAT_IS_newShading2_Region43TM5_"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')#,format = 'pdf')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/AnnualfluxBarJJ_MonthFlux_GOSAT_IS_newShading2_Region43TM5_"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')


#S10 Mean Cycle for individual models
if True:
    if True:#mean plot July to June, JSBACH, OCN, YIBs, CLASSIC, LPJ, POSSIBLE without 2010 and 2011 
        factor = 12/44
        lw = 0.8
        cm = 1/2.54  # centimeters in inches
    
        fig2, ax2 = plt.subplots(2, 1, 
                            figsize = (12.0*cm,16*cm),
                            gridspec_kw={'height_ratios':[1,1]})
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
        indexrange = [6,7,8,9,10,11,0,1,2,3,4,5]

        for k in range(2):
            TRENDYwoLaNina = AllMMTrendyNBP.copy(deep = True)
        
            LModelkind = ['good','bad']
            Modelkind = LModelkind[k]

            def Modellist(Modkind, vari):
                if Modkind == 'good':
                    listM = GoodModelsNames(vari)
                elif Modkind == 'early':
                    listM = EarlyGPPModelNames(vari)
                elif Modkind == 'late':
                    listM = LateRespNames(vari)
                elif Modkind == 'bad':
                    listM = LateRespNames(vari) + EarlyGPPModelNames(vari)
                return listM
                

            MeanCycgoodgpp = TRENDYwoLaNina.groupby('Month')['meangpp'+Modelkind+'Model'].mean().reset_index()
            MeanCycgoodra = TRENDYwoLaNina.groupby('Month')['meanra'+Modelkind+'Model'].mean().reset_index()
            MeanCycgoodrh = TRENDYwoLaNina.groupby('Month')['meanrh'+Modelkind+'Model'].mean().reset_index()
            MeanCycgoodrarh = TRENDYwoLaNina.groupby('Month')['meanrhra'+Modelkind+'Model'].mean().reset_index()
            #STd Dev over Time
            StdCycgoodgpp = TRENDYwoLaNina.groupby('Month')['meangpp'+Modelkind+'Model'].std(ddof=0).reset_index()
            StdCycgoodra = TRENDYwoLaNina.groupby('Month')['meanra'+Modelkind+'Model'].std(ddof=0).reset_index()
            StdCycgoodrh = TRENDYwoLaNina.groupby('Month')['meanrh'+Modelkind+'Model'].std(ddof=0).reset_index()
            StdCycgoodrarh = TRENDYwoLaNina.groupby('Month')['meanrhra'+Modelkind+'Model'].std(ddof=0).reset_index()
            #StdOverModels
            StdMCycgoodgpp = TRENDYwoLaNina.groupby('Month')[Modellist(Modelkind, 'gpp')].mean().reset_index()[Modellist(Modelkind, 'gpp')].std(axis=1,ddof=0).reset_index()
            StdMCycgoodra = TRENDYwoLaNina.groupby('Month')[Modellist(Modelkind, 'ra')].mean().reset_index()[Modellist(Modelkind, 'ra')].std(axis=1,ddof=0).reset_index()
            StdMCycgoodrh = TRENDYwoLaNina.groupby('Month')[Modellist(Modelkind, 'rh')].mean().reset_index()[Modellist(Modelkind, 'rh')].std(axis=1,ddof=0).reset_index()
            StdMCycgoodgpp.rename(columns={0:'meangpp'+Modelkind+'Model'},inplace=True)
            StdMCycgoodra.rename(columns={0:'meanra'+Modelkind+'Model'},inplace=True)
            StdMCycgoodrh.rename(columns={0:'meanrh'+Modelkind+'Model'},inplace=True)

            #Std over models
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange]-StdMCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange]+StdMCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'green', alpha=0.3,zorder = 2)
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange]-StdMCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange]+StdMCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'blue', alpha=0.3,zorder = 2)
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange]-StdMCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange]+StdMCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'brown', alpha=0.3,zorder = 2)

            #Std over time
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange]-StdCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange]+StdCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'green',lw=lw, alpha=0.3,zorder = 2,hatch='\\\\',)
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange]-StdCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange]+StdCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'blue',lw=lw, alpha=0.3,zorder = 2,hatch='\\\\',)
            ax2[k].fill_between(range(1,13),
                            (MeanCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange]-StdCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange])*factor,
                            (MeanCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange]+StdCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange])*factor,
                            color= 'brown',lw=lw, alpha=0.3,zorder = 2,hatch='\\\\',)

                    
            ax2[k].plot(range(1,13),MeanCycgoodgpp['meangpp'+Modelkind+'Model'].iloc[indexrange]*factor,ls='-',lw=lw,marker= '',color= 'green',label='GPP')
            ax2[k].plot(range(1,13),MeanCycgoodra['meanra'+Modelkind+'Model'].iloc[indexrange]*factor,ls='-',lw=lw,marker= '',color= 'blue',label='Plant respiration')
            ax2[k].plot(range(1,13),MeanCycgoodrh['meanrh'+Modelkind+'Model'].iloc[indexrange]*factor,ls='-',lw=lw,marker= '',color= 'brown',label='Soil respiration')
            ax2[k].plot(range(1,13),MeanCycgoodrarh['meanrhra'+Modelkind+'Model'].iloc[indexrange]*factor,ls='--',lw=lw,marker= '',color= 'black',label='Total respiration')
            
            #ax2[k].text(0.03, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[k].transAxes,fontsize=10, weight='bold')

            ax2[k].legend(fontsize=6,loc=1)
            ax2[k].grid(True,which = 'both', axis='x')
            ax2[k].set_xlabel(r'Month',fontsize=6)
            ax2[k].set_xticks(ticks = [3,6,9,12])
            ax2[k].set_xticklabels(['9','12','3','6'])
            ax2[k].tick_params(axis = 'both',labelsize = 6)
            ax2[k].set_ylabel(r'$\rm CO_2~flux~(TgC/month)$',fontsize=7)
    
        ax2[0].text(0.05, 0.93, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10, fontweight='bold')
        ax2[0].text(0.08, 0.93, 'Selection', horizontalalignment='left', verticalalignment='center', transform=ax2[0].transAxes,fontsize=10)
        ax2[1].text(0.05, 0.93, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10, fontweight='bold')
        ax2[1].text(0.08, 0.93, 'Other', horizontalalignment='left', verticalalignment='center', transform=ax2[1].transAxes,fontsize=10)
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
        
    
        plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp10_models_rhra_rh_ra_gpp"+RegionName+str(Numm)+".png", dpi=400,bbox_inches='tight')
        plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp10_models_rhra_rh_ra_gpp"+RegionName+str(Numm)+".pdf", dpi=400,bbox_inches='tight',format = 'pdf')
    
