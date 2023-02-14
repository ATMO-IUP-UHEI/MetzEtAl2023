#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:46:39 2021

@author: eschoema
Function to get Region Name and Coordinates
"""

def getRegion(Numm):
    #Regions not based on TRANSCOM Regions
    region_names = ['North America South', #0 Guerlet
                'North America North', #1 Guerlet
                'Eurasia South',       #2 Guerlet
                'Eurasia North',       #3 Guerlet
                'Asia South',          #4 Guerlet
                'Asia North',          #5 Guerlet
                'South America',       #6 Guerlet
                'South Africa',        #7 Guerlet
                'Australia',           #8 Guerlet
                'Australia',           #9 Detmers
                'Australia',           #10
                'South Africa',        #11
                'South America',       #12
                'Whole Australia',     #13
                'South Africa',        #14 Guerlet Africa moved north
                'South Africa',        #15 Guerlet Africa moved south
                'South America',       #16 Guerlet S Amerika moved north
                'South America',       #17 Guerlet S Amerika moved south
                'Australia',           #18 Guerlet Australia move south
                'Australia',           #19 Guerlet Australia moved completely to south
                'South Africa',        #20 whole South Afrika ACHTUNG Rechteckig!
                'South America',       #21 whole South America ACHTUNG Rechteckig!
                'Australia',           #22 whole Australia ACHTUNG Rechteckig!
                'Australia',           #23 Region around Darwin
                'Australia',           #24 Region around Wollongong
                'Australia',           #25 Region around Cape Grim 5°
                'Australia',           #26 Region around Cape Grim 10°
                'Whole Australia',     #27 whole australia + 15° buffer
                'AIsland',             #28 Region around A Island 15°
                'Reunion',             #29 Region around Reunion 15°
                'Lauder',              #30 Region around Lauder 15°
                'Manaus',              #31 Region around Manaus
                'South America',       #32 Whole South America
                'South Africa',        #33 Whole South Africa
                'world',               #34 no filter
                'South Africa',        #35 Lat: 0-15
                'South Africa',        #36 Lat: 5-20
                'South Africa',        #37 Lat: 10-25
                'South Africa',        #38 Lat: 15-30
                'South Africa',        #39 Lat: 20-35
                'South America',       #40 Lat: 0-15
                'South America',       #41 Lat: 5-20
                'South America',       #42 Lat: 10-25
                'South America',       #43 Lat: 15-30
                'South America',       #44 Lat: 20-35
                'South America',       #45 Lat: 25-40
                'South America',       #46 Lat: 30-45
                'South America',       #47 Lat: 35-50
                'South America',       #48 Lat: 40-55
                'Australia',           #49 east australia
                'Australia',           #50 west australia
                'NH',                  #51 whole northern hemisphere
                'SH',                  #52 whole southern hemisphere
                'Parazoo1',            #53 Parazoo et al 2013 Region 1
                'Parazoo2',            #54 Parazoo et al 2013 Region 2
                'Australia',           #55 Box around Australia corresponding to 949 without transcom
                ]


    #Region =   [0   ,1    ,2  ,3  ,4    ,5    ,6   ,7   ,8   ,9   ,10  ,11  ,12  ,13  ,14  ,15]
    Long_minl = [-120,-124 ,25 ,25 ,82   ,82   ,-72 ,12.5,114 ,130 ,120 ,15 ,-72 ,133 ,12.5,12.5]
    Long_maxl = [-70 ,-47.5,80 ,80 ,137.5,137.5,-40 ,50  ,154 ,140 ,150 ,32 ,-50 ,154 ,50  ,50]
    Lat_minl =  [40  ,52   ,40 ,52 ,40   ,52   ,-30 ,-30 ,-30 ,-30 ,-30 ,-27,-30 ,-39 ,-25 ,-35]
    Lat_maxl =  [52  ,64   ,52 ,64 ,52   ,64   ,-15 ,-10 ,-12 ,-20 ,-20 ,-10,-15 ,-11 ,-5  ,-15]
    
    #Region =    16  ,17 ,18 , 19 ,20 ,21 ,22 , 23,   24,   25, 26  ,27  ,28  , 29 , 30
    Long_minl +=[ -72 ,-72,114,114 ,10 ,-80,114,126  ,146  ,140,135 , 100, -30, 40 , 152]
    Long_maxl +=[ -40 ,-40,154,154 ,50 ,-35,154,136  ,156  ,150,155 , 170, 1  , 71 , 180]
    Lat_minl  +=[ -35 ,-25,-35,-40 ,-35,-55,-40,-17.5,-39.5,-46,-51 , -60, -23, -36, -60]
    Lat_maxl  +=[ -20 ,-10 ,-17,-22 ,0  ,0  ,-12,-7.5 ,-29.5,-36,-31 ,  5 , 8  , -5 , -30]
    
    #Region =    31  ,32  ,33  ,34  , 35,36 ,37 ,38 ,39 ,40 ,41 ,42 ,43 ,44 ,45 ]
    Long_minl +=[-76 ,-84 , 8  ,-180, 8 , 8 , 8 , 8 , 8 ,-83,-83,-83,-83,-83,-83]
    Long_maxl +=[-45 ,-33 , 52 ,180 , 53, 53, 53, 53, 53,-33,-33,-33,-33,-33,-33]
    Lat_minl  +=[-19 ,-56 , -36,-90 ,-15,-20,-25,-30,-35,-15,-20,-25,-30,-35,-40 ]
    Lat_maxl  +=[12  ,  0 , 0  ,90  , 0 ,-5 ,-10,-15,-20, 0 ,-5 ,-10,-15,-20,-25 ]

    #Region =    46 ,47 ,48 ,49 ,50 ,51  ,52  ,53 ,54 ,55]
    Long_minl +=[-83,-83,-83,144,114,-180,-180,-62,-50,122]  
    Long_maxl +=[-33,-33,-33,154,124,180 ,180 ,-51,-40,180]
    Lat_minl  +=[-45,-50,-55,-40,-40, -90,0   ,-15,-16,-50]
    Lat_maxl  +=[-30,-35,-40,-12,-12,0   ,90  ,-8 ,-5 ,-10]



    #Regions based on TRANSCOM Regions, own mask
    Transcom = ['SAT','SA','AU','SATr','SATr_SE','I_P','SAT','SAT','SAT','SAT', #0-9
                'SAT','SAT','SA','SA','SA','SA','SA','AU','AU','AU', #10-19
                'AU','AU','AU','AU','AU','AU','AU','AU','AU','AU', #20-29
                'SAT','SAT','SAT','SAT','SA','SA','SA','Af','SAT_SATr','TA',#30-39
                'EU','NAS','As','NAN','EB','AU','AU','AU','AU','AU', #40-49
                'AU','AU','AU','AU' #50-59
                ]


    coordt = [[-56,0,-84,-32],
          [-36,0,7,57],
          [-48,-12,114,180], 
          [-12,15,-84,-32], 
          [-12,0,-84,-32], #4
          [-12,6,90,160],
          [-56,0,-80,-60], #Longitudinal investigations
          [-56,0,-75,-55],
          [-56,0,-70,-50],
          [-56,0,-65,-45],
          [-56,0,-60,-40], #10
          [-56,0,-55,-35],
          [-36,0,7,30],
          [-36,0,15,35],
          [-36,0,20,40],
          [-36,0,25,45], #15
          [-36,0,30,57],
          [-48,-12,114,135],
          [-48,-12,120,140],
          [-48,-12,125,145],
          [-48,-12,130,150], #20
          [-48,-12,135,160],
          [-48,-12,160,180],
          [-48,-12,114,125],
          [-48,-12,125,135],
          [-48,-12,135,145], #25
          [-48,-12,145,155], #------------------------------------
          [-24,-12,114,180], #Latitudinal investigations
          [-36,-24,114,180],
          [-48,-36,114,180],
          [-15,0,-84,-32],   #30
          [-25,-15,-84,-32],
          [-40,-25,-84,-32],
          [-60,-40,-84,-32],
          [-10,0,7,57],
          [-20,-10,7,57],    #35
          [-36,-20,7,57],
          [-36,40,-20,50],      #whole Africa
          [-56,15,-84,-32],     #whole southern America
          [-12,30,75,160],      #Tropical asia
          [30,85,-30,180],  #40 #Europe + Russia
          [0,49,-180,-10],      #Northern America South
          [0,60,30,150],         #Asia withoupt TrAsia and Russia
          [49,90,-180,-10],      #Northern America North
          [40,90,60,180] ,      #Eurasia Boreal
          [-30,-22,125,145],  #45    # Australien Gain M H subregion
          [-22,-10,114,180], #Latitudinal investigations 2
          [-29,-22,114,180],
          [-50,-29,114,180],  
          [-50,-10,112,180], #49 new main region from 03.11.2021
          [-50,-10,112,180], #50 used in Kombination with Semi-arid Mask
          [-50,-10,112,180], #51 used in Kombination with Wet Mask
          [-50,-10,112,180], #52 used in Kombination with Semi-arid Mask based on Land Cover
          [-50,-10,112,180], #53 used in Kombination with Wet Mask Mask based on Land Cover
          ]

    #Regions based on TRANSCOM Regions, CarbonTracker mask
    TranscomCT = ['AU']
    coordtCT = [[-50,-10,112,180],#same as for 949 but with CT mask
                ]

    if Numm <800:
        Long_min = Long_minl[Numm]
        Long_max = Long_maxl[Numm]
        Lat_min = Lat_minl[Numm]
        Lat_max = Lat_maxl[Numm]
    
        RegionName = region_names[Numm]
    elif Numm < 900:
        Long_min = coordtCT[Numm - 800][2]
        Long_max = coordtCT[Numm - 800][3]
        Lat_min = coordtCT[Numm - 800][0]
        Lat_max = coordtCT[Numm - 800][1]
    
        RegionName = TranscomCT[Numm - 800]
    else:
        Long_min = coordt[Numm - 900][2]
        Long_max = coordt[Numm - 900][3]
        Lat_min = coordt[Numm - 900][0]
        Lat_max = coordt[Numm - 900][1]
    
        RegionName = Transcom[Numm - 900]
    
    return RegionName, Long_min, Long_max, Lat_min, Lat_max

