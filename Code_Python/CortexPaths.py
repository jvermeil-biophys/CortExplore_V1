# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:51:13 2022

@author: JosephVermeil
"""

# %% 0. Imports

import os
import sys
import ctypes
from datetime import date

# %% 1. Paths

COMPUTERNAME = os.environ['COMPUTERNAME']

# 1.1 Init main directories

if COMPUTERNAME == 'ORDI-JOSEPH':
    DirRepo = "C://Users//JosephVermeil//Desktop//CortExplore"
    DirData = "D://MagneticPincherData"
    DirCloud = "C://Users//JosephVermeil//ownCloud//MagneticPincherData_JV"
    DirTempPlots = "C://Users//JosephVermeil//Desktop//TempPlots"
    CloudSaving = 'OwnCloud'
    
    
elif COMPUTERNAME == 'LARISA':
    DirRepo = "C://Users//Joseph//Desktop//CortExplore"
    DirData = "F://JosephVermeil//MagneticPincherData"
    DirCloud = "C://Users//Joseph//ownCloud//MagneticPincherData_JV"
    DirTempPlots = 'C://Users//Joseph//Desktop//TempPlots'
    CloudSaving = 'OwnCloud'
    
    
elif COMPUTERNAME == 'DESKTOP-K9KOJR2':
    DirRepo = "C://Users//anumi//OneDrive//Desktop//CortExplore"
    DirData = "D:/Anumita/MagneticPincherData"
    DirCloud = ""
    DirTempPlots = 'C://Users//anumi//OneDrive//Desktop//TempPlots'
    CloudSaving = ''
    
    
elif COMPUTERNAME =='DATA2JHODR':
    DirRepo = "C://Users//BioMecaCell//Desktop//CortExplore"
    DirData = "D:/Anumita/MagneticPincherData"
    DirCloud = ""
    DirTempPlots = 'C://Users//BioMecaCell//Desktop//TempPlots'
    CloudSaving = ''
    
    
# 1.2 Init sub directories

DirRepoPython = os.path.join(DirRepo, "Code_Python")
DirRepoExp = os.path.join(DirRepo, ".Data_Experimental")

DirDataRaw = os.path.join(DirData, "Raw")
DirDataRawDeptho = os.path.join(DirDataRaw, 'DepthoLibrary')
DirDataRawDepthoInter = os.path.join(DirDataRawDeptho, 'IntermediateSteps')

DirDataExp = os.path.join(DirData, "Data_Experimental")

DirDataAnalysis = os.path.join(DirData, "Data_Analysis")
DirDataTimeseries = os.path.join(DirDataAnalysis, "Data_Timeseries")
DirDataTimeseriesRawtraj = os.path.join(DirDataTimeseries, "Trajectories_raw")
DirDataTimeseriesTraj = os.path.join(DirDataTimeseries, "Trajectories")

DirDataFig = os.path.join(DirData, "Figures")
DirDataFigToday = os.path.join(DirDataFig, "Historique", str(date.today()))

if not CloudSaving == '':
    DirCloudExp = os.path.join(DirCloud, "Data_Experimental")
    DirCloudAnalysis = os.path.join(DirCloud, "Data_Analysis")
    DirCloudTimeseries = os.path.join(DirCloud, "Data_Timeseries")
    DirCloudFig = os.path.join(DirCloud, "Figures")
    DirCloudFigToday = os.path.join(DirCloudFig, "Historique", str(date.today()))
else:
    DirCloudExp, DirCloudFig, DirCloudFigToday = "", "", ""

# 1.3 Add python directory to path

sys.path.append(DirRepoPython)


# %% 2. Useful functions

MainDirs = [DirRepo, DirData, DirTempPlots]
RepoSubdirs = [DirRepoPython, DirRepoExp]
DataSubdirs = [DirDataRaw, DirDataExp, DirDataAnalysis, DirDataFig,
               DirDataTimeseries, DirDataTimeseriesTraj, DirDataTimeseriesRawtraj]
CloudDirs = [DirCloud, DirCloudExp, DirCloudFig, DirCloudAnalysis, DirCloudTimeseries]


def checkDirArchi():
    valid_main = True
    for p in MainDirs:
        if not os.path.exists(p):
            valid_main = False
    
    if not valid_main:
        print('One of the main directories is missing')
        
    else:
        valid_repo, valid_data, valid_cloud = True, True, True
        for p in RepoSubdirs:
            if not os.path.exists(p):
                valid_repo = False
        if not valid_repo:
            print('One of the repository sub-directories is missing')
                
        for p in DataSubdirs:
            if not os.path.exists(p):
                valid_repo = False
        if not valid_repo:
            print('One of the data sub-directories is missing')
        
        if not CloudSaving == '':
            for p in CloudDirs:
                if not os.path.exists(p):
                    valid_cloud = False
        if not valid_cloud:
            print('One of the cloud sub-directories is missing')
            
    if valid_main and valid_repo and valid_data and valid_cloud:
        print('Directories architecture is correct !')
            
            
def makeDirArchi():
    for p in MainDirs:
        if not os.path.exists(p):
            os.makedirs(p)
    
    for p in RepoSubdirs:
        if not os.path.exists(p):
            os.makedirs(p)
    
    if os.path.exists(DirRepoExp):
        FILE_ATTRIBUTE_HIDDEN = 0x02
        ret = ctypes.windll.kernel32.SetFileAttributesW(DirRepoExp, FILE_ATTRIBUTE_HIDDEN)
    
    warningRepoExp = os.path.join(DirRepoExp, 'Warning.txt')
    if not os.path.exists(warningRepoExp):
        f = open(warningRepoExp, "w")
        text = 'WARNING\nDo not modify this file. It is for consultation only.\n'
        text += 'For the modifiable version go to: ' + DirDataExp
        f.write(text)
        f.close()
            
    for p in DataSubdirs:
        if not os.path.exists(p):
            os.makedirs(p)

    if not CloudSaving == '':
        for p in CloudDirs:
            if not os.path.exists(p):
                os.makedirs(p)
        
        warningCloudExp = os.path.join(DirCloudExp, 'Warning.txt')
        if not os.path.exists(warningCloudExp):
            f = open(warningCloudExp, "w")
            text = 'WARNING\nDo not modify this file. It is for consultation only.\n'
            text += 'For the modifiable version go to: ' + DirDataExp
            f.write(text)
            f.close()

                
                
# %% Final Architecture (Ongoing)

# C:/
# ├─ Users/
# │  ├─ User/
# │  │  ├─ Desktop/
# │  │  │  ├─ CortExplore/
# │  │  │  │  ├─ .Data_BeadsCalibration/
# │  │  │  │  ├─ .Data_Experimental/
# │  │  │  │  ├─ .git/
# │  │  │  │  ├─ Code_IJ/
# │  │  │  │  ├─ Code_Matlab/
# │  │  │  │  ├─ Code_Python/
# │  │  │  │  ├─ LICENSE
# │  │  │  │  ├─ README.md
# │  │  ├─ ownCloud/
# │  │  │  ├─ MagneticPincherData_##/
# │  │  │  │  ├─ Data_Analysis/
# │  │  │  │  │  ├─ new_folder/
# │  │  │  │  │  ├─ new_folder/
# │  │  │  │  ├─ Data_Experiemental/
# │  │  │  │  ├─ Data_BeadsCalibration/
# │  │  │  │  ├─ Figures/
# D:/
# ├─ MagneticPincherData/
# │  ├─ Data_Analysis/
# │  │  ├─ GlobalTable.csv
# │  ├─ Data_BeadsCalibration/
# │  ├─ Data_Experimental/
# │  ├─ Data_TimeSeries/
# │  │  ├─ Trajectories/
# │  │  ├─ Trajectories_raw/
# │  │  ├─ yy-mm-dd_M#_P#_C#.csv
# │  ├─ Figures/
# │  ├─ Raw/
# │  │  ├─ DepthoLibrary/
# │  │  ├─ yy.mm.dd/
# │  │  ├─ yy.mm.dd_Deptho/
# new_folder/
# new_folder/
# new_folder/
# new_folder/
# new_file

                
                
                
                
                
                
                
                
                
                
            