# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:20:48 2021

@author: Anumita Jawahar
"""

# %% General imports

# 1. Imports
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

import os
import re
import time
import pyautogui
import matplotlib

from scipy import interpolate
from scipy import signal
from skimage import io, filters, exposure, measure, transform, util
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import linear_sum_assignment

# Local Imports

import sys
import CortexPaths as cp
sys.path.append(cp.DirRepoPython)

import GraphicStyles as gs
import GlobalConstants as gc
import UtilityFunctions as ufun

from BeadTracker import depthoMaker

# 2. Pandas settings
pd.set_option('mode.chained_assignment',None)

# 3. Graphical settings
# gs.set_default_options_jv()

# 6. Others
SCALE_100X = 15.8 # pix/µm
SCALE_63X = 9.9 # pix/µm


# %% EXAMPLE -- All depthos from 21.01.18 3T3 experiments

mainDirPath = 'D://MagneticPincherData//Raw'


date = '21.01.18'
subdir = 'Deptho_M1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)



subdir = 'Deptho_M2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


subdir = 'Deptho_M3'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M3_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

# %% Test deptho 60x Oil Atchoum 21/11/30
mainDirPath = 'D:/Anumita/Data/Raw/'

date = '21.11.30'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step50_60X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = gc.SCALE_60X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 50, d = 'HD', plot = 0)

# %% Test deptho 100x oil Atchoum 21/12/08

mainDirPath = 'D:/Anumita/Data/Raw/'

date = '21.12.08'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step50_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 50, d = 'HD', plot = 0)

#%% Test 2 deptho 100x oil Atchoum 21/12/13 - After optimizing illumination with Joseph


mainDirPath = 'D:/Anumita/Data/Raw/'

date = '21.12.13'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 21-12-20, not so good because beads are floating. Taken at end of the experiment.


mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '21.12.20'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


#%% Deptho from experiment 22-02-03, taken at  the beginning of the experiment


mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.02.03'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-03-01, taken at  the beginning of the experiment


mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.03.01'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-03-22, taken in the middle of the experiment on 
#beads that looked stuck on fibronectin patterns.

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.03.22'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

subdir = 'Deptho_P1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

subdir = 'Deptho_P2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-03-31, taken in the middle of the experiment on 
#beads that looked stuck on fibronectin patterns at PMMH. First mechanics experiment.

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.03.31'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

subdir = 'Deptho_P1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

# subdir = 'Deptho_P2'
# depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
# savePath = os.path.join(mainDirPath, 'DepthoLibrary')

# specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
# beadType = 'M450'
# saveLabel = date + '_P2_M450_step20_100X'
# # convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
# scale = SCALE_100X # pix/µm

# depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-04-05, taken in the beginning of the experiment on 
#beads that looked stuck on fibronectin patterns at PMMH. 

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.04.05'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

subdir = 'Deptho_P1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

# subdir = 'Deptho_P2'
# depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
# savePath = os.path.join(mainDirPath, 'DepthoLibrary')

# specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
# beadType = 'M450'
# saveLabel = date + '_P2_M450_step20_100X'
# # convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
# scale = SCALE_100X # pix/µm

# depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-04-12, First constant field expeirment in PMMH. 
# img/ml PEG+HEPES beads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.04.12'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


#%% Deptho from experiment 22-04-28. Mechanics #3
#  1 mg/ml PEG+HEPES beads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.04.28'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)


#%% Deptho from experiment 22-04-28. Mechanics #3
#  1 mg/ml PEG+HEPES beads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.05.09'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

subdir = 'Deptho_P1'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P1_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

subdir = 'Deptho_P2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-05-31. Mechanics #3
#  1 mg/ml PEG+HEPES beads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.05.31'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

# subdir = 'Deptho_P1'
# depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
# savePath = os.path.join(mainDirPath, 'DepthoLibrary')

# specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
# beadType = 'M450'
# saveLabel = date + '_P1_M450_step20_100X'
# # convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
# scale = SCALE_100X # pix/µm

# depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

subdir = 'Deptho_P2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-07-26. 
#  0.1 mg/ml mPEG-Biotin + Streptavidin Dynabeads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.07.26'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

#%% Deptho from experiment 22-08-26.
#  0.1 mg/ml mPEG-Biotin + Streptavidin Dynabeads

mainDirPath = 'D:/Anumita/MagneticPincherData/Raw/'

date = '22.08.26'
depthoPath = os.path.join(mainDirPath, date + '_Deptho')
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

subdir = 'Deptho_P2'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P2_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)

subdir = 'Deptho_P3'
depthoPath = os.path.join(mainDirPath, date + '_Deptho', subdir)
savePath = os.path.join(mainDirPath, 'DepthoLibrary')

specif = 'all' # can be 'all' or any string that you want to have in the deptho file name
beadType = 'M450'
saveLabel = date + '_P3_M450_step20_100X'
# convention - saveLabel = 'date_manip_beadType_stepSize_otherSpecs'
scale = SCALE_100X # pix/µm

depthoMaker(depthoPath, savePath, specif, saveLabel, scale, beadType = beadType, step = 20, d = 'HD', plot = 0)