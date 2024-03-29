# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:50:53 2021

@author: Joseph Vermeil
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
sys.path.append(cp.DirRepoPythonUser)


import GraphicStyles as gs
import GlobalConstants as gc
import UtilityFunctions as ufun

from BeadTracker import mainTracker


# 2. Pandas settings
pd.set_option('mode.chained_assignment',None)

# 3. Graphical settings
gs.set_default_options_jv()


# 4. Import of the experimental conditions

expDf = ufun.getExperimentalConditions(DirExp = cp.DirRepoExp, save = True, sep = ';', suffix = cp.suffix)

# %% Small things

#### Plot last traj

# fig, ax = plt.subplots(1,1)
# X1, Y1 = listTrajDicts[0]['X'], listTrajDicts[0]['Y']
# X2, Y2 = listTrajDicts[1]['X'], listTrajDicts[1]['Y']
# ax.plot(X1, Y1, 'r-')
# ax.plot(X2, Y2, 'b-')

# fig.show()

#### close all
plt.close('all')

# %% Next Topic !!

# %%% Next experiment day
# %%%% Next manipe
# %%%% Next manipe



# %%% 22.07.15, compressionsLowStart of HoxB8 macrophages, M450, M1 = ctrl & 20um discs, M2 = tko & 20um discs, M3 = tko & glass, M4 = ctrl & glass
# %%%% 22.07.15 one specific cell
dates = '22.07.15'
manips, wells, cells = 1, 1, 1
depthoNames = '22.07.15_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.07.15_M1
dates = '22.07.15'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.07.15_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.07.15_M2
dates = '22.07.15'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.07.15_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.07.15_M3
dates = '22.07.15'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.07.15_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.07.15_M4
dates = '22.07.15'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.07.15_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')



# %% HoxB8 macrophages

# %%% 22.05.05, compressionsLowStart of HoxB8 macrophages, M450, M1 = tko & glass, M2 = ctrl & glass
# %%%% 22.05.05_M1 C1 Seulement
dates = '22.05.05'
manips, wells, cells = 2, 1, 2
depthoNames = '22.05.05_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.05_M1
dates = '22.05.05'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.05.05_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.05_M2
dates = '22.05.05'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.05.05_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 22.05.04, compressionsLowStart of HoxB8 macrophages, M450, M1 = ctrl & 20um discs, M2 = tko & 20um discs, M3 = tko & glass, M4 = ctrl & glass
# %%%% 22.05.04 one specific cell
dates = '22.05.04'
manips, wells, cells = 2, 1, 8
depthoNames = '22.05.04_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.04_M1
dates = '22.05.04'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.05.04_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.04_M2
dates = '22.05.04'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.05.04_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.04_M3
dates = '22.05.04'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.05.04_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.04_M4
dates = '22.05.04'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.05.04_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')



# %%% 22.05.03, compressionsLowStart of HoxB8 macrophages, M450, M1 = ctrl & glass, M2 = tko & glass, M3 = tko & 20um discs, M4 = ctrl & 20um discs
# %%%% 22.05.03_M1 C1 Seulement
dates = '22.05.03'
manips, wells, cells = ['1-1'], 1, 1
depthoNames = '22.05.03_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.03_M1
dates = '22.05.03'
manips, wells, cells = ['1-1', '1-2'], 1, 'all'
depthoNames = '22.05.03_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.03_M2
dates = '22.05.03'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.05.03_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.03_M3
dates = '22.05.03'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.05.03_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.05.03_M4
dates = '22.05.03'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.05.03_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %% Drugs & perturbation

# %%% Next experiment day
# %%%% Next manipe

# %%% 22.03.30, compressionsLowStart of 3T3 LG +++, M450, M1 = Blebbi, M2 = LatA, M3 = Ctrl, M4 = DMSO
# %%%% 22.03.30_M1 C1 Seulement
dates = '22.03.30'
manips, wells, cells = 1, 1, 1
depthoNames = '22.03.30_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.03.30_M1
dates = '22.03.30'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.03.30_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.03.30_M2
dates = '22.03.30'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.03.30_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.03.30_M3
dates = '22.03.30'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.03.30_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.03.30_M4
dates = '22.03.30'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.03.30_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %% Mechanics & Non-linearity

# %%% Next experiment day
# %%%% Next manipe

# %%% 22.03.21, >>>> SINUS <<<< on 3T3 LG+++, M3, M450, various freq and ampli
# %%%% Test !
dates = '22.03.21'
manips, wells, cells = 3, 1, 1
depthoNames = '22.03.21_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 22.03.21, >>>> BROKEN RAMP <<<< on 3T3 LG+++, M4, M450, 3 comp
# %%%% Test !
dates = '22.03.21'
manips, wells, cells = 4, 1, 3
depthoNames = '22.03.21_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 22.02.09, compressionsLowStart of 3T3, M1 = M450, M2 = M450
# %%%% 22.02.09_M1 C1 Seulement
dates = '22.02.09'
manips, wells, cells = 1, 1, 1
depthoNames = '22.02.09_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.02.09_M1
dates = '22.02.09'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.02.09_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.01.12_M2
dates = '22.02.09'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.02.09_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# C3 and C6 are shit xy analysis due to chainof beads plus X motion -> corrected, allowed me to do a LOT of debugging :)

# %%%% 22.01.12 _ Only a few cells
dates = '22.02.09'
manips, wells, cells = 1, 1, [7]
depthoNames = '22.02.09_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# trackAll = False works even for NB = 4 -> Corrected, I think !
# trackAll = True seems to work for NB = 2 but do not when NB = 4



# %%%% 22.01.12_M3
dates = '22.02.09'
manips, wells, cells = 3, 1, 'all'
depthoNames = '22.02.09_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%% 22.01.12, compressionsLowStart of 3T3, M1 = M270, M2 = M450, M4 = M450, pas de M3
# %%%% 22.01.12_M1 C1 Seulement
dates = '22.01.12'
manips, wells, cells = 1, 1, 1
depthoNames = '22.01.12_M1_M270_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.01.12_M1
dates = '22.01.12'
manips, wells, cells = 1, 1, 'all'
depthoNames = '22.01.12_M1_M270_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.01.12_M2
dates = '22.01.12'
manips, wells, cells = 2, 1, 'all'
depthoNames = '22.01.12_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%%% 22.01.12_M4
dates = '22.01.12'
manips, wells, cells = 4, 1, 'all'
depthoNames = '22.01.12_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.12.16, compressions of 3T3, M1 = M450, M2 = M270
# %%%% 21.12.16_M1 C1 Seulement
dates = '21.12.16'
manips, wells, cells = 1, 1, 9
depthoNames = '21.12.16_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%%% 21.12.16_M1
dates = '21.12.16'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.12.16_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.12.16_M2
dates = '21.12.16'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.12.16_M2_M270_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%% 21.12.08, compressions of 3T3, M1 = M270, M2 = M450
# %%%% 21.12.08_M1 C1 Seulement
dates = '21.12.08'
manips, wells, cells = 1, 2, 4
depthoNames = '21.12.08_M1_M270_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%%% 21.12.08_M1
dates = '21.12.08'
manips, wells, cells = 1, 2, 'all'
depthoNames = '21.12.08_M1_M270_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.12.08_M2
dates = '21.12.08'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.12.08_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%% 21.10.25, compressions of 3T3, M1 = M450, M2 = M270
# %%%% 21.10.25_M1 C1 Seulement
dates = '21.10.25'
manips, wells, cells = 1, 1, 1
depthoNames = '21.10.25_M1_M450_100X_step20'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.10.25_M1
dates = '21.10.25'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.10.25_M1_M450_100X_step20'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.10.25_M2
dates = '21.10.25'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.10.25_M2_M270_100X_step20'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.10.18, compressions of 3T3, M1 = M270, M2 = M450
# %%%% 21.10.18_M1
dates = '21.10.18'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.10.18_M1_M270_100X_step20'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.10.18_M2
dates = '21.10.18'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.10.18_M2_M450_100X_step20'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')




# %% MCA project

# %%% Next experiment day
# %%%% Next manipe

# %%% 22.07.27, compressions of 3T3aSFL, F8-3 vs E4-3, M450, M1 = E4-doxy, M2 = E4-control, M3 = F8-doxy, M4 = F8-control
# %%%% One cell

# %%%% 22.07.27_M1

# %%%% 22.07.275_M2

# %%%% 22.07.27_M3

# %%%% 22.07.27_M4

# %%% 22.07.20, compressions of 3T3aSFL, E4-3 vs A11-3, M450, M1 = A11-doxy, M2 = A11-control, M1 = E4-doxy, M2 = E4-control
# %%%% One cell

# %%%% 22.07.20_M1

# %%%% 22.07.20_M2

# %%%% 22.07.20_M3

# %%%% 22.07.20_M4

# %%% 22.07.15, compressions of 3T3aSFL, F8-3 vs A11-3, M450, M1 = F8-doxy, M2 = F8-control, M3 = A11-doxy, M4 = A11-control
# %%%% One cell
dates = '22.07.15'
manips, wells, cells = 1, 1, 1
depthoNames = '22.07.15_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = True, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 22.07.15_M1

# %%%% 22.07.15_M2

# %%%% 22.07.15_M3

# %%%% 22.07.15_M4

# %%% 21.09.09, compressions of 3T3aSFL-A8-2, M450, M1 = doxy, M2 = control
# %%%% 21.09.09_M1
dates = '21.09.09'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.09.09_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.09_M2
dates = '21.09.09'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.09.09_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.09.08, compressions of 3T3aSFL-6FP-2, M450, M1 & M4 = doxy, M2 & M3 = control
# %%%% 21.09.08_M1
dates = '21.09.08'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.09.08_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.08_M2
dates = '21.09.08'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.09.08_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.08_M3
dates = '21.09.08'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.09.08_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.08_M4
dates = '21.09.08'
manips, wells, cells = 4, 1, 'all'
depthoNames = '21.09.08_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default') 



# %%% 21.09.02, compressions of 3T3aSFL-6FP-2, M450, M1 = smifh2, M2 = dmso, M3 = smifh2+doxycyclin, M4 = dmso+doxycyclin
# %%%% 21.09.02_M1
dates = '21.09.02'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.09.02_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.02_M2
dates = '21.09.02'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.09.02_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.02_M3
dates = '21.09.02'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.09.02_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.02_M4
dates = '21.09.02'
manips, wells, cells = 4, 1, 'all'
depthoNames = '21.09.02_M4_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.09.01, compressions of 3T3aSFL with drugs, M450, M1 = dmso, M2 = smifh2
# %%%% 21.09.01_M1
dates = '21.09.01'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.09.01_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.09.01_M2
dates = '21.09.01'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.09.01_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.04.28, constant field of 3T3aSFL-6FP, M450, M1 = doxy, M2 = control
# %%%% 21.04.28_M1

dates = '21.04.28'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.04.28_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.04.28_M2

dates = '21.04.28'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.04.28_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.04.27, constant field of 3T3aSFL-6FP, M450, M1 = control, M2 = doxy
# %%%% 21.04.27_M1

dates = '21.04.27'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.04.27_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.04.27_M2

dates = '21.04.27'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.04.27_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.04.23, constant field of 3T3aSFL, M450, M1 = doxy, M2 = control
# %%%% 21.04.23_M1

dates = '21.04.23'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.04.23_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.04.23_M2

dates = '21.04.23'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.04.23_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.04.21, constant field of 3T3aSFL, M450, M1 = control, M2 = doxy
# %%%% 21.04.21_M1

dates = '21.04.21'
manips, wells, cells = 1, 1, [1, 2, 3, 4, 5]
depthoNames = '21.04.21_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.04.21_M2

dates = '21.04.21'
manips, wells, cells = 2, 1, [1, 3, 6, 8]
depthoNames = '21.04.21_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.02.15, constant field of 3T3aSFL, M450, M1 = control, M2 = doxy

# %%%% 21.02.15_M1_C3

dates = '21.02.15'
manips, wells, cells = 1, 1, 3
depthoNames = '21.02.15_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.02.15_M1

dates = '21.02.15'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.02.15_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.02.15_M2

dates = '21.02.15'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.02.15_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.02.15_M3

dates = '21.02.15'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.02.15_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.02.10, constant field of 3T3aSFL, M450, M1 = control, M2 = doxy

# %%%% 21.02.10_M1_C1

dates = '21.02.10'
manips, wells, cells = 1, 1, 1
depthoNames = '21.02.10_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.02.10_M1

dates = '21.02.10'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.02.10_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.02.10_M2

dates = '21.02.10'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.02.10_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')



# %%% 21.01.21, compressions of 3T3aSFL, M450, M1 = doxy, M2 = control, M3 = doxy

# %%%% 21.01.21_M1 - juste une cellule pour commencer
dates = '21.01.21'
manips, wells, cells = 2, 1, [9]
depthoNames = '21.01.21_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')


# %%%% 21.01.21_M1
dates = '21.01.21'
manips, wells, cells = 1, 1, 'all'
depthoNames = '21.01.21_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.01.21_M2
dates = '21.01.21'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.01.21_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.01.21_M3
dates = '21.01.21'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.01.21_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%% 21.01.18, compressions of 3T3aSFL, M450, M1 = control, M2 = doxy, M3 = control
# %%%% 21.01.18_M1
dates = '21.01.18'
manips, wells, cells = ['1-1', '1-2'], 1, 'all'
depthoNames = '21.01.18_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.01.18_M1
dates = '21.01.18'
manips, wells, cells = ['1-2'], 1, 5
depthoNames = '21.01.18_M1_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.01.18_M2
dates = '21.01.18'
manips, wells, cells = 2, 1, 'all'
depthoNames = '21.01.18_M2_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')

# %%%% 21.01.18_M3
dates = '21.01.18'
manips, wells, cells = 3, 1, 'all'
depthoNames = '21.01.18_M3_M450_step20_100X'

output = mainTracker(dates, manips, wells, cells, depthoNames, expDf, 
                     redoAllSteps = False, MatlabStyle = True, trackAll = False, 
                     sourceField = 'default')



















