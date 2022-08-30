# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:50:53 2021

@author: Joseph Vermeil
"""

# %% General imports

# 1. Imports

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as st
import statsmodels.api as sm
import matplotlib.pyplot as plt

import os
import sys
import matplotlib

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

expDf = ufun.getExperimentalConditions(DirExp = cp.DirRepoExp, save = True, sep = ';')

# %% Small things

#### close all
plt.close('all')



# %% Next Topic !!

# %%% Next experiment day

# %%%% Next manipe

# %%%% Next manipe




# %% Next Topic !!

# %%% Next experiment day

# %%%% Next manipe

# %%%% Next manipe




# %% Example Topic - HoxB8 macrophages

# %%% 22.05.05, compressionsLowStart of HoxB8 macrophages, M450, M1 = tko & glass, M2 = ctrl & glass
# %%%% 22.05.05_M1 C1 only
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




