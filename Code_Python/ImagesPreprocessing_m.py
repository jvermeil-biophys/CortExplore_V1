# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:21:10 2022

@author: anumi

source : https://docs.opencv.org/3.4/db/d5b/tutorial_py_mouse_handling.html
"""

#### Imports

import os
import cv2
import shutil

import numpy as np
import pyjokes as pj

from skimage import io

#### Local Imports

import sys
import CortexPaths as cp
sys.path.append(cp.DirRepoPython)

import GraphicStyles as gs
import GlobalConstants as gc
import UtilityFunctions as ufun

#%% Define parameters # Numi

DirSave = 'D:/Anumita/MagneticPincherData/Raw/22.06.09'
DirExt = 'F:/Cortex Experiments/OptoPincher Experiments/20220906_100xoil_3t3optorhoa_4.5beads_15mT/22.06.09/'
prefix = 'cell'
channel = 'w1TIRF DIC'
microscope = 'metamorph'

#%% Define parameters # Jojo

date = '22.05.03'
DirExt = 'E:/22.05.03_HoxB8/M4_patterns_ctrl'
DirSave = os.path.join(cp.DirDataRaw, date)

# prefix = 'cell'
# channel = 'w1TIRF DIC'
microscope = 'labview'


# %% Functions


def removeFrames(DirSave):
    allFiles = os.listdir(DirSave)
    for f in allFiles:
        if f.endswith('.tif'):
            print('Loading..'+f)
            ic = io.ImageCollection(DirSave+'/'+f, conserve_memory = True)
            stack = io.concatenate_images(ic[:4500])
            print('Saving...'+f)
            io.imsave(DirSave+'/'+f, stack)

def AllMMTriplets2Stack(DirExt, DirSave, prefix, channel):
    allCells = os.listdir(DirExt)
    for currentCell in allCells:
        path = os.path.join(DirExt, currentCell)
        allFiles = os.listdir(path)
        date = ufun.findInfosInFileName(currentCell, 'date')
        allFiles = [path+'/'+string for string in allFiles if channel in string]
        #+4 at the end corrosponds to the '_t' part to sort the array well
        limiter = len(path)+len(prefix)+len(channel)+4 
        allFiles.sort(key=lambda x: int(x[limiter:-4]))
        ic = io.ImageCollection(allFiles, conserve_memory = True)
        stack = io.concatenate_images(ic)
        if '561' in channel:
            date = date.replace('-', '.')
            folderName = date+'_561'
            fluoPath = os.path.join(DirSave, folderName)
            try:
                os.mkdir(fluoPath)
            except:
                pass
            
        io.imsave(fluoPath+'/'+currentCell+'.tif', stack)
        
def renamePrefix(DirExt, currentCell, newPrefix):
    path = os.path.join(DirExt, currentCell)
    allImages = os.listdir(path)
    for i in allImages:
        if i.endswith('.TIF'):
            split = i.split('_')
            split[0] = newPrefix
            newName = '_'.join()
            os.rename(os.path.join(path,i), os.path.join(path, newName))

def Zprojection(currentCell, microscope, kind = 'min'):
    scaleFactor = 4
    path = os.path.join(DirExt, currentCell)
    allFiles = os.listdir(path)
    
    if microscope == 'metamorph':
        allFiles = [path+'/'+string for string in allFiles if channel in string]
        #+4 at the end corrosponds to the '_t' part to sort the array well
        limiter = len(path)+len(prefix)+len(channel)+4 
        allFiles.sort(key=lambda x: int(x[limiter:-4]))
    
    elif microscope == 'labview':
        allFiles = [path+'/'+string for string in allFiles if 'im' in string]
        
    idx = slice(0, len(allFiles), 100)
    
    allFiles = allFiles[idx]
    frame = cv2.imread(allFiles[0])
    imgWidth, imgHeight = frame.shape[1], frame.shape[0]
    ic = io.ImageCollection(allFiles, conserve_memory=True)
    stack = io.concatenate_images(ic)
        
    if kind == 'min':
        Zimg = np.min(stack, axis = 0)
    elif kind == 'max':
        Zimg = np.max(ic, axis = 0)
        
    Zimg = cv2.resize(Zimg, (int(imgWidth/scaleFactor), int(imgHeight/scaleFactor)))
    Zimg = cv2.gs.NORMALize(Zimg, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    return Zimg

def shape_selection(event, x, y, flags, param):
    # grab references to the global variables
    global ref_point, crop, allZimg, iZ
    global ref_point, crop

    # if the left mouse button was clicked, record the starting
    # (x, y) coordinates and indicate that cropping is being performed
    if event == cv2.EVENT_LBUTTONDOWN:
        ref_point = [(x, y)]

    # check to see if the left mouse button was released
    elif event == cv2.EVENT_LBUTTONUP:
        # record the ending (x, y) coordinates and indicate that
        # the cropping operation is finished
        ref_point.append((x, y))

        # draw a rectangle around the region of interest
        cv2.rectangle(allZimg[i], ref_point[0], ref_point[1], (0, 255, 0), 1)

def crop(DirSave, allRefPoints, allCells, microscope):
    count = 0
    for i,j in zip(allRefPoints, allCells):
        path = os.path.join(DirExt, j)
        allFiles = os.listdir(path)
        
        if microscope == 'metamorph':
            allFiles = [path+'/'+string for string in allFiles if channel in string]
            #+4 at the end corrosponds to the '_t' part to sort the array well
            limiter = len(path)+len(prefix)+len(channel)+4 
            allFiles.sort(key=lambda x: int(x[limiter:-4]))
    
        elif microscope == 'labview':
            allFiles = [path+'/'+string for string in allFiles if 'im' in string]
        
        print(gs.BLUE + 'Loading '+j+'...' + gs.NORMAL)
        ic = io.ImageCollection(allFiles, conserve_memory = True)
        stack = io.concatenate_images(ic)
        
        x1, x2, y1, y2 = int(i[0][0]), int(i[1][0]), int(i[0][1]), int(i[1][1])
        
        # To avoid that the cropped region gets bigger than the image itself
        ny, nx = stack.shape[1], stack.shape[2]
        x1, x2, y1, y2 = max(0, x1), min(nx, x2), max(0, y1), min(ny, y2)
        
        cropped = stack[:, y1:y2, x1:x2]
        
        io.imsave(DirSave+'/'+j+'.tif', cropped)
        print(gs.GREEN + j +' saved sucessfully' + gs.NORMAL)
        if count%5 == 0:
            joke = pj.get_joke(language='en', category= 'all')
            print(joke)
            count = count + 1

def moveFiles(DirSave, allCells, filename):
    for i in allCells:
        path = os.path.join(DirExt, i)
        allFiles = os.listdir(path)
        for f in allFiles:
            if filename in f:
                source = os.path.join(path, f)
                destination = DirSave+f+'.txt'
                shutil.copy(source, destination)
                break




# preprocess(DirExt, DirSave, microscope, reset = 0)    

#%% Main function

# def preprocess(DirExt, DirSave, microscope, reset = 0):
    
allCells = os.listdir(DirExt)
ref_point = []
allRefPoints = []
allZimg = []
allZimg_og = []
reset = 0

try:
    os.mkdir(DirSave)
except:
    pass

print(gs.BLUE + 'Constructing all Z-Projections...' + gs.NORMAL)

scaleFactor = 4
invalidCellIndex = []
for i in range(len(allCells)):
    currentCell = allCells[i]
    print(currentCell)
    try:
        Zimg = Zprojection(currentCell, microscope)
        allZimg.append(Zimg)
    except:
        print(currentCell + ' is not a valid cell')
        invalidCellIndex.append(i)

allCells2 = []
for i in range(len(allCells)):
    if not i in invalidCellIndex:
        allCells2.append(allCells[i])
        
allCells = allCells2

allZimg_og = np.copy(np.asarray(allZimg))


print(gs.ORANGE + 'Draw the ROIs to crop...' + gs.NORMAL)

if reset == 1:
    allZimg = np.copy(allZimg_og)
    ref_point = []
    allRefPoints = []

count = 0
for i in range(len(allZimg)):
    if count%25 == 0:
        count = 0
        
    currentCell = allCells[i]
    
    Nimg = len(allZimg)
    ncols = 5
    nrows = ((Nimg-1) // ncols) + 1
    clone = allZimg[i].copy()
    
    cv2.namedWindow(currentCell)
    cv2.moveWindow(currentCell, (count//ncols)*400, count%ncols*200)
    cv2.setMouseCallback(currentCell, shape_selection)
    
    while True:
    # display the image and wait for a keypress
        currentCell = allCells[i]
        cv2.imshow(currentCell, allZimg[i])
        key = cv2.waitKey(20) & 0xFF
        
    # press 'r' to reset the crop
        if key == ord("r"):
            allZimg[i] = clone.copy()    
    
    # if the 'a' key is pressed, break from the loop and move on to the next file
        elif key == ord("a"):
            allRefPoints.append(np.asarray(ref_point)*scaleFactor)
            break
        
    count = count + 1
    
cv2.destroyAllWindows()

#%%
print(gs.BLUE + 'Saving all tiff stacks...' + gs.NORMAL)
crop(DirSave, allRefPoints[17:], allCells[17:], microscope)


#%% Creating .tif stacks of 561n recruitment images

DirSave = 'D:/Anumita/MagneticPincherData/Raw/'
DirExt = 'F:/Cortex Experiments/OptoPincher Experiments/20220322_100xoil_3t3optorhoa_4.5beads_15mT/22.03.22'
prefix = 'cell'
channel = 'w3TIRF 561'


AllMMTriplets2Stack(DirExt, DirSave, prefix, channel)

