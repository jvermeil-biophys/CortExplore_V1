# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:21:10 2022

@author: Anumita Jawahar & Joseph Vermeil

source : https://docs.opencv.org/3.4/db/d5b/tutorial_py_mouse_handling.html
"""

#### Imports

import os
import cv2
import shutil
import traceback

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

# date = '22.07.15'
# DirExt = 'G:/22.07.15_longLinker' #'/M4_patterns_ctrl'
# DirSave = os.path.join(cp.DirDataRaw, date)

# prefix = ''
# channel = ''
# microscope = 'labview'


# %% Functions

def getListOfSourceFolders(Dir, forbiddenWords = ['deptho']):
    res = []
    exclude = False
    for w in forbiddenWords:
        if w.lower() in Dir.lower(): # compare the lower case strings
            exclude = True # If a forbidden word is in the dir name, don't consider it
            
            
    if exclude or not os.path.isdir(Dir):
        return(res) # Empty list
    
    elif ufun.containsFilesWithExt(Dir, '.tif'):
        res = [Dir] # List with 1 element - the name of this dir
        
    else:
        listDirs = os.listdir(Dir)
        for D in listDirs:
            path = os.path.join(Dir, D)
            res += getListOfSourceFolders(path) # Recursive call to the function !
    # In the end this function will have explored all the sub directories of Dir,
    # searching for folders containing tif files, without forbidden words in their names.        
    return(res)

# test

# testDir = 'G:/22.07.15_longLinker'
# A = getListOfSourceFolders(testDir)

#

def copyFieldFiles(ListDirSrc, DirDst, suffix = '_Field.txt'):
    for DirSrc in ListDirSrc:
        ufun.copyFilesWithString(DirSrc, DirDst, suffix)


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
    Zimg = cv2.normalize(Zimg, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
    return Zimg

def shape_selection(event, x, y, flags, param):
    # grab references to the global variables
    global ref_point, crop, allZimg #, iZ

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

  
def shape_selection_V2(event, x, y, flags, param):
    """
    Interactive rectangular selection
    """
    # Grab references to the global variables
    global ix, iy, drawing, ref_point, crop, img, img_copy
    
    # If the left mouse button was clicked, record the starting
    # (x, y) coordinates and indicate that cropping is being performed
    if event == cv2.EVENT_LBUTTONDOWN:
        drawing = True
        ix,iy = x,y
        img_copy = np.copy(img)
        ref_point = [(x, y)]
    
    # If the mouse moves, reinitialize the image and the rectangle to match 
    # the current position
    elif event == cv2.EVENT_MOUSEMOVE: 
        
        if drawing == True:
            img = np.copy(img_copy)
            cv2.rectangle(img,(ix,iy),(x,y),(0,255,0),1)

    # check to see if the left mouse button was released
    elif event == cv2.EVENT_LBUTTONUP:
        drawing = False
        # Record the ending (x, y) coordinates and indicate that
        # the cropping operation is finished
        ref_point.append((x, y))
        # Final rectangle around the region of interest
        cv2.rectangle(img,(ix,iy),(x,y),(0,255,0),1)
        

def cropAndCopy(DirSrc, DirDst, allRefPoints, allCellPaths, microscope):
    count = 0
    for refPts, cellPath in zip(allRefPoints, allCellPaths):
        cellName = cellPath.split('\\')[-1]
        allFiles = os.listdir(cellPath)
        
        if microscope == 'metamorph':
            allFiles = [cellPath+'/'+string for string in allFiles if channel in string]
            #+4 at the end corrosponds to the '_t' part to sort the array well
            limiter = len(cellPath)+len(prefix)+len(channel)+4 
            allFiles.sort(key=lambda x: int(x[limiter:-4]))
    
        elif microscope == 'labview':
            allFiles = [cellPath+'/'+string for string in allFiles if 'im' in string]
        
        print(gs.BLUE + 'Loading '+ cellPath +'...' + gs.NORMAL)
        
        ic = io.ImageCollection(allFiles, conserve_memory = True)
        stack = io.concatenate_images(ic)
        
        x1, x2 = int(refPts[0][0]), int(refPts[1][0])
        y1, y2 = int(refPts[0][1]), int(refPts[1][1])
        
        # To avoid that the cropped region gets bigger than the image itself
        ny, nx = stack.shape[1], stack.shape[2]
        x1, x2, y1, y2 = max(0, x1), min(nx, x2), max(0, y1), min(ny, y2)
        
        try:
            cropped = stack[:, y1:y2, x1:x2]
            io.imsave(DirDst + '/' + cellName + '.tif', cropped)
            print(gs.GREEN + DirDst + '/' + cellName + '.tif' + '\nSaved sucessfully' + gs.NORMAL)
        
        except Exception:
            traceback.print_exc()
            print(gs.RED + DirDst + '/' + cellName + '.tif' + '\nError when saving' + gs.NORMAL)
            continue
        
        if count%5 == 0:
            joke = pj.get_joke(language='en', category= 'all')
            print(joke)
            count = count + 1

def moveFilesWithString(DirSave, allCells, filename):
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

#%% Main function 1/3

# def preprocess(DirExt, DirSave, microscope, reset = 0):
    
allCellsRaw = getListOfSourceFolders(DirExt)[:]
allCells = []
# allCellNames = []
ref_point = []
allRefPoints = []
allZimg = []
allZimg_og = []

#%% Main function 2/3

# reset = 0
checkIfAlreadyExist = True

if not os.path.exists(DirSave):
    os.mkdir(DirSave)

print(gs.BLUE + 'Constructing all Z-Projections...' + gs.NORMAL)

scaleFactor = 4

for i in range(len(allCellsRaw)):
    currentCell = allCellsRaw[i]
    currentCellName = currentCell.split('\\')[-1]
    validCell = True
    print(currentCellName)
        
    if not ufun.containsFilesWithExt(currentCell, '.tif'):
        validCell = False
        print('---> Is not a valid cell')
        
    elif checkIfAlreadyExist and os.path.isfile(os.path.join(DirSave, currentCellName + '.tif')):
        validCell = False
        print('---> Has already been copied')
        
    if validCell:
        allCells.append(currentCell)
        Zimg = Zprojection(currentCell, microscope)
        allZimg.append(Zimg)
        
copyFieldFiles(allCells, DirSave)

# allZimg_og = np.copy(np.asarray(allZimg)) # TBC

#%% Main function 3/3

print(gs.ORANGE + 'Draw the ROIs to crop...' + gs.NORMAL)

# if reset == 1:
#     allZimg = np.copy(allZimg_og)
#     ref_point = []
#     allRefPoints = []

count = 0
for i in range(len(allZimg)):
    if count%24 == 0:
        count = 0
        
    currentCell = allCells[i]
    
    Nimg = len(allZimg)
    ncols = 6
    nrows = ((Nimg-1) // ncols) + 1
    
    # test
    ix,iy = 0, 0
    drawing = False
    img = allZimg[i]
    img_backup = np.copy(img)
    img_copy = np.copy(img)
    
    cv2.namedWindow(currentCell)
    cv2.moveWindow(currentCell, (count//ncols)*400, count%ncols*175)
    
    # cv2.setMouseCallback(currentCell, shape_selection)
    cv2.setMouseCallback(currentCell, shape_selection_V2)
    
    while True:
    # display the image and wait for a keypress
        currentCell = allCells[i]
        cv2.imshow(currentCell, img)
        key = cv2.waitKey(20) & 0xFF
        
    # press 'r' to reset the crop
        if key == ord("r"):
            img = np.copy(img_backup)  
             
    # if the 'a' key is pressed, break from the loop and move on to the next file
        elif key == ord("a"):
            allRefPoints.append(np.asarray(ref_point)*scaleFactor)
            break
        
    count = count + 1
    
cv2.destroyAllWindows()

print(gs.BLUE + 'Saving all tiff stacks...' + gs.NORMAL)

cropAndCopy(DirExt, DirSave, allRefPoints[:], allCells[:], microscope)





#%% Creating .tif stacks of 561n recruitment images

# DirSave = 'D:/Anumita/MagneticPincherData/Raw/'
# DirExt = 'F:/Cortex Experiments/OptoPincher Experiments/20220322_100xoil_3t3optorhoa_4.5beads_15mT/22.03.22'
# prefix = 'cell'
# channel = 'w3TIRF 561'


# AllMMTriplets2Stack(DirExt, DirSave, prefix, channel)

