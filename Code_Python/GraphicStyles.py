# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:52:35 2022

@author: JosephVermeil
"""

# %% 0. Imports

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


import itertools
import matplotlib

from cycler import cycler



# %% 1. Markers
my_default_marker_list = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']
markerList10 = ['o', 's', 'D', '>', '^', 'P', 'X', '<', 'v', 'p']

# %% 2. Colors

# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']
my_default_color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                         '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
my_default_color_cycle = cycler(color=my_default_color_list)
plt.rcParams['axes.prop_cycle'] = my_default_color_cycle

pairedPalette = sns.color_palette("tab20")
pairedPalette = pairedPalette.as_hex()
pairedPalette

# clist = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a']
# sns.color_palette(clist)
colorList10 = my_default_color_list
sns.color_palette(my_default_color_list)


bigPalette1 = sns.color_palette("tab20b")
bigPalette1_hex = bigPalette1.as_hex()
bigPalette1

bigPalette2 = sns.color_palette("tab20c")
bigPalette2_hex = bigPalette2.as_hex()
bigPalette2

customPalette_hex = []
for ii in range(2, -1, -1):
    customPalette_hex.append(bigPalette2_hex[4*0 + ii]) # blue
    customPalette_hex.append(bigPalette2_hex[4*1 + ii]) # orange
    customPalette_hex.append(bigPalette2_hex[4*2 + ii]) # green
    customPalette_hex.append(bigPalette1_hex[4*3 + ii]) # red
    customPalette_hex.append(bigPalette2_hex[4*3 + ii]) # purple
    customPalette_hex.append(bigPalette1_hex[4*2 + ii]) # yellow-brown
    customPalette_hex.append(bigPalette1_hex[4*4 + ii]) # pink
    customPalette_hex.append(bigPalette1_hex[4*0 + ii]) # navy    
    customPalette_hex.append(bigPalette1_hex[4*1 + ii]) # yellow-green
    customPalette_hex.append(bigPalette2_hex[4*4 + ii]) # gray
    
# customPalette = sns.color_palette(customPalette_hex)
colorList30 = customPalette_hex

customPalette_hex = []
for ii in range(3, -1, -1):
    customPalette_hex.append(bigPalette2_hex[4*0 + ii]) # blue
    customPalette_hex.append(bigPalette2_hex[4*1 + ii]) # orange
    customPalette_hex.append(bigPalette2_hex[4*2 + ii]) # green
    customPalette_hex.append(bigPalette1_hex[4*3 + ii]) # red
    customPalette_hex.append(bigPalette2_hex[4*3 + ii]) # purple
    customPalette_hex.append(bigPalette1_hex[4*2 + ii]) # yellow-brown
    customPalette_hex.append(bigPalette1_hex[4*4 + ii]) # pink
    customPalette_hex.append(bigPalette1_hex[4*0 + ii]) # navy    
    customPalette_hex.append(bigPalette1_hex[4*1 + ii]) # yellow-green
    customPalette_hex.append(bigPalette2_hex[4*4 + ii]) # gray
    

# customPalette = sns.color_palette(customPalette_hex)
colorList40 = customPalette_hex

# TEST to get a darker list
# colorList40_darker = []
# for i in range(len(colorList40)):
#     c = colorList40[i]
#     print(c)
#     c_darker = jvu.lighten_color(c, 1.25)
#     colorList40_darker.append(c_darker)

# colorList40_darker = colorList40_darker.as_hex()

#### Test the colors

# N = len(my_default_marker_list)
# X = np.arange(1, N+1)
# Y = np.arange(1, N+1)
# fig, ax = plt.subplots(1, 1, figsize = (3, 3))
# for i in range(N):
#     for j in range(N):
#         ax.plot([X[i]], [Y[-1-j]], color = my_default_color_list[i], marker = my_default_marker_list[j], 
#                 ls = '', markersize = 10, markeredgecolor = 'k')
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# plt.show()

def colorTester():
    N = len(colorList40)
    M = len(markerList10)
    X = np.arange(1, N+1)
    Y = np.arange(1, M+1)
    fig, ax = plt.subplots(1, 1, figsize = (0.3*N, 0.3*M))
    for i in range(N):
        for j in range(M):
            ax.plot([X[i]], [Y[-1-j]], color = colorList40[i], marker = markerList10[j], 
                    ls = '', markersize = 10, markeredgecolor = 'k')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.show()

# %% 3. Global constants



#### Console text styles

NORMAL  = '\033[1;0m'
RED  = '\033[0;31m' # red
GREEN = '\033[1;32m' # green
ORANGE  = '\033[0;33m' # orange
BLUE  = '\033[0;36m' # blue
CYAN  = '\033[1;36m' # blue
YELLOW = '\033[1;33m' # yellow
PURPLE = '\033[1;35m' # purple

#### Test the styles

def printTester01():
    print(NORMAL + 'normal' + NORMAL)
    print(RED + 'nothing rhyme with red' + NORMAL)
    print(ORANGE + 'nothing rhyme with orange' + NORMAL)
    print(YELLOW + 'nothing rhyme with yellow' + NORMAL)
    print(GREEN + 'nothing rhyme with green' + NORMAL)
    print(CYAN + 'nothing rhyme with cyan' + NORMAL)
    print(BLUE + 'nothing rhyme with blue' + NORMAL)
    print(PURPLE + 'nothing rhyme with purple' + NORMAL)
    print('\n')


def printTester02():
    print("\033[0;37;48m Normal text\n")
    print("\033[2;37;48m Underlined text\033[0;37;48m \n")
    print("\033[1;37;48m Bright Colour\033[0;37;48m \n")
    print("\033[3;37;48m Negative Colour\033[0;37;48m \n")
    print("\033[5;37;48m Negative Colour\033[0;37;48m\n")
    print("\033[1;37;40m \033[2;37:40m TextColour BlackBackground          TextColour GreyBackground                WhiteText ColouredBackground\033[0;37;40m\n")
    print("\033[1;30;40m Dark Gray      \033[0m 1;30;40m            \033[0;30;47m Black      \033[0m 0;30;47m               \033[0;37;41m Black      \033[0m 0;37;41m")
    print("\033[1;31;40m Bright Red     \033[0m 1;31;40m            \033[0;31;47m Red        \033[0m 0;31;47m               \033[0;37;42m Black      \033[0m 0;37;42m")
    print("\033[1;32;40m Bright Green   \033[0m 1;32;40m            \033[0;32;47m Green      \033[0m 0;32;47m               \033[0;37;43m Black      \033[0m 0;37;43m")
    print("\033[1;33;48m Yellow         \033[0m 1;33;48m            \033[0;33;47m Brown      \033[0m 0;33;47m               \033[0;37;44m Black      \033[0m 0;37;44m")
    print("\033[1;34;40m Bright Blue    \033[0m 1;34;40m            \033[0;34;47m Blue       \033[0m 0;34;47m               \033[0;37;45m Black      \033[0m 0;37;45m")
    print("\033[1;35;40m Bright Magenta \033[0m 1;35;40m            \033[0;35;47m Magenta    \033[0m 0;35;47m               \033[0;37;46m Black      \033[0m 0;37;46m")
    print("\033[1;36;40m Bright Cyan    \033[0m 1;36;40m            \033[0;36;47m Cyan       \033[0m 0;36;47m               \033[0;37;47m Black      \033[0m 0;37;47m")
    print("\033[1;37;40m White          \033[0m 1;37;40m            \033[0;37;40m Light Grey \033[0m 0;37;40m               \033[0;37;48m Black      \033[0m 0;37;48m")
    print("\n")
    
    
    
    
    