# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 19:57:04 2017

@author: Jack
"""

import numpy as np
import random

#initialize parameters
Lrev = 10
Hrev = 10
res = 1
L = Lrev/res
H = Hrev/res

rev = np.zeros((H,L ))
path = np.zeros((H,L))
aggregatePaths =  np.zeros((H,L))

random.seed(a=1)


#initialize checkerboard of percentages
for i in range(1,L+1):
    for j in range(1,H+1):
        if (j % 2 == 0) and (i % 2 == 0):
            rev[i-1,j-1]=random.uniform(0,1)
       
        if (j % 2 <> 0) and (i % 2 <> 0):
            rev[i-1,j-1]=random.uniform(0,1)         
    
np.savetxt('rev.csv', rev)




#generate path

#generate starting spot

y = random.randint(0,H)


#check percentage of two diagonal squares to the 'right'
for x in range (0,H-1):
    path[y,x]= 1
    if y>0 and y<H:
        up = rev[y+1, x+1]
        down = rev[y-1, x-1]
        diff = up-down
        if random.uniform(-1,1) < diff:
            y= y+1
        else:
            y= y-1
    
    elif y == 0:
        y=y+1
    elif y == H:
        y=y-1
path[y,x+1]= 1
#record where you went
np.savetxt('path.csv', path)


