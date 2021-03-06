# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 19:57:04 2017

@author: Jack
"""

import numpy as np
import random
import matplotlib.pyplot as plt


#initialize parameters
Lrev = 100
Hrev = 100
res = 1
L = Lrev/res
H = Hrev/res

rev = np.zeros((H,L ))
path = np.zeros((H,L))
aggregatePaths =  np.zeros((H,L))

random.seed(7983)

#initialize random percentages
for i in range(1,L+1):
    for j in range(1,H+1):
#        if (j % 2 == 0) and (i % 2 == 0):
            rev[i-1,j-1]=random.uniform(0,1)
       
#        if (j % 2 <> 0) and (i % 2 <> 0):
#            rev[i-1,j-1]=random.uniform(0,1)   



#turn k into a layer
#introduce a permeable layer 
#for i in range (0, L/2):
#    rev[i-1,:]=0.01
##introduce impermeable layer 
#for i in range (L/2, L+1):
#    rev[i-1,:]=0.99
    
np.savetxt('rev.csv', rev)




#generate path

#generate starting spot


trials=10000
#check percentage of two diagonal squares to the 'right'
for trial in range (0,trials) :

    random.seed(a=trial+6)
    y = L/2 #start from land surface
    for x in range (0,L-1):
        path[y,x]= 1
        if y>0 and y<(H-1):
            up = rev[y-1, x+1]
            down = rev[y+1, x+1]
            
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


    aggregatePaths = aggregatePaths + path
    path = np.zeros((H,L))

np.savetxt('aggregatePath.csv', aggregatePaths)
plt.imshow(aggregatePaths, cmap='hot', interpolation='nearest')
plt.show()