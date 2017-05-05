# -*- coding: utf-8 -*-
"""
Created on Wed May 03 20:30:13 2017

@author: Jack
"""

import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colorbar
import matplotlib.colors



#######From user ImportanceOfBeingErnest on stack overflow#################################################
def cuboid_data(center, size=(1,1,1)):
    # code taken from
    # http://stackoverflow.com/questions/30715083/python-plotting-a-wireframe-3d-cuboid?noredirect=1&lq=1
    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    o = [a - b / 2 for a, b in zip(center, size)]
    # get the length, width, and height
    l, w, h = size
    x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in bottom surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in upper surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  # x coordinate of points in outside surface
         [o[0], o[0] + l, o[0] + l, o[0], o[0]]]  # x coordinate of points in inside surface
    y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in bottom surface
         [o[1], o[1], o[1] + w, o[1] + w, o[1]],  # y coordinate of points in upper surface
         [o[1], o[1], o[1], o[1], o[1]],          # y coordinate of points in outside surface
         [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]    # y coordinate of points in inside surface
    z = [[o[2], o[2], o[2], o[2], o[2]],                        # z coordinate of points in bottom surface
         [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],    # z coordinate of points in upper surface
         [o[2], o[2], o[2] + h, o[2] + h, o[2]],                # z coordinate of points in outside surface
         [o[2], o[2], o[2] + h, o[2] + h, o[2]]]                # z coordinate of points in inside surface
    return x, y, z

def plotCubeAt(pos=(0,0,0), c="b", alpha=0.1, ax=None):
    # Plotting N cube elements at position pos
    if ax !=None:
        X, Y, Z = cuboid_data( (pos[0],pos[1],pos[2]) )
        ax.plot_surface(X, Y, Z, color=c, rstride=1, cstride=1, alpha=0.1)

def plotMatrix(ax, x, y, z, data, cmap="jet", cax=None, alpha=0.1):
    # plot a Matrix 
    norm = matplotlib.colors.Normalize(vmin=data.min(), vmax=data.max())
    colors = lambda i,j,k : matplotlib.cm.ScalarMappable(norm=norm,cmap = cmap).to_rgba(data[i,j,k]) 
    for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                for k, zi, in enumerate(z):
                    plotCubeAt(pos=(xi, yi, zi), c=colors(i,j,k), alpha=alpha,  ax=ax)



    if cax !=None:
        cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')  
        cbar.set_ticks(np.unique(data))
        # set the colorbar transparent as well
        cbar.solids.set(alpha=alpha)              

#######################################################################################################


#main


append = 0


L= 5
H= 5
W= 5

earth = np.random.rand(L,W,H)
path = np.zeros((L,W,H))
aggregatePaths =  np.zeros((L,W,H))


def advance_step(l,w,h,t):
    global L
    global H
    global W
    global earth
    global path
    global aggregatePaths
    
    north =0
    south =0
    east =0
    west = 0
    down = 0
    if  l != 0 and w!= 0 and w != W-1 and l!= L-1 and h!=0:
                north = earth[l+1,w,h]
                south =earth[l-1,w,h]
                east = earth[l,w+1,h]
                west =earth[l,w-1,h]
                down = earth[l,w,h-1] #error
                print 1
                nextIndex =  max(north, south, east, west, down)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
                elif nextIndex == south:
                    path[l-1,w,h] += 1
                    l=l-1
                elif nextIndex == east:
                    path[l,w+1,h] += 1
                    w=w+1
                elif nextIndex == west:
                    path[l,w-1,h] += 1
                    w= w-1
                elif nextIndex == down:
                    path[l,w,h-1] += 1
                    h=h-1
              ##check corneres
             
              #l=0,w=0
    elif  l == 0 and w== 0 and h!=0:
                north = earth[l+1,w,h]
               
                east = earth[l,w+1,h]
         
                down = earth[l,w,h-1]
                print 2
                nextIndex =  max(north,  east,  down)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
              
                elif nextIndex == east:
                    path[l,w+1,h] += 1
                    w=w+1
              
                elif nextIndex == down:
                    path[l,w,h-1] += 1
                    h=h-1
              #l=0,w=w
    elif  l == 0  and w == W-1 and h!=0:
                north = earth[l+1,w,h]
             
                print 3
                west =earth[l,w-1,h]
                down = earth[l,w,h-1]
                
                nextIndex =  max(north,  west, down)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
                
                elif nextIndex == west:
                    path[l,w-1,h] += 1
                    w= w-1
                elif nextIndex == down:
                    path[l,w,h-1] = 1
                    h=h-1
              #l=l,w=0
    elif   w== 0 and l== L-1 and h!=0:
                print 4
                south =earth[l-1,w,h]
                east = earth[l,w+1,h]
       
                down = earth[l,w,h-1]
                
                nextIndex =  max(south, east, down)
                
                if nextIndex == south:
                    path[l-1,w,h] += 1
                    l=l-1
                elif nextIndex == east:
                    path[l,w+1,h] += 1
                    w=w+1
               
                elif nextIndex == down:
                    path[l,w,h-1] += 1
                    h=h-1
              #l=l,w=w
    elif w == W-1 and l== L-1 and h!=0:
                print 5
                south =earth[l-1,w,h]
              
                west =earth[l,w-1,h]
                down = earth[l,w,h-1]
                
                nextIndex =  max( south,  west, down)
                
                if nextIndex == south:
                    path[l-1,w,h]+= 1
                    l=l-1
                
                elif nextIndex == west:
                    path[l,w-1,h] += 1
                    w= w-1
                elif nextIndex == down:
                    path[l,w,h-1] += 1
                    h=h-1
               #l=0,w=0, h=0
    elif  l == 0 and w== 0 and h==0 :
                north = earth[l+1,w,h]
              
                east = earth[l,w+1,h]
                print 6
             
                
                nextIndex =  max(north,  east)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
              
                elif nextIndex == east:
                    path[l,w+1,h]+= 1
                    w=w+1
               
                  
              #l=0,w=w,h=0
    elif  l == 0 and w == W-1 and h==0:
                north = earth[l+1,w,h]
            
                print 7
                west =earth[l,w-1,h]
               
                
                nextIndex =  max(north,  west)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
            
                elif nextIndex == west:
                    path[l,w-1,h]+= 1
                    w= w-1
             
              #l=l,w=0,h=0
    elif  w== 0  and l== L-1 and h==0:
                print 8
                south =earth[l-1,w,h]
                east = earth[l,w+1,h]
                
          
                
                nextIndex =  max(south, east)
           
                if nextIndex == south:
                    path[l-1,w,h] += 1
                    l=l-1
                elif nextIndex == east:
                    path[l,w+1,h] += 1
                    w=w+1
                
              #l=l,w=w,h=0
    elif   w == W-1 and l== L-1 and h==0:
                print 9
                south =earth[l-1,w,h]
               
                west =earth[l,w-1,h]
              
                
                nextIndex =  max( south,  west)
               
                if nextIndex == south:
                    path[l-1,w,h]+= 1
                    l=l-1
               
                elif nextIndex == west:
                    path[l,w-1,h] += 1
                    w= w-1
    elif   h==0 and l != 0 and w!= 0 and w != W-1 and l!= L-1:
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
               
                    print 14
                    nextIndex =  max(north, south, east, west)
                    if nextIndex == north:
                        path[l+1,w,h]+= 1
                        l=l+1
                    elif nextIndex == south:
                        path[l-1,w,h] += 1
                        l=l-1
                    elif nextIndex == east:
                        path[l,w+1,h] += 1
                        w=w+1
                    elif nextIndex == west:
                        path[l,w-1,h] += 1
                        w= w-1
              #check edges
    elif  l == 0 :
                    north = earth[l+1,w,h]
                    print 10
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       down = earth[l,w,h-1]
                    nextIndex =  max(north,  east, west, down)
                    if nextIndex == north:
                        path[l+1,w,h]+= 1
                        l=l+1
               
                    elif nextIndex == east:
                        path[l,w+1,h]+= 1
                        w=w+1
                    elif nextIndex == west:
                        path[l,w-1,h] += 1
                        w= w-1
                    elif nextIndex == down:
                        path[l,w,h-1]+= 1
                        h=h-1
                    
    elif  w== 0 :
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    print 11
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       down = earth[l,w,h-1]
                
                    nextIndex =  max(north, south, east,  down)
                    if nextIndex == north:
                        path[l+1,w,h] += 1
                        l=l+1
                    elif nextIndex == south:
                        path[l-1,w,h]+= 1
                        l=l-1
                    elif nextIndex == east:
                        path[l,w+1,h] += 1
                        w=w+1
               
                    elif nextIndex == down:
                        path[l,w,h-1]+= 1
                        h=h-1
                    
    elif   w == W-1:
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                    print 12
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       down = earth[l,w,h-1]
                
                    nextIndex =  max(north, south,  west, down)
                    if nextIndex == north:
                        path[l+1,w,h] += 1
                        l=l+1
                    elif nextIndex == south:
                        path[l-1,w,h] += 1
                        l=l-1
               
                    elif nextIndex == west:
                        path[l,w-1,h] += 1
                        w= w-1
                    elif nextIndex == down:
                        path[l,w,h-1] += 1
                        h=h-1
                    
    elif   l== L-1 :
      
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       down = earth[l,w,h-1]

                    print 13
                    nextIndex =  max( south, east, west, down)
             
                    if nextIndex == south:
                        path[l-1,w,h] += 1
                        l=l-1
                    elif nextIndex == east:
                        path[l,w+1,h] += 1
                        w=w+1
                    elif nextIndex == west:
                        path[l,w-1,h] += 1
                        w= w-1
                    elif nextIndex == down:
                        path[l,w,h-1]+= 1
                        h=h-1
                    
          
            #need corneres
        
    if path[l,w,h]>1:
                path[l,w,h]=path[l,w,h]-1
                path[l,w,h-1]+=1  #eeror
                print 'doubled back', l ,w, h
                if nextIndex == north:
                    path[l-1,w,h] += -1
                  
                elif nextIndex == south:
                    path[l+1,w,h] += -1
                
                elif nextIndex == east:
                    path[l,w-1,h] += -1
                   
                elif nextIndex == west:
                    path[l,w+1,h] += -1
                   
                elif nextIndex == down:
                
                    path[l,w,h+1] += -1
                
                t=t-1
                if(h!=0):
                    h=h-1
                else:
                    h=h
                    earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                    path =  np.concatenate((np.zeros((L,W,1)),path),2)
                    aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                    append+=1
                path[l,w,h]+=1
    return [l,w,h,t]




timesteps =10

#simulate rainfall event, run a trace from everywhere on the land surfac
for i in range(L-2,L-1):
    for j in range(W-2,W-1):
        h=H-1
        l=i
        w=j
        path[l,w,h]=1
        for t in range (0,timesteps):
           #central points
            print  l,w,h
           
            l,w,h,t =      advance_step(l,w,h,t) 
                
        aggregatePaths =  aggregatePaths + path
        print np.sum(path)
    
        
        path = np.zeros((L,W,H))    

#print np.sum(aggregatePaths)/(timesteps*L*W)
   

X=range(0,L)
Y=range(0,W)
Z=range (0,H+append)

fig = plt.figure(figsize=(10,4))
ax = fig.add_axes([0.1, 0.1, 0.7, 0.8], projection='3d')
ax_cb = fig.add_axes([0.8, 0.3, 0.05, 0.45])
ax.set_aspect('equal')
plotMatrix(ax, X, Y, Z, aggregatePaths, cmap="jet", cax = ax_cb)

plt.savefig('plt'+".png")
plt.show()

