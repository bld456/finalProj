## Jack Lange ##
# Simple Water infiltration model #
#5/5/17#


#interpretation:
# 'earth' represents some volume of the land surface and bedrock.
# 'tracer' refers to water following a path of least resistance into the ground
# 'erosion' of these infiltration pathways is depicted by increaseing the probability of pathways taken by the tracer
# The model simulates a rainfall event by following  the path of a tracer from each node on the land surface into the earth. 
# results:
# When the model is run without tracers causing erosion, minimal correlation is seen between the number of visits to a node and the random number assigned to the node.
# when erosion is added, a distinct correlation is seen between node visits and the probability associated with each node 
# when a fracture with probability = 1 is added to the model a vast majority of node visits are visits along the fracture (wher eprobability =1) indicating that it is a major flow conduit. 
# (Fracture is not subject to erosion, though the rest of 'earth' is.)
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colorbar
import matplotlib.colors



#######From user ImportanceOfBeingErnest on stack overflow#################################################
#methods for making  3d plot

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

def plotMatrix(ax, x, y, z, data, cmap="YlOrBr", cax=None, alpha=0.1):
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




#define global variable
append = 0 #number of times the 'earth' has been ammended in the H direction

t=0 #timestep counter
L= 5 #initial earth length
H= 5 #initial earth height
W= 5 #initial earth width

earth = np.random.rand(L,W,H) #generate a cube of 'earth'. Each node in the cube assigned a random number
path = np.zeros((L,W,H)) # a matrix that will track the path of each tracer throug hthe earth
aggregatePaths =  np.zeros((L,W,H)) #a matrix to compile the  paths of each tracer 


def advance_step(l,w,h,t):
#method to advance the tracer one step spatialy. The tracer move one spatial step for every temporal step. The direction which the tracer
# moves is dictated by the random numbers associated with each node surrounding the tracer's current position. The tracer wil move to the
# node with the highest number assiciated with it. 
# the tracer is confined to move north, south, east and west in the L-W plane. The tracer can not leave the L-W boundaries.
# The tracer is able to move down into a lower L-W plane. If the tracer reaches the bottom of the 'earth' a new L-W plane is generated and 
# placed at the bottom of the 'earth'.
# Backtracking is not permitted. If the tracer backtracks it is moved down one plane in the H direction. This can be interpreted as 
# tracer pooling in one position and increasing pressure, allowing the tracer to penetrate into a deeper level.
    global L
    global H
    global W
    global earth
    global path
    global aggregatePaths
    global append

    north =0
    south =0
    east =0
    west = 0
    down = 0
    if  l != 0 and w!= 0 and w != W-1 and l!= L-1 and h!=0: #for interior points in th W-L plane, and points that are not on the bottom of 'earth'
                north = earth[l+1,w,h]
                south =earth[l-1,w,h]
                east = earth[l,w+1,h]
                west =earth[l,w-1,h]
                down = earth[l,w,h-1] #error
          
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
            
                    
    elif  l == 0 and w== 0 and h!=0: # for the corner L=0, W=0
                north = earth[l+1,w,h]
               
                east = earth[l,w+1,h]
         
                down = earth[l,w,h-1]
               
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
            
    elif  l == 0  and w == W-1 and h!=0: # For the corner L=0, W=W
                north = earth[l+1,w,h]
             
                
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
             
    elif   w== 0 and l== L-1 and h!=0:# For the corner L=L, W=0
               
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
             
    elif w == W-1 and l== L-1 and h!=0: # for the corner L=L, W=W
           
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
  
    elif  l == 0 and w== 0 and h==0 :              # for the bottom corner, l=0,w=0, h=0
                north = earth[l+1,w,h]
              
                east = earth[l,w+1,h]
                #need to append a new layer on the bottom because the tracer has reached the bottom
                earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                append+=1
                h+=1
                down = earth[l,w,h-1]
             
                
                nextIndex =  max(north,  east)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
              
                elif nextIndex == east:
                    path[l,w+1,h]+= 1
                    w=w+1
                elif nextIndex == down:
                    path[l,w,h-1]+= 1
                    h=h-1
                  
         
    elif  l == 0 and w == W-1 and h==0: #for the bottom corner l=0,w=w,h=0
                north = earth[l+1,w,h]
            
                
                west =earth[l,w-1,h]
                 #need to append a new layer on the bottom because the tracer has reached the bottom
                earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                append+=1
                h+=1
                down = earth[l,w,h-1]
                
                nextIndex =  max(north,  west,down)
                if nextIndex == north:
                    path[l+1,w,h] += 1
                    l=l+1
            
                elif nextIndex == west:
                    path[l,w-1,h]+= 1
                    w= w-1
                elif nextIndex == down:
                    path[l,w,h-1]+= 1
                    h=h-1
                  
             
          
    elif  w== 0  and l== L-1 and h==0: #for the corner l=l,w=0,h=0
             
                south =earth[l-1,w,h]
                east = earth[l,w+1,h]
                 #need to append a new layer on the bottom because the tracer has reached the bottom
                earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                append+=1
                h+=1
                down = earth[l,w,h-1]
          
                
                nextIndex =  max(south, east,down)
           
                if nextIndex == south:
                    path[l-1,w,h] += 1
                    l=l-1
                elif nextIndex == east:
                    path[l,w+1,h] += 1
                    w=w+1
                elif nextIndex == down:
                    path[l,w,h-1]+= 1
                    h=h-1
             
    elif   w == W-1 and l== L-1 and h==0:    #l=l,w=w,h=0
               
                south =earth[l-1,w,h]
               
                west =earth[l,w-1,h]
                #need to append a new layer on the bottom because the tracer has reached the bottom
                earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                append+=1
                h+=1
                down = earth[l,w,h-1]
                
                nextIndex =  max( south,  west, down)
               
                if nextIndex == south:
                    path[l-1,w,h]+= 1
                    l=l-1
               
                elif nextIndex == west:
                    path[l,w-1,h] += 1
                    w= w-1
                elif nextIndex == down:
                    path[l,w,h-1]+= 1
                    h=h-1
    elif   h==0 and l != 0 and w!= 0 and w != W-1 and l!= L-1: #case of interior, bottom 
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
                   #need to append a new layer on the bottom because the tracer has reached the bottom
                    earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                    path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                    aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                    append+=1
                    h+=1
                    down = earth[l,w,h-1]
                    
                    nextIndex =  max(north, south, east, west, down)
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
                    elif nextIndex == down:
                        path[l,w,h-1]+= 1
                        h=h-1
              
    elif  l == 0 :                          #L=0 edge
                    north = earth[l+1,w,h]
                
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       append+=1
                       h+=1
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
                        if h>0:
                            h=h-1
                  
                             
                                
    elif  w== 0 :                               #W=0 edge
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       append+=1
                       h+=1
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
                        if h>0:
                            h=h-1
                    
    elif   w == W-1:                            #W=w edge
                    north = earth[l+1,w,h]
                    south =earth[l-1,w,h]
                 
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2) 
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       append+=1
                       h+=1
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
                        if h>0:
                            h=h-1
                    
    elif   l== L-1 :                    #L=L edge
      
                    south =earth[l-1,w,h]
                    east = earth[l,w+1,h]
                    west =earth[l,w-1,h]
                    if h!=0:
                        down = earth[l,w,h-1] #error
                    else:
                       aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                       earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                       path =  np.concatenate((np.zeros((L,W,1)),path),2)
                       append+=1
                       h+=1
                      
                       down = earth[l,w,h-1]

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
                        if h>0:
                            h=h-1
                    
          
      
        
    if path[l,w,h]>1:  #Redirect back-tracking tracers downwards (maintianing that the number of timesteps is equal to the number of nodes that the tracer visits)
                t=t-1 #de-incriment time
                path[l,w,h]=path[l,w,h]-1 #delete the backtrack
           
                #find where the tracer came from
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
                
                #put the tracer one level lower
                if(h!=0):
                    h=h-1
                    path[l,w,h]+=1
             
                else:
                    h=h
                    earth =  np.concatenate((np.random.rand(L,W,1),earth),2)
                    path =  np.concatenate((np.zeros((L,W,1)),path),2)
                    aggregatePaths =  np.concatenate((np.zeros((L,W,1)),aggregatePaths),2)
                    append+=1
                    path[l,w,h]+=1
                  
    return [l,w,h,t]

#condensed code for creating the 3d travel plot
def Threedplot():
    
    X=range(0,L)
    Y=range(0,W)
    Z=range (0,H+append)
    
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8], projection='3d')
    ax_cb = fig.add_axes([0.8, 0.3, 0.05, 0.45])
    ax.set_aspect('equal')
    plotMatrix(ax, X, Y, Z, aggregatePaths, cmap="jet", cax = ax_cb)
    ax.set_xlabel('N-S axis, north towards zero')
    ax.set_ylabel('E-W axis, east towards zero')
    ax.set_zlabel('h axis, down towards zero')
    ax.set_title('number of visits at each node')


    #condensed code for the correlation plot
def correlation()  :  
    a,b,c = earth.shape
    aggPlt= np.reshape(aggregatePaths, a*b*c)
    earthPlt= np.reshape(earth, a*b*c)
    

    plt.figure()
    plt.plot(earthPlt, aggPlt, 'ro')
    plt.xlabel('probability')
    plt.ylabel('number of visits')
    
    




timesteps =15 #number of moves each tracer makes

# run a trace from everywhere on the land surface to simulate a rainfall event. In futur eversions, multiple rainfall events could be simulated, and the tracers would move
# simmultaneously.
#no erosion
for i in range(0,L-1):
    for j in range(0,W-1):
        h=H+append-1
        l=i
        w=j
        path[l,w,h]=1
        t=0
        #step through time for the moving object
        while  t < timesteps:
           #central points
           
        
            l,w,h,t =      advance_step(l,w,h,t) 
            t+=1    
        aggregatePaths =  aggregatePaths + path #add the most recent path to the complete colection
        path = np.zeros((L,W,H+append))    


   
#create 3d path visualization
Threedplot()
plt.show

#plot correlation
correlation()
plt.title('without erosion')
plt.show()


#now lets do the same experiment, with erosion

for i in range(0,L-1):
    for j in range(0,W-1):
        h=H+append-1
        l=i
        w=j
        path[l,w,h]=1
        t=0
        #step through time for the moving object
        while  t < timesteps:
           #central points
           
        
            l,w,h,t =      advance_step(l,w,h,t) 
            t+=1    
        aggregatePaths =  aggregatePaths + path
        #probability enhancment: enhance each traveled spot by 10% (erosion)
        a,b,c= path.shape   
        pathR= np.reshape(path, a*b*c)
        earthR= np.reshape(earth, a*b*c)
        for k in range(0, a*b*c):
            if pathR[k] == 1:
                earthR[k]= earthR[k]*1.1
        earth = np.reshape(earthR, (a,b,c))
        path = np.zeros((L,W,H+append))     

#plot the erosion results
Threedplot()
plt.show
correlation()
plt.title('with erosion')
plt.show()


#repeate experinent with vertical fracture

earth[:,2 ,:]=1 #fracture the land
for i in range(0,L-1):
    for j in range(0,W-1):
        h=H+append-1
        l=i
        w=j
  
        path[l,w,h]=1
        t=0
            #step through time for the moving object
        while  t < timesteps:
           #central points
           
        
            l,w,h,t =      advance_step(l,w,h,t) 
            t+=1    
            aggregatePaths =  aggregatePaths + path
            #probability enhancment: enhance each traveled spot by 10%
            a,b,c= path.shape   
            pathR= np.reshape(path, a*b*c)
            earthR= np.reshape(earth, a*b*c)
            earth[:,2 ,:]=1 #include fracture in appended layers
        for k in range(0, a*b*c):
            if pathR[k] == 1 and earthR[k]!=1:
                earthR[k]= earthR[k]*1.1
        earth = np.reshape(earthR, (a,b,c))
        path = np.zeros((L,W,H+append))   

#plot fracture results       
Threedplot()
plt.show
correlation()
plt.title('with fracture')
plt.show()



#Future versions: turn nested for loops into a method that can be called with various paramenters such as fractures, impermeable layers
#Future versions: allow the 'earth' to expand in the L and W direction  in addition to the H direction
