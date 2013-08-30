from Loronoi import *


## ONE POINT
#points = np.array([ [0,0,0] ])

## TWO POINTS
#points = np.array([  [0,0,0],  [0,1,0]  ]) 

## THREE POINTS
#points = np.array([  [0,2,0], [0,1,0],  [0,-1,0]  ])    #colinear
#points = np.array([  [0,2,0], [5,-3,0.1],  [0,1,0]  ])  #circumcenter external
#points = np.array([  [0,2,0], [-1,0,-1],  [0,0,1]  ])  #circumcenter internal

## FOUR POINTS
#points = np.array([  [0,-1,0], [1,0,1],  [0,1,0.3], [1,1,0]  ])  #random
points = np.array([  [0.3,-1,0.1], [1,0,1],  [0,1,1], [1,1,0]  ])  #internal circumcenters
#points = np.array([  [0,-1,0], [1,0,0],  [0,1.0/sqrt(2),1.0/sqrt(2)], [1.0/sqrt(2),0,1.0/sqrt(2)]  ])  #random

## RANDOM POINTS
#N=10
#randoms = []
#for i in range(N):
#    randoms.append( (1+5*random())*randomUnitVector()) 
# 
#points = np.array(randoms)



lor = Loronoi(points)

lor.plot()

while True:
    rate(60)
    #scene.center.y += 0.01
