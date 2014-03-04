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

minideltakey = 0.05
megadeltakey = 0.15
#prose = label()
scene.right = cross( scene.forward, scene.up )
scene.autoscale = False

centerball = sphere(pos=scene.center, color=(1,1,1), radius=0.01)
centerhelix = helix(pos=scene.center, color=(0.5,0.5,0.5), radius=0.05, thickness=0.01, axis=scene.forward*0.001, coils=1)
scene.lights = [ distant_light( direction=-scene.up, color=color.gray(0.8) ),
                         distant_light( direction=-scene.forward, color=color.gray(0.5) ) ]
while True:
    rate(60)
    if scene.kb.keys: # event waiting to be processed?
        string = scene.kb.getkey() # get keyboard info
        alts = string.split("+")
        if len(alts) > 1:
            alt = alts[0]
            deltakey = megadeltakey
            s = alts[1]
        else:
            s = alts[0]
            if s.upper() == s:
                alt = "shift"
                s = s.lower() 
                deltakey = megadeltakey
            else:
                alt = ""
                deltakey = minideltakey

        if (s == "left"):
            deltaforward = -scene.right
            scene.forward += deltakey*deltaforward
            scene.forward /= mag(scene.forward)

            scene.right = cross( scene.forward, scene.up )

        elif (s == "right"):
            deltaforward = scene.right
            scene.forward += deltakey*deltaforward
            scene.forward /= mag(scene.forward)
            
            scene.right = cross( scene.forward, scene.up )

        elif (s == "up"):
            deltaforward = scene.up

            scene.forward += deltakey*deltaforward
            scene.forward /= mag(scene.forward)
            
            scene.up = cross( scene.right, scene.forward )

        elif (s == "down"):
            deltaforward = -scene.up

            scene.forward += deltakey*deltaforward
            scene.forward /= mag(scene.forward)
            
            scene.up = cross( scene.right, scene.forward )

        elif (s == "w"):
            scene.center += deltakey*scene.up

        elif (s == "a"):
            scene.center -= deltakey*scene.right
        
        elif (s == "s"):
            scene.center -= deltakey*scene.up

        elif (s == "d"):
            scene.center += deltakey*scene.right
        
        elif (s == "q"):
            scene.center -= deltakey*scene.forward
        
        elif (s == "e"):
            scene.center += deltakey*scene.forward

#       # these next guys don't work in my vpython for some reason.
#        elif (s == "z"):
#            scene.scale *= 1.1*deltakey
#        
#        elif (s == "c"):
#            scene.scale *= 0.9*deltakey

        scene.lights = [ distant_light( direction=-scene.up, color=color.gray(0.8) ),
                         distant_light( direction=-scene.forward, color=color.gray(0.5) )
                       ]
#        elif ((s == 'backspace' or s == 'delete') and
#            len(prose.text)) > 0: 
#            prose.text = prose.text[:-1] # erase letter
#        elif s == 'shift+delete': 
#            prose.text = '' # erase all text
#        else: #if len(s) == 1: 
#            prose.text += s # append new character

        centerball.pos = scene.center
        centerhelix.pos = scene.center
        centerhelix.axis = scene.forward*0.001


