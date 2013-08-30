from visual import *
import numpy as np
from scipy.spatial import Voronoi
from random import random

norm2tolerance = 1E-10 # tolerance in a mag2.
angletolerance = 1E-9 # tolerance before two angles are seen as identical
infinity = 5 # how far out infinity goes

colors = [ (1,0,0), (0,0,1),  (0,1,0),(1,0,1),(0,1,1),(1,1,0),
           (0.5,1,0.5), (1,0.5,0.5),(0.5,0.5,1), (1,1,1),
           (0.7,1,0), (0,1,0.7), (1,0.7,0),(1,0,0.7), (0.7,0,1),(0,0.7,1) ]

def error(string):
    print 'ERROR!'
    print string
    exit(1) #bad

def unit( v1 , errormsg="vector input to unit() is too small in magnitude"):
    v1 = vector(v1)
    v1norm2 = mag2( v1 )
    if ( v1norm2 < norm2tolerance ):
        error(errormsg)
    return v1 / sqrt(v1norm2) 


def randomUnitVector():
    ##gives random unit vectors uniformly distributed on a sphere
    z1 = 2*random() - 1
    rho1= sqrt( 1 - z1**2)
    phi1 = 2 * pi * random()
    x1 = rho1 * cos( phi1 )
    y1 = rho1 * sin( phi1 )
    return vector( x1, y1, z1 )

def randomPerpVector( v1 ):
    ##gives a random unit vector perpendicular to the given vector
    u1 = unit(v1, "input to randomPerpVector appears to be zero vector")
    u2 = randomUnitVector() 
    count1 = 1
    count1max = 10
    while ( 1 - fabs(dot( u1, u2)) < 10*norm2tolerance  and count1 < count1max ):
        u2 = randomUnitVector()
        count1 += 1
    
    if ( count1 == count1max ):
        error(" too many iterations, couldn't find perpendicular vector ")

    return unit(cross(u1, u2))



def Face(v1, v2, v3, color1=(0.5,0.5,0.5) ):
    return faces( pos=np.array([ v1, v2, v3 ]), color=color1)

def squareFace(origin, v1, v2, color1=(0.5,0.5,0.5)):
    return faces( pos=np.array([ origin+0.5*(v1+v2), origin+0.5*(-v1+v2), origin+0.5*(-v1-v2),
                        origin+0.5*(-v1-v2), origin+0.5*(v1-v2), origin+0.5*(v1+v2)  ]),
                  color=np.array( [color1]*6 )  )

def squareDoubleFace(origin, v1, v2, color1=(0.5,0.5,0.5)):
    return faces( pos=np.array([ origin+0.5*(v1+v2), origin+0.5*(-v1+v2), origin+0.5*(-v1-v2),
                        origin+0.5*(-v1-v2), origin+0.5*(v1-v2), origin+0.5*(v1+v2),
                        origin+0.5*(v1+v2), origin+0.5*(-v1-v2), origin+0.5*(-v1+v2), 
                        origin+0.5*(-v1-v2), origin+0.5*(v1+v2), origin+0.5*(v1-v2) ]), 
                  color=np.array( [color1]*12 )  )


def circumcircle( pointlist , normdir=-1 ):
    N = len(pointlist)
    if N <> 3:
        error("need 3 points for circumcircle")
    else:
        p1 = vector(pointlist[0])
        p2 = vector(pointlist[1])
        p3 = vector(pointlist[2])
        if normdir == -1:
            normdir = cross( p2 - p1, p3 - p1 )
            normmag2 = mag2(normdir)
            if ( normmag2 < norm2tolerance):
                error(" 3 points colinear, no circumcircle (except at infinity?) ")
            else:
                normdir /= sqrt(normmag2)
        center = 1.0*(p1+p2+p3)/3
        R1 = vector(p1 - center)
        R2 = vector(p2 - center)
        R3 = vector(p3 - center)

######## in the plane of the points, create some new x,y directions:
        ex = unit(vector( p2 - p1 ))##arbitrary "x" direction in the plane
        ey = cross( normdir, ex ) #arbitrary "y" direction
######## create some strange vectors

        Rx = vector( dot(R1, ex), dot(R2, ex), dot(R3, ex) )
        Ry = vector( dot(R1, ey), dot(R2, ey), dot(R3, ey) )
        M2 = vector( mag2(R1), mag2(R2), mag2(R3) )
        Rxy = sum( cross( Rx, Ry ) )
        ux = -0.5*sum(  cross( Ry, M2 ) ) / Rxy
        uy =  0.5*sum(  cross( Rx, M2 ) ) / Rxy
        
        circumcenterc = ux * ex + uy * ey ##circumcenter with respect to center of system
        raysout = np.array([ unit(cross(normdir, pointlist[0]-pointlist[1])),    
                             unit(cross(normdir, pointlist[1]-pointlist[2])),
                             unit(cross(normdir, pointlist[2]-pointlist[0]))  ])

#        angles = [ acos( dot( unit(R1 - circumcenterc), unit(R2 - circumcenterc)) ), 
#                            acos( dot(unit(R2 - circumcenterc), unit(R3 - circumcenterc)) ),
#                            acos( dot(unit(R3 - circumcenterc), unit(R1 - circumcenterc)) ) ] 
#        biggestangleindex = 0
#        if angles[1] > angles[0]:
#            biggestangleindex = 1
#        if angles[2] > angles[biggestangleindex]:
#            biggestangleindex = 2
#        biggestangle = angles.pop(biggestangleindex)
#        if fabs(biggestangle  - sum( angles ) ) < angletolerance:
############# flip Rav2cc corresponding to angle
#            print "flipping the bird"
#            Rav2cc[biggestangleindex] *= -1
#
#
#        

        return [ circumcenterc + center, raysout   ]

## PLOTTING
# We can plot it as follows. First the points and the Voronoi vertices:

## plt.plot(points[:,0], points[:,1], 'o')
## plt.plot(vor.vertices[:,0], vor.vertices[:,1], '*')
## plt.xlim(-1, 3); plt.ylim(-1, 3)

#Plotting the finite line segments goes as for the convex hull, but now we have to guard for the infinite edges:

## # # for simplex in vor.ridge_vertices:
##     simplex = np.asarray(simplex)
##     if np.all(simplex >= 0):
##         plt.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], 'k-')

### # # The ridges extending to infinity require a bit more care:
## center = points.mean(axis=0)
## for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
##     simplex = np.asarray(simplex)
##     if np.any(simplex < 0):
##         i = simplex[simplex >= 0][0] # finite end Voronoi vertex
##         t = points[pointidx[1]] - points[pointidx[0]] # tangent
##         t /= np.linalg.norm(t)
##         n = np.array([-t[1], t[0]]) # normal
##         midpoint = points[pointidx].mean(axis=0)
##         far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
##         plt.plot([vor.vertices[i,0], far_point[0]],
##                  [vor.vertices[i,1], far_point[1]], 'k--')
## plt.show()




class infiniteArrow:
    def __init__( self, origin, ray ):
        self.origin = vector(origin)
        self.ray = vector(ray)

class infiniteAnglePlane:
    def __init__( self, origin, ray1, ray2 ):
        self.origin = vector(origin) 
        self.ray1 = unit(ray1, "ray1 too small in infiniteAnglePlane") 
        self.ray2 = unit(ray2, "ray2 too small in infiniteAnglePlane")

        ##if an error here, the rays were not different enough
        self.normal = unit(cross(self.ray1, self.ray2), 
                "ray1 and ray2 insufficiently differet in infiniteAnglePlane")  

class infinitePyramid:
    def __init__( self, origin, rays ):
        self.origin = vector(origin)
        self.rays = np.array( rays )
        self.N = len( self.rays )
        

class infinitePlane:
    def __init__( self, origin, normal, color=(0.5,0.5,0.5), color2 = -1 ):
        self.color = color
        if color2 == -1:
            self.color2 = color
        else:
            self.color2 = color2
        self.origin = vector(origin)
        self.normal = unit(normal, "input normal to infintePlane appears to be zero")
        self.u1 = randomPerpVector( self.normal )
        self.u2 = cross( self.normal, self.u1 )
    def plot( self ):
        self.f = squareFace( self.origin, infinity*self.u1, infinity*self.u2, self.color )
        self.f2 = squareFace( self.origin,infinity*self.u2, infinity*self.u1,  self.color2 )





class infiniteHalfPlane:
    global norm2tolerance
    def __init__( self, origin, indir, alongdir, color=(0.5,0.5,0.5), color2=-1):
        self.color = color
        if color2 == -1:
            self.color2 = color
        else:
            self.color2 = color2
        self.origin = vector(origin)
        self.alongdir = unit(alongdir, "alongdir insufficient in infiniteHalfPlane")
        self.indir = vector(indir)
        self.indir -= dot(self.alongdir, self.indir)*self.alongdir
        self.indir = unit(self.indir, "indir insufficient in magnitude after projection")
    def plot( self ):
        self.f = squareFace( self.origin+0.25*infinity*self.indir, 
                                    0.5*infinity*self.indir, 
                                    infinity*self.alongdir, self.color )
        self.f2 = squareFace( self.origin+0.25*infinity*self.indir, 
                                    infinity*self.alongdir, 
                                    0.5*infinity*self.indir,  self.color2 )


class Loronoi:
    def __init__( self, points ):
######## reassign points to a numpy array
        self.points = np.array( points )

######## check and see if all points are unique.  if not, raise an error.
        for i in range(len(self.points)):
            for j in range(i+1, len(self.points)):
                dist2 = mag2( self.points[i] - self.points[j] )
                if dist2 < norm2tolerance:
                    error("points "+str(i)+", "+str(j)+" are too close.")


######## grab some standard properties of the set:
        self.center = self.points.mean(axis=0)
        self.N = len(points)
        self.planes =  [] 
        self.circumcenters =  [] 
        self.circumcirclepointindices =  [] 
        self.raysout = []

######## init ()
        if (self.N == 1):
            pass 
############
############ end N == 1
############
        elif (self.N == 2):
            normal = vector(self.points[1] - self.points[0])
            self.planes += [  infinitePlane( self.center, normal, colors[1], colors[0] )  ]
############
############ end N == 2
############
        elif (self.N == 3):
############
############ CHECK FOR POINTS BEING COLINEAR
############
            alongdir = cross( self.points[1] - self.points[0], 
                              self.points[2] - self.points[0] )
            alongmag2 = mag2(alongdir)
            if ( alongmag2 < norm2tolerance):
################
################ POINTS ARE COLINEAR
################
                self.planes += [ infinitePlane( 0.5*(self.points[0]+self.points[1]), 
                                                self.points[1] - self.points[0],
                                                colors[1], colors[0] ),
                                 infinitePlane( 0.5*(self.points[2]+self.points[1]), 
                                                self.points[2] - self.points[1],
                                                colors[2], colors[1] ) ]
################
            else:
################
################ IF POINTS ARE NOT COLINEAR THEN YOU CAN USE CIRCUMCIRCLE TYPE STUFF
################
                alongdir /= sqrt( alongmag2 ) 
                circumcenter, raysout = circumcircle( self.points, alongdir )
                self.circumcenters = [ circumcenter ]
                self.raysout = [ raysout ]

                self.planes += [ infiniteHalfPlane( self.circumcenters[0], 
                                    raysout[0], #cross(alongdir, self.points[0]-self.points[1]),
                                    alongdir, colors[0], colors[1]   ),
                                 infiniteHalfPlane( self.circumcenters[0], 
                                    raysout[1], #cross(alongdir, self.points[1]-self.points[2]),
                                    alongdir, colors[1], colors[2]  ) ,
                                 infiniteHalfPlane( self.circumcenters[0], 
                                    raysout[2], #cross(alongdir, self.points[2]-self.points[0]),
                                    alongdir, colors[2] , colors[0]  ) ]
################
############
############ end N == 3
############
        elif self.N >= 4:
            self.vor = Voronoi(points)
            print "vertices"
            print self.vor.vertices
            print "regions"
            print self.vor.regions

            print "ridge_points"
            print self.vor.ridge_points
            print "ridge_vertices"
            print self.vor.ridge_vertices


            if self.N == 4:
                for i in range(4):
                    circle = range(4)
                    circle.remove(i)
                    circumcenter, raysout = circumcircle( self.points[ circle ] )
                    self.circumcirclepointindices.append( circle )
                    self.circumcenters.append( circumcenter)
                    self.raysout.append( raysout )
#            circle = range(4)
#            circle.remove(0)
#            self.circumcenters.append( circumcircle( [ self.points[c] for c in circle ] ) )
#
#            for i in range(4):
#                for j in range(i+1,4):
#                    self.circumcenters.append( circumcircle([ self.points[i], self.points[j] ] ))
#                for c in circle:
#                    print "i = "+str(i)+": circum radius = "+str(mag(self.circumcenters[i] - self.points[c]))
                

################
############
############ end N >= 4
############
######## 
######## END init()
######## 
####
#### class Loronoi
####
    def which_region( point ):
######## 
######## determines which region the point belongs to.
######## 
        if (self.N == 1):  #if only one original point, any other point belongs to region 0
            return [0]
        elif (self.N == 2):  #if two original points...
            norm0 = vector( point - self.points[0] )
            norm0 = mag2(norm0)
            norm1 = vector( point - self.points[1] )
            norm1 = mag2(norm1)
            diff = norm1 - norm0
            if ( fabs(diff) < norm2tolerance ):
                return [0,1]  #split the difference between the two
            elif (diff > 0):
                return [1]
            else:
                return [0]
################
############
########
        elif (self.N == 3):  #if three original points
            pass 
############
######## 
######## END which_region()
######## 
####
#### class Loronoi
####
    def plot(self):
        for i in range(len(self.points)):
            sphere(pos = self.points[i], color=colors[i], radius=0.1)

        for c in self.circumcenters:
            cone(pos = c, radius=0.1,length=0.1, height=0.1)
            

        for pl in self.planes:
            pl.plot()

        if ( self.N <= 3 ):
            pass
        elif ( self.N >= 4 ):
            for v in self.vor.vertices:
                box(pos = v, color=(0.5,0.5,0.5), width=0.1,height=0.1,length=0.1)

            for i in range(len(self.raysout)):
                c = self.circumcenters[i]
                for j in range(3):
                    cylinder(pos = c, radius=0.01, axis=0.5*infinity*self.raysout[i][j],
                            color = tuple(0.5*(vector(colors[self.circumcirclepointindices[i][j]])
                                               +vector(colors[self.circumcirclepointindices[i][mod(j+1,3)]])) )
                            )



#            for pointidx, simplex in zip(self.vor.ridge_points, self.vor.ridge_vertices):
#                if np.all( simplex >= 0 ): ##make sure all indices in r are >=0 
#                    print simplex
#                    Face( self.vor.vertices[simplex[0]], 
#                            self.vor.vertices[simplex[1]], 
#                            self.vor.vertices[simplex[2]], 
#                            colors[pointidx[0]] )
#                    Face( self.vor.vertices[simplex[0]], 
#                            self.vor.vertices[simplex[2]], 
#                            self.vor.vertices[simplex[1]], 
#                            colors[pointidx[1]] )
#                for c in self.circumcenters:
#                    axis1 = vector(c-v)
#                    axis1 /= mag(axis1)
#                    cylinder(pos = v, radius=0.01, axis=0.5*infinity*axis1)


############
########
######## END plot()
########
####

#### END class Loronoi
####
