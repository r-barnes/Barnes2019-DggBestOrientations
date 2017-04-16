#!/usr/bin/env python

import numpy as np
import rtree
import fiona
import shapely.geometry
import shapely.speedups
import shapely.ops
import pyproj
import math
from matplotlib import pyplot as plt

def TransformLatLon(latr,lonr,plat,plon):
  """Take a point at (latr,lonr) in a system with a pole at (90,*) and rotate it
     into a system with a pole at (plat,plon) while first rotating it ptheta"""
  latr = latr.copy()
  lonr = lonr.copy()
  latr *= np.pi/180
  lonr *= np.pi/180
  plat = plat*np.pi/180
  plon = plon*np.pi/180
  xr   = np.cos(lonr)*np.cos(latr)
  yr   = np.sin(lonr)*np.cos(latr)
  zr   = np.sin(latr)
  x    = np.cos(plat)*np.cos(plon)*xr + np.sin(plon)*yr + np.sin(plat)*np.cos(plon)*zr
  y    = -np.cos(plat)*np.sin(plon)*xr + np.cos(plon)*yr - np.sin(plat)*np.sin(plon)*zr
  z    = -np.sin(plat)*xr + np.cos(plat)*zr
  lat  = np.arcsin(z)     *180/np.pi
  lon  = np.arctan2(y, x) *180/np.pi
  return lat,lon

def wgs_to_mercator(lat, lon):
  prj_wgs = pyproj.Proj(init='epsg:4326')
  prj_mer = pyproj.Proj(init='epsg:3857')
  lat     = np.clip(lat,-85,85)
  #lon     = np.clip(lon,-179.9,180)
  x, y = pyproj.transform(prj_wgs, prj_mer, lon-0.0001, lat-0.001, radians=False) #Yes, `lon,lat` is the correct order here
  return x, y

def IntersectsPoly(x,y,poly):
  return poly.contains(shapely.geometry.Point(x,y))

def PlotPoints(px,py):
  fig   = plt.figure()
  ax    = fig.add_subplot(111)
  for i in features[0:8]:
    pol_x,pol_y = i.exterior.xy
    ax.plot(pol_x,pol_y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round')
  ax.scatter(px,py,color='green',marker='o')
  plt.show()

def PlotFeature(i):
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  pol_x,pol_y = features[i].exterior.xy
  ax.plot(pol_x,pol_y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round')
  plt.show()

def PlotPolyPoint(x,y):
  fig   = plt.figure()
  ax    = fig.add_subplot(111)
  isecs = ridx.intersection((x,y,x,y), objects="raw")
  isecs = [p for p in isecs if IntersectsPoly(x,y,p)]
  for i in isecs:
    px,py = i.exterior.xy
    ax.plot(px, py, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round')
  ax.scatter(x,y,color='green',marker='o')
  plt.show()

shapely.speedups.enable()


#https://en.wikipedia.org/wiki/Regular_icosahedron#Spherical_coordinates
#The locations of the vertices of a regular icosahedron can be described using
#spherical coordinates, for instance as latitude and longitude. If two vertices
#are taken to be at the north and south poles (latitude ±90°), then the other
#ten vertices are at latitude ±arctan( 1 / 2 ) ≈ ±26.57°. These ten vertices are
#at evenly spaced longitudes (36° apart), alternating between north and south
#latitudes. This scheme takes advantage of the fact that the regular icosahedron
#is a pentagonal gyroelongated bipyramid, with D5d dihedral symmetry—that is, it
#is formed of two congruent pentagonal pyramids joined by a pentagonal
#antiprism.

#Generate lat-lon pairs of where the vertices are
vertices = [(90,0),(-90,0)]
for theta in enumerate(np.arange(-180,180,36)):
  if theta[0]%2==0: #Arbitrarily call even points the northernly ones
    vertices.append( (26.57,theta[1]) )
  else:
    vertices.append( (-26.57,theta[1]) )

vertices = np.array(vertices)
olats    = vertices[:,0]
olons    = vertices[:,1]


#http://openstreetmapdata.com/data/land-polygons
shapefilename     = '/home/rick/projects/dgg_best_poles/simplified-land-polygons-complete-3857/simplified_land_polygons.shp'
features          = [x for x in fiona.open(shapefilename)]                       #Read shapefile
features          = [shapely.geometry.shape(x['geometry']) for x in features]    #Convert to shapely shapes
features          = [x.simplify(40000) for x in features]                        #Simplify shapes for speed
features.sort(key = lambda poly: len(poly.exterior.xy[0]), reverse=True)

ridx      = rtree.index.Index([ (i,x.bounds,x) for i,x in enumerate(features) if x.exterior is not None ]) #Build spatial index
ufeatures = shapely.ops.cascaded_union(features)

#FeatureID (when sorted) | Landform
#0                       | Eurasia
#1                       | North & South America
#2                       | Greenland
#3                       | Antarctica
#4                       | Baffin Island
#5                       | Ellesmere
#6                       | Australia
#7                       | Africa
#8                       | England

found = []
for lat in np.arange(10,52,0.5):
  print(lat)
  for lon in np.arange(0,72,0.5):
    x,y = wgs_to_mercator(*TransformLatLon(olats,olons,lat,lon))
    pts = zip(x,y)
    good = len(olats)
    for p in pts:
      isecs = ridx.intersection((p[0],p[1],p[0],p[1]), objects="raw")
      isecs = [poly for poly in isecs if IntersectsPoly(p[0],p[1],poly)]
      if len(isecs)==0:
        good -= 1
    found.append( (good,lat,lon) )

found.sort(key=lambda x: x[0])



fout = open('/z/out','w')
for x in found:
  fout.write("{0} {1} {2}\n".format(*x))
fout.close()

i       = 419
wpts    = wgs_to_mercator(*TransformLatLon(olats,olons,found[i][1],found[i][2]))
uvertex = shapely.geometry.MultiPoint([shapely.geometry.Point(*x) for x in list(zip(*wpts))])
ufeatures.distance(uvertex)
PlotPoints(*wpts)
