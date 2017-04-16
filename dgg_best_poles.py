#!/usr/bin/env python

import numpy as np
import numpy.linalg as LA
import rtree
import fiona
import shapely.geometry
import shapely.speedups
import shapely.ops
import pyproj
import math
from matplotlib import pyplot as plt

def CountPoints(geom):
  """Count the number of points used to define a geometry"""
  if geom.type == 'Polygon':
    exterior_coords = len(geom.exterior.xy[0])
    interior_coords = 0
  elif geom.type == 'MultiPolygon':
    exterior_coords = 0
    interior_coords = 0
    for part in geom:
      epc = CountPoints(part)  # Recursive call
      exterior_coords += epc['ext']
      interior_coords += epc['int']
  else:
    raise ValueError('Unhandled geometry type: ' + repr(geom.type))
  return {'ext': exterior_coords,
          'int': interior_coords}

def PairwiseCircular(iterable): #TODO: Cut
  """
  Iterate through an itertable returning adjacent pairs, including first-last
  :param iterable   An iterable
  :returns: Pairs of sequential, adjacent entries from the iterable, as well as
            the first-last pair
  """
  it    = iter(iterable)
  a     = next(it, None)
  first = a
  for b in it:
    yield (a, b)
    a = b
  yield (a,first)

def Pairwise(iterable):
  """
  Iterate through an itertable returning adjacent pairs
  :param iterable   An iterable
  :returns: Pairs of sequential, adjacent entries from the iterable
  """
  it    = iter(iterable)
  a     = next(it, None)
  for b in it:
    yield (a, b)
    a = b

def LatLonToXYZ(lat,lon,radius):
  """
  Convert a latitude-longitude pair to 3D-Cartesian coordinates
  :param lat    Latitude in degrees
  :param lon    Longitude in degrees
  :param radius Radius in arbitrary units
  :returns: x,y,z coordinates in arbitrary units
  """
  lat = np.radians(lat)
  lon = np.radians(lon)
  x   = radius * np.cos(lon) * np.cos(lat)
  y   = radius * np.sin(lon) * np.cos(lat)
  z   = radius * np.sin(lat)
  return x,y,z

def XYZtoLatLon(x,y,z):
  """
  Convert 3D-Cartesian coordinates to a latitude-longitude pair
  :param x      x-coordinate in arbitrary units
  :param y      y-coordinate in arbitrary units
  :param z      z-coordinate in arbitrary units
  :returns A (lat,lon) pair in degrees
  """
  radius = np.sqrt(x*x+y*y+z*z)
  lat    = np.degrees(np.arcsin(z/radius))
  lon    = np.degrees(np.arctan2(y, x))   
  return lat,lon

def Haversine(lat1, lon1, lat2, lon2):
  """
  Calculate the Great Circle distance on Earth between two latitude-longitude
  points
  :param lat1 Latitude of Point 1 in degrees
  :param lon1 Longtiude of Point 1 in degrees
  :param lat2 Latitude of Point 2 in degrees
  :param lon2 Longtiude of Point 2 in degrees
  :returns Distance between the two points in kilometres
  """
  Rearth = 6371
  lat1   = np.radians(lat1)
  lon1   = np.radians(lon1)
  lat2   = np.radians(lat2)
  lon2   = np.radians(lon2)
  #Haversine formula 
  dlon = lon2 - lon1 
  dlat = lat2 - lat1 
  a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
  c = 2 * np.arcsin(np.sqrt(a)) 
  return Rearth*c

#http://stackoverflow.com/a/1302268/752843
def NearestPointOnGC(alat1,alon1,alat2,alon2,plat,plon):
  """
  Calculate the location of the nearest point on a Great Circle to a query point
  :param lat1 Latitude of start of arc in degrees
  :param lon1 Longtiude of start of arc in degrees
  :param lat2 Latitude of end of arc in degrees
  :param lon2 Longtiude of end of arc in degrees
  :param plat Latitude of query point in degrees
  :param plon Longitude of query point in degrees
  :returns: A (lat,lon) pair in degrees of the closest point
  """
  Rearth    = 6371 #km
  #Convert everything to Cartesian coordinates
  a1        = np.array(LatLonToXYZ(alat1,alon1,Rearth))
  a2        = np.array(LatLonToXYZ(alat2,alon2,Rearth))
  p         = np.array(LatLonToXYZ(plat, plon, Rearth))
  G         = np.cross(a1,a2) #Plane of the Great Circle containing A and B
  F         = np.cross(p,G)   #Plane perpendicular to G that passes through query pt
  T         = np.cross(G,F)   #Vector marking the intersection of these planes
  T         = Rearth*T/LA.norm(T) #Normalize to lie on the Great Circle
  tlat,tlon = XYZtoLatLon(*T)
  return tlat,tlon

def DistanceToGCArc(alat,alon,blat,blon,plat,plon):
  """
  Calculate the distance from a query point to the nearest point on a
  Great Circle Arc
  :param lat1 Latitude of start of arc in degrees
  :param lon1 Longtiude of start of arc in degrees
  :param lat2 Latitude of end of arc in degrees
  :param lon2 Longtiude of end of arc in degrees
  :param plat Latitude of query point in degrees
  :param plon Longitude of query point in degrees
  :returns: The distance in kilometres from the query point to the great circle
            arc
  """
  tlat,tlon = NearestPointOnGC(alat,alon,blat,blon,plat,plon) #Nearest pt on GC
  abdist    = Haversine(alat,alon,blat,blon)  #Length of arc
  atdist    = Haversine(alat,alon,tlat,tlon)  #Distance arc start to nearest pt
  tbdist    = Haversine(tlat,tlon,blat,blon)  #Distance arc end to nearest pt
  #If the nearest point T on the Great Circle lies within the arc, then the
  #length of the arc is approximately equal to the distance from T to each end
  #of the arc, accounting for floating-point errors
  PRECISION = 1e-3 #km 
  #We set the precision to a relatively high value because super-accuracy is not
  #to needed here and a failure to catch this can lead to vast under-estimates
  #of distance
  if np.abs(abdist-atdist-tbdist)<PRECISION: #Nearest point was on the arc
    return Haversine(tlat,tlon,plat,plon)
  #Okay, the nearest point wasn't on the arc, so the nearest point is one of the
  #ends points of the arc
  apdist = Haversine(alat,alon,plat,plon)
  bpdist = Haversine(blat,blon,plat,plon)
  return min(apdist,bpdist)

def Distance3dPointTo3dPolygon(lat,lon,geom):
  """
  Calculate the closest distance between a polygon and a latitude-longitude
  point, using only spherical considerations
  :param lat  Latitude of query point in degrees
  :param lon  Longitude of query point in degrees
  :param geom A `shapely` geometry whose points are in latitude-longitude space
  :returns: The minimum distance in kilometres between the polygon and the
            query point
  """
  if geom.type == 'Polygon':
    dist = math.inf
    xy   = features[0]['geometry'][0].exterior.xy
    #Polygons are closed rings, so the first-last pair is automagically delt with
    for p1, p2 in Pairwise(zip(*xy)):
      dist = min(dist,DistanceToGCArc(p1[1],p1[0],p2[1],p2[0],lat,lon))
  elif geom.type == 'MultiPolygon':
    dist = min(*[Distance3dPointTo3dPolygon(lat,lon,part) for part in geom])
  return dist

def TransformLatLon(latr,lonr,plat,plon):
  """Take a point at (latr,lonr) in a system with a pole at (90,*) and rotate it
     into a system with a pole at (plat,plon) while first rotating it ptheta"""
  latr       = latr.copy()
  lonr       = lonr.copy()
  latr       = np.radians(latr)
  lonr       = np.radians(lonr)
  plat       = np.radians(plat)
  plon       = np.radians(plon)
  xr, yr, zr = LatLonToXYZ(latr,lonr,1)
  x          =  np.cos(plat)*np.cos(plon)*xr + np.sin(plon)*yr + np.sin(plat)*np.cos(plon)*zr
  y          = -np.cos(plat)*np.sin(plon)*xr + np.cos(plon)*yr - np.sin(plat)*np.sin(plon)*zr
  z          = -np.sin(plat)*xr + np.cos(plat)*zr
  lat, lon   = XYZtoLatLon(x,y,z)
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
