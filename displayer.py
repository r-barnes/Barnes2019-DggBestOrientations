#!/usr/bin/env python3
import pyproj
import numpy as np
from matplotlib import pyplot as plt
import os
import pickle
import sys

unprojected_landmassfile = 'data/land-polygons-complete-4326/land_polygons.shp'
plottable_landmassfile   = 'data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp'
landmassfile             = 'data/land-polygons-split-3857/land_polygons.shp'

#storedir = '/home/rbarnes1/scratch/dgg_best'
storedir  = 'temp'

def PlotPoints(px,py):
  """Plot points along with some of the largest land masses, for content"""
  fig   = plt.figure()
  ax    = fig.add_subplot(111)
  for f in plottable_landmasses[0:8]:
    pol_x,pol_y = f.exterior.xy
    ax.plot(pol_x,pol_y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round')
  ax.scatter(px,py,color='green',marker='o')
  plt.show()



def TransformLatLon(latr,lonr,plat,plon,ptheta):
  """Take a point at (latr,lonr) in a system with a pole at (90,*) and rotate it
     into a system with a pole at (plat,plon) while first rotating it ptheta"""
  latr = latr.copy()
  lonr = lonr.copy()
  plat = np.radians(plat)
  plon = np.radians(plon)
  lonr += 180                               #Move lonr to the [0,360] system
  lonr += ptheta                            #Rotate points by longitude
  lonr = np.fmod(360+np.fmod(lonr,360),360) #Map everything to [0,360]
  lonr -= 180                               #Move back to the [-180,180] system
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

def GeneratePlottableLandmasses():
  if not os.path.isfile(os.path.join(storedir,'plottable_landmasses.pickle')):
    plottable_landmasses = [x for x in fiona.open(plottable_landmassfile)]
    plottable_landmasses = [sg.shape(x['geometry']) for x in plottable_landmasses]
    plottable_landmasses.sort(key=lambda x: x.area, reverse=True)
    with open(os.path.join(storedir,'plottable_landmasses.pickle'), 'wb') as f:
      pickle.dump(landmasses, f, protocol=-1)

sys.argv[1:] = map(float,sys.argv[1:])

plottable_landmasses = GeneratePlottableLandmasses()
plottable_landmasses = pickle.load(open(os.path.join(storedir,'plottable_landmasses.pickle'), 'rb'))

#print(np.array([olats,olons]).transpose())
#print(np.array(TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3])).transpose())
#print(np.array(wgs_to_mercator(*TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3]))).transpose())

PlotPoints(*wgs_to_mercator(*TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3])))