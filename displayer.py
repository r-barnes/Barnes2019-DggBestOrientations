#!/usr/bin/env python3
import pyproj
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import pickle
import sys

unprojected_landmassfile = 'data/land-polygons-complete-4326/land_polygons.shp'
plottable_landmassfile   = 'data/simplified-land-polygons-complete-3857/simplified_land_polygons.shp'
landmassfile             = 'data/land-polygons-split-3857/land_polygons.shp'

#storedir = '/home/rbarnes1/scratch/dgg_best'
storedir  = 'temp'

def PlotPoints(lat,lon):
  """Plot points along with some of the largest land masses, for content"""
  import time, calendar, datetime, numpy
  from mpl_toolkits.basemap import Basemap
  import matplotlib.pyplot as plt
  # draw map with markers for float locations
  m = Basemap(projection='mill',lon_0=0)
  m.drawmapboundary(fill_color='#99ffff')
  m.fillcontinents(color='#cc9966',lake_color='#99ffff')
  x,y=m(lon,lat)
  m.scatter(x,y,
            marker='o',color='k',
            zorder=10)
  plt.title('Hammer projection, data on top',fontsize=12)
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

sys.argv[1:] = map(float,sys.argv[1:])

#print(np.array([olats,olons]).transpose())
#print(np.array(TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3])).transpose())
#print(np.array(wgs_to_mercator(*TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3]))).transpose())

PlotPoints(*TransformLatLon(olats,olons,sys.argv[1],sys.argv[2],sys.argv[3]))