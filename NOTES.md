perf record ./src/salamander.exe parms/z_example_param.param

http://webpages.sou.edu/~sahrk/dgg/orientation/dymorient.html
The Dymaxion orientation of R. Buckminster Fuller can be specified by placing one icosahedron vertex at 5.245390W longitude, 2.3008820N latitude, and placing an adjacent vertex at an azimuth of 7.466580 from the first vertex. Note that this placement is the only known icosahedron placement with no icosahedron vertices falling on land. 


http://www.rwgrayprojects.com/rbfnotes/maps/graymap4.html
Fuller's World Map: Coordinates


Below is a table I composed for the positioning of the icosahedron as used to create Fuller's icosahedron based world map. These coordinates are very slightly modified to those I first saw in an IBM Technical Report by Mr. Scott. The modifications mentioned here were to make the distances between connected vertices all the same. The coordinates I saw in the IBM TR were off ever so slightly and may have been due to rounding of the numbers. You will also have to adjust the coordinates listed below if you wish to use them and you need them to be more accurate than I list here. A year or so later, I saw some notes by Chris Kitrick at the BFI which also had the coordinates for the 12 vertices. I do not know what coordinates Bucky Fuller really used, nor the coordinates which Shoji Sadao used, but they must be very close to those I list below.

Fuller's Positioning of the Icosahedron
VERTEX  X Y Z LONGITUDE LATITUDE
1   0.420152  0.078145   0.904083  +10.53620 +64.7
2   0.995005  -0.091348  0.040147   -5.24539  +2.300882
3   0.518837  0.835420   0.181332  +58.15771 +10.447378
4   -0.414682 0.655962   0.630676 +122.3     +39.1
5   -0.515456 -0.381717  0.767201 -143.47849 +50.103201
6   0.355781  -0.843580  0.402234  -67.13233 +23.717925
7   0.414682  -0.655962 -0.630676  -57.7     -39.1
8   0.515456  0.381717  -0.767201  +36.5215  -50.1032
9   -0.355781 0.843580  -0.402234 +112.86767 -23.717925
10  -0.995009 0.091348  -0.040147 +174.7546   -2.3009
11  -0.518837 -0.835420 -0.181332 -121.84229 -10.447345
12  -0.420152 -0.078145 -0.904083 -169.4638  -64.7
Go To: Table of Contents or Home Page
Usage Note: My work is copyrighted. You may use my work but you may not include my work, or parts of it, in any for-profit project without my consent.

rwgray@rwgrayprojects.com


perf record -e LLC-loads,LLC-load-misses yourExecutable


 5.245390W longitude, 2.3008820N latitude, and placing an adjacent vertex at an azimuth of 7.466580 from the first vertex. Note that this placement is the only known icosahedron placement with no icosahedron vertices falling on land. 






  //https://www.mapbox.com/blog/cheap-ruler/
  //http://www.focusonmath.org/sites/focusonmath.org/files/assets/MT2004-08-20a%281%29.pdf
  //From: https://www.gpo.gov/fdsys/pkg/CFR-2005-title47-vol4/pdf/CFR-2005-title47-vol4-sec73-208.pdf. Valid for distances <295 miles
  // const double ml      = (lat1+lat2)/2;
  // const double cs      = std::cos(ml);
  // const double cs2     = 2*cs*cs  - 1;
  // const double cs3     = 2*cs*cs2 - cs;
  // const double cs4     = 2*cs*cs3 - cs2;
  // const double cs5     = 2*cs*cs4 - cs3;
  // const double kpd_lat = 111.13209    - 0.56605*cs2 + 0.00120*cs4;
  // const double kpd_lon = 111.41513*cs - 0.09455*cs3 + 0.00012*cs5;
  // const double ns      = kpd_lat*(lat1-lat2);
  // const double ew      = kpd_lon*(lon1-lon2);
  // return std::sqrt(ns*ns+ew*ew);

  //Law of Cosines method
  //const double Rearth = 6371; //km
  //return Rearth*std::acos( std::sin(lat1)*std::sin(lat2) + std::cos(lat1)*std::cos(lat2)*std::cos(lon2-lon1) );









set env LD_PRELOAD /usr/lib64/libsqlite3.so.0

set env LD_LIBRARY_PATH /share/apps/compute/gcc-6.2.0/lib64/:/home/rbarnes1/os/anaconda3/lib:/opt/boost/intel/mvapich2_ib/lib/:/home/rbarnes1/os/lib:/home/rbarnes1/os/lib:/home/rbarnes1/lib::/usr/lib64:/home/rbarnes1/os/anaconda3/lib





g++6 -I/home/rbarnes1/os/include -Ilibs -I/home/rbarnes1/os/anaconda3/include -march=native -mtune=native --std=c++11 -O3 -g -Wall -ffast-math -fopenmp -o dgfinder.exe cpp_attempt.cpp  GeoStuff.cpp  Icosa.cpp  PointCloud.cpp  Point.cpp  Polygon.cpp  SpIndex.cpp POI.cpp -L/home/rbarnes1/os/anaconda3/lib -lgdal -L/home/rbarnes1/os/anaconda3/lib -lgeos_c -L/home/rbarnes1/os/anaconda3/lib -lxerces-c -lpthread -L/home/rbarnes1/os/anaconda3/lib -lnetcdf -L/home/rbarnes1/os/anaconda3/lib -lhdf5 -L/home/rbarnes1/os/anaconda3/lib -lmfhdf -ldf -L/home/rbarnes1/os/anaconda3/lib -lgeotiff -L/home/rbarnes1/os/anaconda3/lib -ltiff -lz -lpthread -lm -lrt -ldl -L/home/rbarnes1/os/anaconda3/lib -lcurl -L/home/rbarnes1/os/anaconda3/lib -lkea -L/home/rbarnes1/os/anaconda3/lib -lhdf5 -lhdf5_hl -lhdf5_cpp -ljpeg -lgeos -lproj -llzma -lssl -lcrypto -I/opt/boost/intel/mvapich2_ib/include/ -DENV_XSEDE -lGeographic
