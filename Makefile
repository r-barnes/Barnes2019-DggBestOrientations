export CXX=g++
export gdconfig=gdal-config
export GDAL_LIBS=`${gdconfig} --libs`
#export GDAL_LIBS=/z/gdal/lib/libgdal.a ../../richdem_x86_64/install/lib/libproj.a
export GDAL_CFLAGS=`${gdconfig} --cflags`
export ARCH_FLAGS=-march=native -mtune=native #-m32
export CXXFLAGS=$(GDAL_CFLAGS) $(ARCH_FLAGS) --std=c++11 -g -Wall -ffast-math -fopenmp

main:
	$(CXX) $(CXXFLAGS) -o dgfinder.exe cpp_attempt.cpp Polygon.cpp SpIndex.cpp PointCloud.cpp $(GDAL_LIBS)