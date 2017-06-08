# Declaration of variables
CC = g++
CC_FLAGS = -w

#Libraries
gdconfig=gdal-config
GDAL_LIBS=`${gdconfig} --libs` -lGeographic
GDAL_CFLAGS=`${gdconfig} --cflags` -Ilibs
ARCH_FLAGS=-march=native -mtune=native #-m32
CXXFLAGS=$(GDAL_CFLAGS) $(ARCH_FLAGS) --std=c++11 -g -Wall -ffast-math -fopenmp -DENV_LAPTOP 
 
# File names
EXEC = dgfinder.exe
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(CXXFLAGS) $(OBJECTS) -o $(EXEC) $(GDAL_LIBS)
 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(CXXFLAGS) $< -o $@ 
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)