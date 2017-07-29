# Declaration of variables
CC = g++
CC_FLAGS = -w

#Libraries
gdconfig=gdal-config
GDAL_LIBS=`${gdconfig} --libs` -lGeographic
GDAL_CFLAGS=`${gdconfig} --cflags`
ARCH_FLAGS=-march=native -mtune=native #-m32
CXXFLAGS=--std=c++14 -g -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual -fopenmp -DENV_LAPTOP -DCOARSE_RESOLUTION #-Wmissing-declarations 

#ENABLE THIS LINE FOR PRODUCTION
#FLAGS=$(GDAL_CFLAGS) $(ARCH_FLAGS) $(CXXFLAGS) -O3 -DDOCTEST_CONFIG_DISABLE

#ENABLE THIS LINE TO PERFORM CODE TESTS
FLAGS=$(GDAL_CFLAGS) $(ARCH_FLAGS) $(CXXFLAGS) --coverage -O0

# File names
EXEC = dgfinder.exe
#SOURCES = $(subst cpp_attempt.cpp,,$(wildcard *.cpp))
SOURCES = $(subst test.cpp,,$(wildcard *.cpp))
#SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(FLAGS) $(OBJECTS) -o $(EXEC) $(GDAL_LIBS)
 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(FLAGS) $< -o $@ 
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)