cmake_minimum_required (VERSION 3.9)

project(dgg_poi_optimizer LANGUAGES CXX)

set(CMAKE_MODULE_PATH /usr/share/cmake/geographiclib/)

find_package(GeographicLib REQUIRED)
find_package(GDAL REQUIRED)
find_package(OpenMP)

if(NOT OpenMP_CXX_FOUND)
  message(WARNING "Compiling without OpenMP!")
endif()

add_executable(dgfinder.exe
  GeoStuff.cpp
  IndexedShapefile.cpp
  main.cpp
  Orientation.cpp
  OrientationIndex.cpp
  PointCloud.cpp
  Point.cpp
  Polygon.cpp
  Progress.cpp
  random.cpp
  Solid.cpp
  SpIndex.cpp
  Timer.cpp
)
target_compile_features(dgfinder.exe PRIVATE cxx_std_14)
target_compile_options(dgfinder.exe PRIVATE -march=native -Wall -Wextra -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual)
target_compile_definitions(dgfinder.exe PRIVATE ENV_LAPTOP COARSE_RESOLUTION DOCTEST_CONFIG_DISABLE)
target_link_libraries(dgfinder.exe PRIVATE OpenMP::OpenMP_CXX ${GDAL_LIBRARY} ${GeographicLib_LIBRARIES})
target_include_directories(dgfinder.exe PRIVATE ${GDAL_INCLUDE_DIR})

add_executable(poi.exe
  GeoStuff.cpp
  IndexedShapefile.cpp
  Orientation.cpp
  OrientationIndex.cpp
  Point.cpp
  PointCloud.cpp
  pole_of_inaccessibility.cpp
  Polygon.cpp
  Progress.cpp
  random.cpp
  Solid.cpp
  SpIndex.cpp
  Timer.cpp
)
target_compile_features(poi.exe PRIVATE cxx_std_14)
target_compile_options(poi.exe PRIVATE -march=native -Wall -Wextra -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual)
target_compile_definitions(poi.exe PRIVATE ENV_LAPTOP COARSE_RESOLUTION DOCTEST_CONFIG_DISABLE)
target_link_libraries(poi.exe PRIVATE OpenMP::OpenMP_CXX ${GDAL_LIBRARY} ${GeographicLib_LIBRARIES})
target_include_directories(poi.exe PRIVATE ${GDAL_INCLUDE_DIR})


#ENABLE THIS LINE TO PERFORM CODE TESTS
#FLAGS=$(GDAL_CFLAGS) $(ARCH_FLAGS) $(CXXFLAGS) --coverage -O0
