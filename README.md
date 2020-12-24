Source Code: Optimal orientations of discrete global grids
==========================================================

**Title of Manuscript**:
Optimal orientations of discrete global grids

**Authors**: Richard Barnes

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**DOI Number of Manuscript**
TODO

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2019-DggBestOrientations)

This repository contains the algorithms and program described in the manuscript
above, along with information on acquiring the various datasets used, and code
to perform correctness tests. OpenMP is used to parallelize the code.



Abstract
--------

Spatial analyses involving binning require that every bin have the same area,
but this is impossible using a rectangular grid laid over the Earth or over any
projection of the Earth. Discrete global grids use hexagons, triangles, and
diamonds to overcome this issue, overlaying the Earth with equally-sized bins.
Such discrete global grids are formed by tiling the faces of a polyhedron.
However, the orientations of these polyhedra have been chosen to satisfy only
simple criteria such as equatorial symmetry or minimizing the number of vertices
intersecting landmasses. Here, efficient algorithms are used to optimize the
orientation of polyhedra with respect to a number of quantities of interest. The
optimized orientations are used to reenvision Fuller's Dymaxion map.



Compilation
-----------

Clone this repo using

    git clone git@github.com:r-barnes/2017-DggBestOrientations.git

You must also obtain certain prerequisites:

    sudo apt install TODO gdal geographiclib boost-geometry

[Anaconda3](#TODO) was used to install the aforementioned libraries on XSEDE.
The file `Makefile.xsede` should provide useful hints for getting everything
working; however, given the uniqueness of potential target environments for the
code, I do not give details here. While I am happy to provide support for my
code, it is best to work with your supercomputer's help staff to get this
compiled on your machine of choice.

To compile the program for production run:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

The result of either compilation is a program called `dgfinder.exe`.

Running this program will either run extensive code tests or perform the optimal
orientation search.



Converting the ADD data
=======================

The ADD high-res data can be acquired from here:

https://data.bas.ac.uk/collections/e74543c0-4c4e-4b41-aa33-5bb2f67df389/

The exact dataset used is here:

https://data.bas.ac.uk/items/ad7d345a-0650-4f44-b7eb-c48e1999086b/

It must be converted into a WGS84 projection and appropriate datasets selected,
as follows:
```bash
ogr2ogr -progress -t_srs '+proj=longlat +datum=WGS84 +no_defs' -where "surface='ice coastline' OR surface='grounding line'" -f "ESRI Shapefile" add_ice_coast_and_grounding_line.shp add_coastline_high_res_line_v7.3.shp
```
