2017-DggBestOrientations
========================

**Title of Manuscript**:
Optimal orientations of discrete global grids

**Authors**: Richard Barnes

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**DOI Number of Manuscript**
TODO

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2017-DggBestOrientations)
 * [Journal's GitHub Repository](TODO)

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

    make

To compile the program for debugging, please comment/uncomment the indicated
lines in `Makefile`.

The result of either compilation is a program called `dgfinder.exe`.

Running this program will either run extensive code tests or perform the optimal
orientation search.
