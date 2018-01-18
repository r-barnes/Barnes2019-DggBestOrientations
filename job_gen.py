#!/usr/bin/env python3
import itertools

job="""#!/bin/bash  
#SBATCH --job-name="dgg-{shape}-{findwhat}"
#SBATCH --output="dgg-{shape}-{findwhat}.%j.%N.out"  
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH -t {timelen}

cd /home/rbarnes1/dgg_best_poles2

export LD_LIBRARY_PATH=$HOME/os/anaconda3/lib:$HOME/os/lib:$LD_LIBRARY_PATH

ulimit -c unlimited

./dgfinder.exe optimize ellipsoidal {shape} {findwhat}
"""

shapes = [
  "regular_icosahedron",
  "regular_dodecahedron",
  "regular_tetrahedron",
  "regular_octahedron",
  "cuboctahedron"
]

goals = [
  ("min_avg_dist","00:20:00"),
  ("max_avg_dist","00:20:00"),
  ("min_edge","20:00:00"),
  ("max_edge","20:00:00")
]

for s,g in itertools.product(shapes, goals):
  filename = s+"_"+g[0]+".job"
  with open(filename, 'w') as fout:
    fout.write(job.format(
      shape    = s,
      findwhat = g[0],
      timelen  = g[1]
    ))
  print("sbatch " + filename)