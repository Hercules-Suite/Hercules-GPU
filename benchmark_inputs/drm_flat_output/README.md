# Single-file DRM recorder checks

These checks use synthetic geometry and displacements; they do not establish
physical accuracy of a simulation or calibration.

Run the input/geometry regression without MPI or CUDA:

```sh
c++ -std=c++11 -Wall -Wextra -pedantic benchmark_inputs/drm_flat_output/points_test.cc -o /tmp/drm_points_test
/tmp/drm_points_test
# Optionally check a real CSV with the production point parser:
/tmp/drm_points_test /path/to/drm_points_geodetic_1m.csv
```

Run the actual MPI-IO writer on a machine with MPI (CUDA is not required):

```sh
mpicxx -std=c++11 -Wall -Wextra -pedantic benchmark_inputs/drm_flat_output/mpi_writer_test.cc -o /tmp/drm_mpi_test
mpirun -np 1 /tmp/drm_mpi_test /tmp/drm-one.bin
mpirun -np 4 /tmp/drm_mpi_test /tmp/drm-four.bin
cmp /tmp/drm-one.bin /tmp/drm-four.bin
```

This checks an analytical payload in global point order, interleaved rank
ownership, an empty rank, multiple batches, a partial final batch, and
truncation of an older file. Both runs must produce identical bytes.

Before production use, build Hercules on the target CUDA/MPI system and run
a small existing elastic case with a few `output_stations` rows also supplied
to `drm_points_file`, using equal station/DRM print rates. Compare every XYZ
displacement snapshot with the corresponding station displacement columns
(allow the station ASCII formatting error). Repeat for the classic and v3
communication paths and the supported host/CUDA-aware MPI configurations;
compare one-rank and multiple-rank output. Check the zero initial snapshot,
last partial batch, and a checkpoint restart into a **different** output path.
For a terrain or poroelastic application, add the corresponding station
comparison before claiming that combination validated.
