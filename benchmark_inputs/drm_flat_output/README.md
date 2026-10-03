# Structured DRM boxes and station binary checks

These checks use synthetic geometry and displacements; they do not establish
physical accuracy of a simulation or calibration.

Run the input/geometry regression without MPI or CUDA:

```sh
c++ -std=c++11 -Wall -Wextra -pedantic benchmark_inputs/drm_flat_output/points_test.cc -o /tmp/drm_points_test
/tmp/drm_points_test
c++ -std=c++11 -Wall -Wextra -pedantic benchmark_inputs/drm_flat_output/coordinates_test.cc -o /tmp/drm_coordinates_test
/tmp/drm_coordinates_test /tmp/drmcoordinates-test.0
```

The coordinate regression runs the production ASCII writer. It checks exact
float64 round trips for rotated/fractional grids with synthetic terrain/depth
shifts, binary sample order (including duplicate face edges), locale-independent
formatting, truncation of an older file, and failed output paths. On Linux it
also checks write/close failure with `/dev/full`.

Run the actual MPI-IO writer on a machine with MPI (CUDA is not required):

```sh
mpicxx -std=c++11 -Wall -Wextra -pedantic benchmark_inputs/drm_flat_output/mpi_writer_test.cc -o /tmp/drm_mpi_test
mpirun -np 1 /tmp/drm_mpi_test /tmp/drm-one.bin 1
mpirun -np 4 /tmp/drm_mpi_test /tmp/drm-four.bin 4
cmp /tmp/drm-one.bin /tmp/drm-four.bin
cmp /tmp/drm-one.bin.box0 /tmp/drm-four.bin.box0
cmp /tmp/drm-one.bin.box1 /tmp/drm-four.bin.box1
cmp /tmp/drm-one.bin.box0.coords /tmp/drm-four.bin.box0.coords
cmp /tmp/drm-one.bin.box1.coords /tmp/drm-four.bin.box1.coords
```

This checks an analytical payload in global point order, interleaved rank
ownership, an empty rank, multiple batches, a partial final batch, and
truncation of an older file. It also writes two differently sized/rotated boxes
simultaneously with different output rates, checks every XYZ value against an
affine field, and verifies their independent file sizes and partial batches.
Rank zero exports the shifted sample coordinates through the production ASCII
writer; reading those rows back must recover every binary sample's location.
Both runs must produce identical bytes for all five files. The geometry test
checks the full row/column ordering of all ten faces, repeated shared edges,
center offsets, rotations, exterior-to-interior dimensions, and table parsing.

Both MPI executables require `OUTPUT_PATH EXPECTED_RANKS`, with `EXPECTED_RANKS`
equal to `1` or `4`. They compare it with `MPI_COMM_WORLD` before opening output
and abort on a mismatch, even if compiled with `NDEBUG`. A launcher that starts
four independent singleton worlds must not pass as a four-rank test. Successful
setup prints the verified world size once from rank zero. The output path and
the DRM `.box0`, `.box1`, and `.coords` companions are overwritten by valid runs.

The GitHub Actions job uses Ubuntu 24.04 with `openmpi-bin`/`libopenmpi-dev`,
explicit `mpic++.openmpi` and `mpirun.openmpi` commands, `--oversubscribe` for
four ranks on a small runner, and the Open MPI `ob1`/`self,tcp` host transports.
Use a compatible MPI compiler and launcher on other systems. Ubuntu's MPICH
4.2.0-5build3 package links an external PMIx library incompatible with its
Hydra launcher; see [MPICH issue #7064](https://github.com/pmodels/mpich/issues/7064#issuecomment-2301026290).
The CI job also deliberately launches each test with one rank while requiring
four, checks the specific failure diagnostic, and verifies that an existing
output file remains unchanged. It then checks the full one-/four-rank payloads.
Earlier green jobs without a world-size check did not establish multi-rank
coverage: independent singleton processes could race on identical output files.

Test the production station interpolation, sampling dispatch, and close paths
against an analytical field without a CUDA build:

```sh
python3 benchmark_inputs/drm_flat_output/extract_station_test.py /tmp/station_functions.inc
c++ -std=c++11 -Wall -Wextra -pedantic -I/tmp benchmark_inputs/drm_flat_output/station_test.cc -o /tmp/station-test
/tmp/station-test
```

The harness extracts the actual station functions from `psolve.cc` and the
tetrahedral displacement routine from `topography.cc`. It checks default ASCII
and binary output for regular hexahedra, cut-cell tetrahedra, a supplied acoustic
displacement state, duplicate station locations, and skipped sampling steps.
Binary mode leaves text file handles null and checks that ancillary output is
not called. This does not test acoustic reconstruction itself or GPU transfers.
The host version uses an in-memory sink. To also test the production MPI writer
with interleaved ownership, an empty rank, and OpenMP interpolation:

```sh
mpicxx -std=c++11 -Wall -Wextra -pedantic -fopenmp -DSTATION_TEST_MPI -I/tmp benchmark_inputs/drm_flat_output/station_test.cc -o /tmp/station-mpi-test
OMP_NUM_THREADS=2 mpirun -np 1 /tmp/station-mpi-test /tmp/station-one.bin 1
OMP_NUM_THREADS=2 mpirun -np 4 /tmp/station-mpi-test /tmp/station-four.bin 4
cmp /tmp/station-one.bin /tmp/station-four.bin
```

Before production use, build Hercules on the target CUDA/MPI system and run
a small existing elastic case twice with the same `output_stations` and
`output_stations_print_rate`: once with the default `station_binary=no`, then
with `station_binary=yes` and `station_output_binary_file` set. Compare every
XYZ binary snapshot with the corresponding ASCII station displacement columns
(allow the station ASCII formatting error). Repeat for the classic and v3
communication paths and the supported host/CUDA-aware MPI configurations;
compare one-rank and multiple-rank output. Check the zero initial snapshot,
last partial batch, and a checkpoint restart into a **different** output path.
For a terrain or poroelastic application, add the corresponding station
comparison before claiming that combination validated.
Also test two `output_drm_boxes` rows together with station binary output. Compare
each generated face at coincident station locations, including common face
edges, a nonzero strike, and the interior/exterior bottom separation. The
`drmcoordinates.N` rows must match the solver interpolation points, including
building/topography shifts and repeated shared edges/corners. Verify the ASCII
and binary pair for each box and restart segment. The
synthetic MPI regression does not compile or execute the CUDA solver path.
