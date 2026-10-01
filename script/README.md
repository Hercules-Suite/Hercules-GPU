# Hercules input and output utilities

Run these commands from the repository root, or use an absolute script path
from a case directory. Input and output paths resolve from the current working
directory. Each script provides `--help`. Geometry and saved-snapshot intervals
are supplied manually; these utilities do not parse `parameters.in`.

| Script | Purpose | Dependencies |
| --- | --- | --- |
| [convert_3d_velocity_model.py](convert_3d_velocity_model.py) | Convert a CSV into count-prefixed material files | Python standard library |
| [validate_3d_velocity_model.py](validate_3d_velocity_model.py) | Check counts, profiles, material ranges, and optional Q files | NumPy |
| [generate_topography.py](generate_topography.py) | Generate the documentation's Gaussian hill | NumPy |
| [validate_topography.py](validate_topography.py) | Check spacing, elevations, base coordinate, and perimeter relief | NumPy |
| [taper_topography.py](taper_topography.py) | Taper relief to a flat collar before lateral PML boundaries | Python standard library |
| [read_stations.py](read_stations.py) | Show the station header, time range, and peak displacement | NumPy |
| [read_planes.py](read_planes.py) | Inspect plane layout and scan all displacements for finite values and peak vertical motion | NumPy |
| [compare_plane_outputs.py](compare_plane_outputs.py) | Compare every XYZ displacement in two completed sets of plane files, with bounded memory | NumPy |
| [diagnose_plane_output.py](diagnose_plane_output.py) | Locate displacement peaks, extract full point histories, and join optional material/CPML metadata | NumPy |
| [read_drm.py](read_drm.py) | Inspect ten structured DRM box faces or single-file station output, time range, and peak XYZ displacements | Python standard library |
| [read_mesh.py](read_mesh.py) | Check paired mesh coordinates, material records, and element IDs | NumPy |
| [animations_mpi.py](animations_mpi.py) | Render vertical-velocity PNG frames from plane displacements | NumPy, xarray, PyGMT/GMT, joblib |

## Prepare material and topography files

```bash
python3 script/convert_3d_velocity_model.py model.csv --output-dir input/materialfiles
python3 script/validate_3d_velocity_model.py input/materialfiles

python3 script/generate_topography.py \
  --lx 40000 --ly 40000 --spacing 100 --height 800 --sigma 3000 \
  --center-x 20000 --center-y 20000 --out input/topography/topography.in
python3 script/validate_topography.py input/topography/topography.in \
  --lx 40000 --ly 40000 --base 1250

mkdir -p input/topography_pml
python3 script/taper_topography.py input/topography/topography.in \
  --out input/topography_pml/topography.in --lx 40000 --ly 40000 \
  --buffer 6000 --flat-margin 1500 --pml-thickness 780 --z-ref auto
```

The hill and mesh dimensions above are the documentation's illustrative case.
Choose dimensions, base coordinate, and PML settings for your own case.
Generation creates output directories and refuses to replace an existing file
unless `--force` is supplied. Tapering requires the output directory to exist
and likewise uses `--force` for replacement. The CSV converter retains its
existing behavior of replacing material files in the chosen directory.

Material profiles use **Y outer, X inner** ordering. Topography uses **X outer,
Y inner**, with `flat_index = x_index * n_east + y_index`, as in the solver's
`point_elevation` routine. The relocated taper now follows this topography
order; the previous implementation interpreted rectangular grids incorrectly.
These validation utilities check file structure and numeric ranges, not the
physical suitability of a model.

## Inspect solver output

```bash
python3 script/read_stations.py output/stations/station.0
python3 script/read_planes.py output/planes/planedisplacements.0 \
  --nx 101 --ny 81 --dt 0.02
python3 script/read_mesh.py output/mesh/mesh_coordinates.0 output/mesh/mesh_data.0
```

`read_stations.py` reads time and the first three displacement columns and
prints the header for interpreting optional velocity, acceleration, pressure,
and stress columns. It does not assign fixed offsets to those optional fields.

Plane data have shape `(time, along-strike point, down-dip point, XYZ)` in C
order. `--dt` is the interval between **saved snapshots**, in seconds; use
`--first-time` when the first stored snapshot is not at zero. The plane reader
scans the entire file one snapshot at a time using a read-only memory mapping.
Binary readers assume native byte order matching the machine that wrote the
files. Mesh records are packed: 32 bytes per corner and 20 bytes per material
record, with eight corner records per element.

Compare completed legacy and MPI-IO output from a case's `outputfiles` directory:

```bash
python3 /path/to/Hercules_GPU/script/compare_plane_outputs.py planes planes_mpiio \
  --shapes 361x361 71x71 71x71 71x71 71x71 71x71 71x71 \
  --expected-snapshots 1000
```

This scans all seven files and defaults to exact numerical equality. Optional
`--atol` (metres) and `--rtol` apply a per-value acceptance criterion; NaN and
infinity always fail. The report distinguishes byte identity from numerical
agreement and locates maximum errors. It does not load whole files into RAM.
See [comparison criteria and options](../doc/Extended%20README.md#compare-plane-output-backends).

Locate large surface-plane motion and extract small reports from a completed run:

```bash
python3 /path/to/Hercules_GPU/script/diagnose_plane_output.py planes/planedisplacements.0 \
  --nx 361 --ny 361 --dt 0.1 --expected-snapshots 1000 \
  --spacing 250 250 --xy-origin 0 0 --domain 90000 90000 \
  --cpml-thickness 439.453125 --points 0,183 88,0 1,90 \
  --output-dir plane_diagnostic
```

The geometry options assume a horizontal, unrotated plane in solver X/Y metres.
This scans twice and writes `peak_points.csv`, `histories.csv`,
`frame_peaks.csv`, and `summary.json` into a new directory. The optional
`--metadata /path/to/planes/plane_point_metadata.*.csv --plane-id 0` joins
containing-element materials and actual CPML flags exported at setup by
`HERCULES_PLANE_METADATA=1`. Use metadata from the same mesh and plane
configuration. Without it, acoustic/solid/CPML classification stays unknown;
the geometric boundary band is only a candidate region. See the
[diagnostic workflow](../doc/Extended%20README.md#diagnose-large-plane-displacements).

Each generated DRM box writes one `drmdisplacements.N` containing native
float64 triples, ordered `(time, concatenated face sample, XYZ)`. Pass exterior
dimensions and spacing to inspect its ten structured faces:
`python3 script/read_drm.py output/drmdisplacements.0 --box 100 100 50 1 --dt 0.01`.
Face order is exterior front/back/left/right/bottom, then the five interior
faces. Shared edges/corners repeat to keep each grid rectangular. The helper
`box_layout(lx, ly, depth, spacing)` returns `(name, rows, columns, point_offset)`
tuples for slicing these faces. Center and strike affect coordinates, not shapes;
all output vectors retain global solver XYZ components. The saved interval is
`output_drm_boxes_print_rate * simulation_delta_time_sec`.

Each box also writes `drmcoordinates.N` once at setup. This ASCII file has
comment headers (`#`) followed by three columns `X_m Y_m Z_m`, in solver
north/east/down metres after rotation and building/topography depth shifts.
Read it with `xyz = numpy.loadtxt("output/drmcoordinates.0", ndmin=2)`.
Row i matches binary sample i, including repeated shared edges/corners;
use the same `box_layout` slices to reshape coordinates into face grids.
Coordinates use 17 significant digits and are overwritten on every run.

With `station_binary=yes`, `station_output_binary_file` stores
`(time, station, XYZ)` in `output_stations` row order. Run
`python3 script/read_drm.py output/station_displacements.0 --points 2 --dt 0.01`
for two stations, using `output_stations_print_rate * simulation_delta_time_sec`
as `--dt`. Neither format has headers or embedded point lists. Retain the
input tables. Use `--first-time` for restart segments. Read completed output;
DRM boxes buffer snapshots according to `drm_output_batch_size`, while station
binary output writes every saved snapshot.

For rendering and `mpirun` affinity settings, see the dedicated
[plane velocity animation guide](../doc/plane_velocity_animations.md).
