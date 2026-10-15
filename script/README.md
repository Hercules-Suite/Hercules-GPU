# Hercules script usage guide

Run the examples from the repository root, or replace `script/` with its absolute
path when working in a case directory. Relative input/output paths resolve from
the current working directory. Every top-level script accepts `--help`:

```bash
python3 script/read_planes.py --help
```

Use a Python 3 environment with the dependencies listed below. Most readers need
NumPy (`python3 -m pip install numpy`). Animation rendering additionally needs a
working PyGMT/GMT installation, xarray, and joblib; its `--check-only` mode needs
only NumPy. No utility below launches a Hercules simulation. Use the repository's
[`run_phaseA_phaseB.sh`](../run_phaseA_phaseB.sh) for meshing and solving.

These tools do not parse `parameters.in`. Supply the geometry and timing from
the run that produced the files. Solver X is north, Y is east, and Z/depth is
positive downward. Coordinates and displacements are in metres unless stated
otherwise. Topography values are elevations, positive upward from the reference
datum. Binary readers use native byte order unless an `--endian` option is
listed. Read completed output files; buffering and incomplete runs can leave
partial files or fewer complete snapshots than expected.

## Script index

| Script | Purpose | Dependencies |
| --- | --- | --- |
| [convert_3d_velocity_model.py](#convert_3d_velocity_modelpy) | Convert CSV to material files | Python standard library |
| [validate_3d_velocity_model.py](#validate_3d_velocity_modelpy) | Check material counts, profiles, ranges, and Q files | NumPy |
| [generate_topography.py](#generate_topographypy) | Generate a Gaussian hill | NumPy |
| [validate_topography.py](#validate_topographypy) | Check topography spacing, elevations, datum, and perimeter | NumPy |
| [taper_topography.py](#taper_topographypy) | Make a flat collar before lateral CPML | Python standard library |
| [read_stations.py](#read_stationspy) | Summarize an ASCII station history | NumPy |
| [read_planes.py](#read_planespy) | Inspect a binary displacement plane | NumPy |
| [read_mesh.py](#read_meshpy) | Inspect paired mesh coordinate/material files | NumPy |
| [animations_mpi.py](#animations_mpipy) | Render vertical-velocity PNG frames | NumPy; rendering: PyGMT/GMT, xarray, joblib |
| [read_drm.py](#read_drmpy) | Inspect structured DRM or flat station displacement binaries | Python standard library |

The [local investigation helpers](#local-investigation-helpers) section also
documents four ignored scripts and the ignored `tests/` directory retained in
the maintainer's workspace. They are **not distributed in GitHub checkouts**.

## Input preparation

### convert_3d_velocity_model.py

Convert a row-wise CSV into the count-prefixed files read by Hercules:

```bash
python3 script/convert_3d_velocity_model.py model.csv \
  --output-dir inputfiles/materialfiles
python3 script/validate_3d_velocity_model.py inputfiles/materialfiles
```

The positional CSV and `--output-dir` are required; there are no other controls.
Required columns are `x,y,z,vs,vp,rho`. X/Y/Z use metres, Vs/Vp use m/s, and
**rho uses tonne/m³**, as in the 3-D material input format (not the kg/m³ used in
many solver outputs). Optional `qs,qp` must both be present or both absent and
contain positive dimensionless values. For example, this complete 2 × 2 grid
has two depth samples per column; values are illustrative:

```csv
x,y,z,vs,vp,rho
0,0,0,500,1500,2.0
0,0,1000,1000,2500,2.2
500,0,0,500,1500,2.0
500,0,1000,1000,2500,2.2
0,500,0,500,1500,2.0
0,500,1000,1000,2500,2.2
500,500,0,500,1500,2.0
500,500,1000,1000,2500,2.2
```

X and Y must each have at least two uniformly spaced coordinates, with every
X/Y combination present. Profiles may have different sample counts; depths
must be nonnegative and unique within each profile. Require finite values,
Vs >= 0, Vp > 0, and rho > 0. Input rows can be unordered.

The converter creates the output directory and writes `x.in`, `y.in`,
`index.in`, `z.in`, `vs.in`, `vp.in`, `rho.in`, and optional `qs.in`/`qp.in`.
Every file starts with its value count. Profiles use **Y outer, X inner**
ordering; `index.in` holds zero-based offsets into concatenated depth/material
arrays. Values are sorted within each profile by increasing depth.

Existing generated files are **overwritten without a `--force` option**.
Unrelated files, including old Q files when the new CSV has no Q columns, are
not removed. Use a fresh output directory to avoid mixing models. The command
prints the profile and sample counts; it does not configure `parameters.in`.

The converter and validator hold their input arrays in memory. Use a streaming
exporter and bounded-memory verification for models with billions of samples;
these convenience scripts are intended for models that fit available memory.
Hercules supports signed 64-bit sample counts and profile offsets, but the
complete material arrays are replicated on each solver rank. See the
[material memory estimate](../doc/Extended%20README.md#3d-velocity-model-files).

### validate_3d_velocity_model.py

```bash
python3 script/validate_3d_velocity_model.py inputfiles/materialfiles
```

The only argument is the directory containing the seven required material files
listed above. The validator checks finite values, count headers, increasing
uniform X/Y axes, exactly one ordered profile per X/Y location, increasing
depths, equal material-array lengths, and the same value ranges as the converter.
If either Q file exists, both must exist and have one positive Q per sample.

It prints `valid model: ... columns, ... material samples` on success and exits
nonzero on an invalid model. It does not write files, check domain coverage, or
establish the physical suitability of the material values.

### generate_topography.py

Generate a regular Gaussian-hill elevation grid:

```bash
python3 script/generate_topography.py \
  --lx 40000 --ly 40000 --spacing 100 --height 800 --sigma 3000 \
  --center-x 20000 --center-y 20000 \
  --out inputfiles/topography/topography.in
```

| Argument | Meaning / default |
| --- | --- |
| `--lx`, `--ly` | Required north/east domain lengths, m; positive |
| `--spacing` | Required common grid spacing, m; must divide both lengths |
| `--height` | Required peak elevation, m |
| `--sigma` | Required positive Gaussian width, m |
| `--center-x`, `--center-y` | Hill centre, m; defaults to the domain midpoint |
| `--out` | Required output filename |
| `--force` | Allow replacement of an existing output file; off by default |

The output directory is created automatically. The file contains one spacing
value followed by `(lx/spacing+1)*(ly/spacing+1)` elevations. It uses **X outer,
Y inner** ordering: `flat_index = x_index * n_east + y_index`. The command
prints the generated dimensions and sample count. This generates an elevation
file only; select `topographybase_zcoord` and the material model separately.
The Gaussian is not explicitly tapered to an exactly flat boundary.

### validate_topography.py

```bash
python3 script/validate_topography.py inputfiles/topography/topography.in \
  --lx 40000 --ly 40000 --base 1250
```

The file, `--lx`, `--ly`, and `--base` are required. Supply domain lengths and
`topographybase_zcoord` in metres. The validator requires finite values,
positive lengths/base/spacing, spacing that divides the lengths, the expected
elevation count, and a base strictly greater than the highest elevation.

It prints the grid size, elevation range, and perimeter elevation range, without
writing files. A nonzero perimeter range is **reported, not rejected**. Even a
zero perimeter range does not prove that the full CPML collar is flat. This
check does not validate bathymetry, material coverage, or source clearance.
`--base` is a coordinate datum inside the domain; it does not increase the
domain depth configured in `parameters.in`.

### taper_topography.py

Smoothly reduce relief to a constant elevation before entering lateral CPML:

```bash
mkdir -p inputfiles/topography_pml
python3 script/taper_topography.py inputfiles/topography/topography.in \
  --out inputfiles/topography_pml/topography.in \
  --lx 40000 --ly 40000 --buffer 6000 --flat-margin 1500 \
  --pml-thickness 780 --z-ref 0
python3 script/validate_topography.py inputfiles/topography_pml/topography.in \
  --lx 40000 --ly 40000 --base 1250
```

| Argument | Meaning / default |
| --- | --- |
| `infile` | Required input `topography.in` |
| `--out` | Output filename; default `topography_tapered.in` |
| `--lx`, `--ly` | North/east domain lengths, m; each defaults to 40000 |
| `--buffer` | Total band width measured inward from each boundary, m; 6000 |
| `--flat-margin` | Constant-elevation collar within that band, m; 1500 |
| `--pml-thickness` | Thickness used for the collar-width check only, m; 780 |
| `--z-ref` | Flat target elevation, m; `0.0`, or `auto` for the mean of the first interior ring at the taper band |
| `--force` | Permit replacing the input or an existing output; off by default |

Let `r` be distance to the nearest lateral boundary. Elevations are unchanged
for `r >= buffer`, equal to `z_ref` for `r <= flat_margin`, and follow a raised
cosine transition between those limits. Require
`0 < flat_margin < buffer < min(lx,ly)/2`. Set the flat margin at least as wide
as CPML, preferably with additional margin; the tool only warns when it is too
narrow. `--pml-thickness` does not change the solver's CPML configuration.

The output directory must already exist. The file retains its spacing and
X-outer/Y-inner ordering; elevations are written to four decimal places. Thus,
even untapered interior values may be rounded. Output reports include changed
point count, largest elevation change, target level, and perimeter flatness.
The relocated utility handles rectangular grids in the solver's ordering;
older versions interpreted those grids incorrectly.

Select the generated directory with
`topography_directory = inputfiles/topography_pml` and regenerate the mesh for
the changed geometry. Choose a reference appropriate to the case: `auto` can
move a zero-elevation water surface if nearby land raises the ring average.
Tapering this file does not modify the 3-D material water mask or seabed and
does not increase source-to-boundary clearance. See the
[CPML/topography discussion](../doc/Extended%20README.md#cpml-execution-and-topography).

## Output inspection and rendering

### read_stations.py

```bash
python3 script/read_stations.py outputfiles/stations/station.0
```

The only argument is an ASCII `station.N` file. The reader echoes leading
comment headers and loads the numeric table. It requires at least the first
four columns: time (s), then X/Y/Z displacement (m). Those columns must be finite
and time must increase strictly.

Output includes sample count, saved-time interval range, time span, column
count, and peak displacement **vector magnitude** `sqrt(ux²+uy²+uz²)` in metres.
This is not a report of three separate component peaks or peak acceleration.
Optional velocity, acceleration, pressure, and stress columns remain described
by the echoed header; the utility does not validate or summarize them. It reads
the whole table into memory and does not write output files. For flat binary
station output, use `read_drm.py --points` below.

### read_planes.py

```bash
python3 script/read_planes.py outputfiles/planes/planedisplacements.0 \
  --nx 361 --ny 361 --dt 0.1
```

The filename, `--nx`, `--ny`, and `--dt` are required. `nx` is the along-strike
point count and `ny` the down-dip point count, not element counts. `--first-time`
defaults to 0 s; set it to the first saved time for a restart segment.

`--dt` means **time between saved snapshots**, usually
`simulation_delta_time_sec * output_planes_print_rate`. For example, a 0.0005 s
solver timestep and print rate 200 give `--dt 0.1`.

The headerless binary contains native float64 displacement triples with C-order
shape `(time, nx, ny, 3)`; the down-dip index varies fastest. Each snapshot has
`nx*ny*3*8` bytes. Require positive counts/interval and at least one complete
snapshot. The read-only memory mapping scans every snapshot for finite XYZ
values. Output reports shape, time span, and peak absolute **global Z**
displacement over the entire file, including for inclined planes. No files are
written. Incorrect dimensions with the same product cannot be detected from
the binary alone; the reader also cannot detect missing complete frames without
knowing the expected count.

### read_mesh.py

```bash
python3 script/read_mesh.py \
  outputfiles/mesh_coordinates.0 outputfiles/mesh_data.0
```

Both positional files must come from the same run and MPI rank. Use the actual
directory where the solver wrote them. They are packed native-endian records:

| File | Record layout | Size |
| --- | --- | --- |
| `mesh_coordinates.RANK` | int64 element ID, float64 X/Y/Z in m | 32 bytes per corner |
| `mesh_data.RANK` | int64 element ID, float32 Vs/Vp/rho | 20 bytes per element |

There must be eight consecutive corner records for each material record, with
matching element IDs. Material velocities use m/s and density uses kg/m³.
The reader rejects incomplete records, mismatched counts/IDs, and nonfinite
coordinates/properties. It prints element/corner counts and coordinate ranges,
loads both files into memory, and writes nothing. It does not read the solver's
`mesh.e` etree/archive or plot the mesh, and it does not reject finite sentinel
material values solely because they are negative.

### animations_mpi.py

Check the plane layout first, then render a PNG sequence:

```bash
python3 script/animations_mpi.py \
  --plane-file outputfiles/planes/planedisplacements.0 \
  --nx 361 --ny 361 --dx 250 --dy 250 --dt 0.1 --check-only

python3 script/animations_mpi.py \
  --plane-file outputfiles/planes/planedisplacements.0 \
  --nx 361 --ny 361 --dx 250 --dy 250 --dt 0.1 \
  --workers 4 --output-dir Animations_velocity_run1
```

Required: `--plane-file`, `--nx`, `--ny`, and saved-snapshot interval `--dt`.
Use at least two points per axis and two snapshots. The binary layout matches
`read_planes.py`: headerless native float64 in C order, shaped
`(NT, nx, ny, 3)`, with XYZ components contiguous and down-dip points varying
faster than along-strike points. For spacing/count pairs `500 201` and
`325 126`, supply `--nx 201 --ny 126 --dx 500 --dy 325`; the shape is
`(NT, 201, 126, 3)`. Read the file on a machine with matching byte order.

Rendering assumes **strike=0 (modulo 360), dip=0**; the script cannot verify
those angles from the file. Plots show solver Y rightward and X upward, as
offsets from the plane origin in km, without a transpose or mirror. Geographic
origin conversion and source-location consistency remain separate checks.

Set `--dt` to `simulation_delta_time_sec * output_planes_print_rate`, in
seconds; no additional print-rate factor is applied. For restarted output,
set `--first-time` to the time of the first snapshot actually stored. This
changes time labels; the binary has no timestamps and the script cannot detect
restarts. Geometry and timing are supplied manually, not read from
`parameters.in`.

| Optional control | Default / effect |
| --- | --- |
| `--dx`, `--dy` | 500 and 325 m; along-strike/down-dip spacing |
| `--first-time` | 0 s; time label for the first stored snapshot |
| `--output-dir` | `Animations_velocity`; creates the directory |
| `--workers` | Workers per rank; `FRAME_WORKERS_PER_RANK`, else `FRAME_WORKERS_PER_NODE`, else 1 |
| `--velocity-limit` | Symmetric colour limit in cm/s; `VELOCITY_LIMIT_CM_S`, else 10 |
| `--progress-every` | Progress interval in frames; `FRAME_PROGRESS_EVERY`, else 100; 0 disables |
| `--dpi` | PNG resolution, 300 |
| `--plot-upsample` | Display refinement per axis, 4; use 1 for native samples |
| `--plot-interpolation` | `linear` (default) or `cubic` |
| `--title` | `Hercules vertical velocity` |
| `--check-only` | Check layout and velocities at first/middle/last frames without plotting |

Velocity is obtained by differentiating saved global Z displacements, using
central differences internally and one-sided differences at endpoints, then
converting to **cm/s, positive downward**. Reported values are not clipped;
colours saturate at the display limit. Progress reports include extrema and
the percentage beyond the colour limits. Differentiating decimated snapshots
can differ from station velocities saved at every solver step, especially at
higher frequencies.

Output is `frame_00000.png`, `frame_00001.png`, etc., not an encoded video.
Existing matching frames are replaced, and unrelated/old frames are retained;
use a fresh directory for each run. `--check-only` writes no frames and samples
at most three snapshots, so it is not a full-file validation.

#### Display settings

The renderer uses a continuous `vik` colour palette. The defaults are
equivalent to:

```bash
--plot-upsample 4 --plot-interpolation linear --dpi 300
```

Upsampling adds intermediate display values to reduce visible blocks between
plane samples. `--plot-upsample` must be a positive integer: use `1` for the
native grid or `8` for a finer display. A 201 by 126 plane becomes 801 by 501
display points at factor four. Display point counts grow approximately with
the square of this factor, increasing memory and rendering cost per worker.
Increasing `--dpi` alone does not add intermediate grid values.

The default `linear` interpolation avoids overshoot. `cubic` gives smoother
gradients but can overshoot locally; GMT's clipping modifier constrains the
interpolated grid to the input's global range. Interpolation applies only to
the display after velocity is computed. Original samples, matrix orientation,
and reported extrema and saturation percentages are unchanged. A finer
display does not increase the simulation's physical resolution.

Solid red or blue regions can indicate saturation, which interpolation cannot
remove. The default `--velocity-limit 10` saturates beyond +/-10 cm/s; use, for
example, `--velocity-limit 50` to distinguish amplitudes within +/-50 cm/s.
Keep one fixed limit across frames for a comparable animation.

#### MPI placement and CPU affinity

MPI is optional. Rank/size are detected from Open MPI, PMI/PMIx, or Slurm
environment variables; no `mpi4py` dependency or inter-rank communication is
used. Each rank renders a disjoint frame range, without a shared marker file;
numbered PNG names preserve time order.

Check `mpirun --version` first: the following examples use **Open MPI** options;
Intel MPI and MPICH have different affinity controls. On Linux, child
processes inherit CPU affinity, so binding a Python rank to one core can
constrain all its joblib workers to that core. See the
[Linux affinity documentation](https://man7.org/linux/man-pages/man2/sched_setaffinity.2.html).

For **an existing allocation** of two nodes with at least four cores per rank:

```bash
mpirun -np 2 --map-by ppr:1:node:PE=4 --bind-to core --report-bindings \
  python3 script/animations_mpi.py \
  --plane-file outputfiles/planes/planedisplacements.0 \
  --nx 361 --ny 361 --dx 250 --dy 250 --dt 0.1 \
  --workers 4 --output-dir Animations_velocity_run2
```

`ppr:1:node:PE=4` gives each rank and its workers a shared group of four physical
cores under Open MPI's default core-based PE interpretation. Workers are not
individually pinned. For two such ranks per node, use `ppr:2:node:PE=4`, allocate
at least eight cores per node, and set `-np` to the intended total rank count.

For **two exclusively allocated nodes, each with at least 100 available
cores**, one rank per node can let 100 workers use the available CPU set:

```bash
mpirun -np 2 --map-by ppr:1:node --bind-to none --report-bindings \
  python3 script/animations_mpi.py \
  --plane-file outputfiles/planes/planedisplacements.0 \
  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 \
  --workers 100 --output-dir Animations_velocity_run3
```

This allows up to 200 rendering workers across the two nodes; each rank caps
its worker count at its assigned frame count. `--bind-to none` adds no Open MPI
binding; scheduler-enforced CPU sets still apply. For this layout,
set `-np` to the node count and match workers to each rank's available cores
and memory. Substitute the geometry and saved-snapshot interval for your file.
Mapping and binding control placement within an existing allocation; they do
not request nodes or CPUs. See the
[Open MPI mapping and binding reference](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/mpirun.1.html).

`--report-bindings` reports each launched rank's binding. On Linux, inspect the
CPU set visible to Python using the same placement policy:

```bash
mpirun -np 2 --map-by ppr:1:node --bind-to none --report-bindings \
  python3 -c 'import os, socket; cpus = sorted(os.sched_getaffinity(0)); print(socket.gethostname(), "allowed logical CPUs:", len(cpus), cpus, flush=True)'
```

Repeat with the `PE=4` policy when checking that alternative. Counts and IDs
refer to logical CPUs; simultaneous multithreading can make the count larger
than the physical core count. This checks the parent rank's allowed CPU set,
not individual worker placement or measured rendering performance.

#### GMT sessions and troubleshooting

Each worker keeps a read-only mapping, private temporary GMT directory, and
colour table for its assigned chunk. Its PyGMT session ends before the
directory is removed, including on rendering failure, and the corresponding
exit callback is unregistered. A reused worker starts a new session for its
next chunk. Both GMT and Python temporary files, including download locks,
use that private directory.

The script sets `GMT_DATA_UPDATE_INTERVAL = off` before importing PyGMT,
disabling remote catalogue refresh. These plots use local displacement files
and the installed `vik` palette, with no remote data requirement. See the
[GMT catalogue setting](https://docs.generic-mapping-tools.org/6.4/gmt.conf.html#gmt-data-update-interval).

Older versions could report `Not available in classic mode` or
`Module 'end' failed with status code 30` at exit after deleting the session
directory too early. Concurrent catalogue refreshes could also produce
`gmt_data_server.txt.download` or `Bad record counter` errors. Those exit or
catalogue diagnostics do not change already saved PNGs; use the current script
for session teardown and catalogue isolation. Rendering failures still
propagate as errors.

## Local investigation helpers

The following scripts and `script/tests/` are ignored and retained only in the
maintainer's local workspace. Their usage is recorded here for existing local
copies; a fresh GitHub checkout does not contain them. Keep these helpers and
their generated reports out of commits. They describe numerical diagnostics,
not physical validation criteria. New reports can live in `.local-diagnostics/`.

### compare_plane_outputs.py

Compare complete plane output directories with NumPy, for example legacy and
MPI-IO output for the seven-plane Izmir layout:

```bash
python3 script/compare_plane_outputs.py \
  outputfiles/planes outputfiles/planes_mpiio \
  --shapes 361x361 71x71 71x71 71x71 71x71 71x71 71x71 \
  --expected-snapshots 1000 --atol 0 --rtol 0
```

The two directories and `--shapes` are required. Supply one `NXxNY` shape per
plane ID, starting at 0, for matching `planedisplacements.N` files. Extra or
missing planes, incomplete frames, and inconsistent counts fail layout checks.
`--expected-snapshots` is optional but catches two equally truncated runs.

Defaults: `--atol 0` m, `--rtol 0` (dimensionless), `--chunk-points 262144`
XYZ points per read, and `--endian native` (`little`/`big` are alternatives for
both directories). A component passes when
`abs(candidate-reference) <= atol + rtol*abs(reference)`; NaN/Inf always fail.
Choose tolerances for the comparison's scientific purpose. Output reports byte
identity separately from numeric agreement, component errors, and zero-based
snapshot/point locations. It writes only to the terminal. Exit status is 0 for
PASS, 1 for numeric/nonfinite failure, and 2 for an input/read error.

### diagnose_plane_output.py

Locate signed peaks, classify their containing elements, and export histories:

```bash
python3 script/diagnose_plane_output.py \
  outputfiles/planes/planedisplacements.0 --nx 361 --ny 361 --dt 0.1 \
  --expected-snapshots 1000 --spacing 250 250 --xy-origin 0 0 \
  --domain 90000 90000 --cpml-thickness 439.453125 \
  --metadata outputfiles/planes/plane_point_metadata.*.csv --plane-id 0 \
  --points 2,93 2,92 1,67 \
  --output-dir .local-diagnostics/plane-report-run1
```

Requires NumPy and the local `compare_plane_outputs.py` beside it. Required
arguments are the binary file, `--nx`, `--ny`, saved interval `--dt`, and a **new**
`--output-dir`. The directory and its parents are created; existing output
directories are rejected. The completed binary is scanned twice, one snapshot
at a time, and file changes/nonfinite values are rejected.

| Optional control | Default / meaning |
| --- | --- |
| `--first-time` | 0 s |
| `--expected-snapshots` | Unset; enforce the expected complete-run count when supplied |
| `--top` | 5 distinct peak points per component/group |
| `--threshold` | 1 m; count points exceeding this diagnostic threshold |
| `--points` | Extra zero-based `I,J` histories; none by default |
| `--spacing DX DY` | Plane spacings, m; horizontal, unrotated planes only |
| `--xy-origin X Y` | Solver coordinates, m; 0 0 |
| `--domain LX LY` | Domain lengths, m, with lateral bounds 0..LX and 0..LY |
| `--cpml-thickness` | Geometric lateral-band width, m |
| `--metadata` | One or more rank CSVs from the same run/plane configuration |
| `--plane-id` | 0; selects rows from metadata |
| `--endian` | `native`, `little`, or `big`; default `native` |

Supply `--spacing`, `--domain`, and `--cpml-thickness` together or omit all three.
For a rotated/inclined plane, omit this horizontal geometry. Obtain the CSVs by
setting `HERCULES_PLANE_METADATA=1` before the simulation. Metadata gives actual
medium/CPML membership; distance-to-edge alone is only a geometric indicator.
Without metadata, material classification remains unknown.

Outputs are `summary.json`, `peak_points.csv`, `frame_peaks.csv`, and
`histories.csv`. Histories contain all snapshots for selected/peak points;
coordinates and signed XYZ displacements retain solver conventions. Exit status
is 0 on success and 2 on handled input/diagnostic errors.

### diagnose_acoustic_stations.py

Compare pressure, displacement, velocity, and acceleration at the same acoustic
stations across completed runs. Requires NumPy; `--plot` also needs Matplotlib.
Each `CASE_DIRECTORY` must contain `output/stations/station.N` and
`output/stations/station_mediums.txt`. Enable station velocity, acceleration,
and pressure output in the solver; this reader identifies columns by their
headers and checks acoustic ownership, density, time axes, and optional XYZ.

Create a JSON probe list such as `probes.json`, using actual station IDs and
coordinates from the runs. This example describes station 18; `cpml_by_case`
overrides membership when a changed layer includes the point:

```json
[
  {
    "station": 18,
    "label": "water_surface",
    "x_m": 500.0,
    "y_m": 23250.0,
    "z_m": 2812.5,
    "surface": true,
    "cpml": 0,
    "cpml_by_case": {"thicker": 1}
  }
]
```

```bash
python3 script/diagnose_acoustic_stations.py \
  --cases control=cases/control thicker=cases/thicker \
  --probes probes.json --expected-samples 200000 --late-start 80 \
  --gravity-cases control thicker --history-stride 200 --plot \
  --output-dir .local-diagnostics/station-report-run1
```

Required: `--cases NAME=CASE_DIRECTORY ...`, `--probes`, and `--output-dir`.
Case names must be unique letters/digits/underscores/hyphens. The first case is
the comparison reference. Probe entries require `station` and `label`; XYZ,
`surface`, `cpml`, and `cpml_by_case` are optional. CPML flags are 0, 1, or null
for unknown. Do not infer membership from a different mesh.

`--expected-samples` is optional. `--late-start` defaults to 80 s and must leave
at least three samples in the late window. `--gravity-cases` defaults to none;
list only cases that enabled gravity. The surface-relation residual is evaluated
only for declared surface probes outside CPML. `--history-stride N` optionally
exports every Nth raw sample plus the last; all metrics still use every sample.

The tool creates a new output directory and rejects an existing one. It writes
`summary.json`, `metrics.csv`, optional `station.N.npz` histories, and optional
`station.N.png` plots. Data are not filtered or detrended. Integration residuals
use saved samples and can reflect output sampling/rounding; agreement does not
independently validate acoustic reconstruction. Handled errors exit with status 1.

### check_acoustic_gravity_mode.py

Check flat-water dispersion and surface-gradient recovery against an analytical
mode using extracted production Hex8 stiffness/mass expressions:

```bash
mkdir -p .local-diagnostics/gravity-mode
python3 script/check_acoustic_gravity_mode.py \
  --depth 175.78125 --wavelengths 250 500 1000 2000 \
  --layers 8 16 32 64 --output .local-diagnostics/gravity-mode/report.json
```

Requires NumPy, the repository's `quake/forward/` source files, and a host C++
compiler (`CXX`, otherwise `c++`); CUDA and MPI are not required. Only `--output`
is required. Defaults are the values shown above: depth/wavelengths in metres,
and integer vertical layer counts. Every wavelength/layer combination is tested
with element size `depth/layers`.

The output parent must exist; an existing report is rejected. The JSON records
source hashes, compiler/NumPy details, frequencies, errors, and recovered-motion
metrics; a summary is also printed. Extracted C++ is compiled in a temporary
directory and removed afterward. Scope is uniform flat water, a periodic
horizontal mode, rigid bottom, and lossless acoustic physics **without CPML or
solid coupling**. Passing this check does not validate the coupled GPU solver.

### Local regression scripts in tests/

These ignored `unittest` modules are developer checks, not case launchers.
From the repository root, run all existing local modules or one selected file:

```bash
python3 -m unittest discover -s script/tests -p 'test_*.py' -v
python3 -m unittest discover -s script/tests -p 'test_compare_plane_outputs.py' -v
```

| Local file | Scope | Additional requirements |
| --- | --- | --- |
| `test_compare_plane_outputs.py` | Comparison tolerances, layouts, nonfinite/mutating files | NumPy, local comparison helper |
| `test_diagnose_plane_output.py` | Peaks, histories, geometry and metadata handling | NumPy, local plane helpers |
| `test_diagnose_acoustic_stations.py` | Station headers, ownership, metrics, changed CPML membership | NumPy, local station helper |
| `test_acoustic_gravity_mode.py` | Flat-water mode convergence and recovered displacement | NumPy, local gravity helper, C++ |
| `test_acoustic_plane_reconstruction.py` | Extracted reconstruction and CPML-history expressions | C++ |
| `test_acoustic_free_surface.py` | Shifted free-surface enrollment and assembly | C++, reconstruction test module |
| `test_acoustic_neighbors.py` | Remote-leaf geometry/material lookup helpers | C++17, reconstruction test module |
| `test_plane_point_metadata.py` | Metadata flags, coordinates, and containing material | C++ |
| `test_plane_output_mode.py` | Output-mode parsing, aliases, I/O-pool fallback | C++ |

C++ extraction checks need the repository sources and honor `CXX` (default
`c++`); they compile host-only temporary executables, not CUDA/MPI simulations.
Tests report pass/fail through `unittest` and use temporary fixtures rather than
production output. A clean checkout has no `script/tests/`; an empty discovery
run is not evidence that these checks passed.

## Structured DRM box output and flat binary stations

### read_drm.py

This standard-library reader handles two headerless native-float64 formats.
Choose **exactly one** of `--box LX LY DEPTH SPACING` or `--points N`, and supply
the positional filename and saved-snapshot interval `--dt`. `--first-time`
defaults to 0 s and can describe a restart segment.

For a structured DRM box:

```bash
python3 script/read_drm.py outputfiles/drmdisplacements.0 \
  --box 100 100 50 1 --dt 0.01
```

Use the exterior box dimensions and spacing in metres from `output_drm_boxes`.
Lengths must be integer multiples of spacing; LX and LY must each exceed twice
the spacing, and depth must exceed one spacing. Dimensions are local box axes;
centre and strike affect coordinates, not array shapes. Saved interval is
`output_drm_boxes_print_rate * simulation_delta_time_sec`.

Each box writes one binary with layout `(time, concatenated face sample, XYZ)`.
There are ten faces: exterior front/back/left/right/bottom, then the same five
interior faces. Shared edges/corners are repeated to keep grids rectangular.
The reader prints each face's `(time, rows, columns, 3)` shape and point offset.
The importable helper `box_layout(lx, ly, depth, spacing)` returns
`(name, rows, columns, point_offset)` tuples for slicing those faces.

Each box also writes `drmcoordinates.N` at setup: comment headers followed by
`X_m Y_m Z_m Vs_m_s Vp_m_s density_kg_m3 Qs Qp`. Coordinates are solver-local
X/Y/down metres, including rotation and building/topography shifts. The appended
columns are initial owning-element speeds (m/s), density (kg/m³), and dimensionless
material/BKT input Q values; `nan` means unavailable. They are not evolving
nonlinear or effective poroelastic attenuation. See the
[material semantics and mesh regeneration requirement](../doc/Extended%20README.md#structured-drm-box-output).
Coordinate row i matches binary sample i, including repeated edges/corners.
Older three-column-only readers need `usecols=(0, 1, 2)` for positions. The binary
reader interface is unchanged. Read and slice both arrays with the same layout:

```python
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path("script").resolve()))  # repository-root working directory
from read_drm import box_layout

layout = box_layout(100, 100, 50, 1)
npoints = sum(rows * columns for _, rows, columns, _ in layout)
metadata = np.loadtxt("outputfiles/drmcoordinates.0", comments="#", ndmin=2)
assert metadata.shape == (npoints, 8)
xyz = metadata[:, :3]
vs, vp, density, qs, qp = metadata[:, 3:].T
u = np.memmap("outputfiles/drmdisplacements.0", dtype=np.float64, mode="r")
assert u.size % (npoints * 3) == 0
u = u.reshape(-1, npoints, 3)
name, rows, columns, offset = layout[0]
face_xyz = xyz[offset:offset + rows * columns].reshape(rows, columns, 3)
face_u = u[:, offset:offset + rows * columns].reshape(-1, rows, columns, 3)
```

For the single-file station format enabled by `station_binary=yes`:

```bash
python3 script/read_drm.py outputfiles/station_displacements.0 \
  --points 2 --dt 0.01
```

Use `station_output_binary_file` as the filename, the number of
`output_stations` rows as `--points`, and
`output_stations_print_rate * simulation_delta_time_sec` as `--dt`. Layout is
`(time, station, XYZ)` in the station input order. This file contains displacement
only; use `read_stations.py` for per-station ASCII histories.

In both modes the reader requires at least one complete snapshot, scans all
components for finite values, and prints layout, time range, and maximum
absolute displacement for each global XYZ component in metres. It streams one
snapshot at a time and writes no files. Wrong geometry with the same point
count cannot be detected from the binary. Keep the run's input tables and read
completed output: DRM snapshots may be buffered by `drm_output_batch_size`,
whereas station binary output writes every saved snapshot. The solver overwrites
coordinate files on each run; preserve matching coordinates with each binary.
