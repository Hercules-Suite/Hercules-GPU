# Plane velocity animations

[`animations_mpi.py`](../script/animations_mpi.py) generates vertical-velocity PNG
frames from completed `planedisplacements.N` files. Run it from the simulation
output directory, selecting one binary and specifying its dimensions and
snapshot interval manually. The script does not read `parameters.in` or any
other configuration file:

```bash
python /path/to/Hercules_GPU/script/animations_mpi.py \
  --plane-file planes_fullDomain/planedisplacements.0 \
  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 --check-only

# Open MPI: two allocated nodes, at least four cores per rank.
mpirun -np 2 --map-by ppr:1:node:PE=4 --bind-to core --report-bindings \
  python /path/to/Hercules_GPU/script/animations_mpi.py \
  --plane-file planes_fullDomain/planedisplacements.0 \
  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 \
  --workers 4 --output-dir Animations_velocity_run1
```

The plotting environment needs NumPy, xarray, joblib, and PyGMT/GMT. Worker count
is **per MPI rank**; match it to the CPUs allocated to each task. Use a fresh
output directory for each run. `--plane-file`, `--nx`, `--ny`, and `--dt` are
required. Select any plane by its binary filename, with that plane's dimensions.
`--dx` and `--dy` set the along-strike and down-dip spacing in metres (defaults:
500 and 325). Relative paths resolve from the current working directory.

The headerless native-float64 layout is `[time, strike point, dip point, XYZ]`
in C order, with down-dip points varying fastest. The helper assumes
strike=0 (modulo 360), dip=0; it cannot check these angles from the binary.
A plane with spacing/count pairs `500 201` and
`325 126` uses `(NT, 201, 126, 3)`, with solver X upward and solver Y rightward
in the plot. Neither transpose nor mirror is applied. Axes are distances from
the plane origin; geographic origin conversion and source-location consistency
remain separate checks.

`--dt` is the interval between **saved snapshots**, in seconds, equal to
`simulation_delta_time_sec * output_planes_print_rate` for the run. Supply
that product directly; the script applies no additional print-rate factor.
The example's `0.02` is appropriate only for snapshots saved every 0.02 s.
The script differentiates component 2, retaining downward-positive velocity,
then converts m/s to cm/s. It uses central differences except for one-sided
endpoint differences. Numerical velocities remain unclipped while the GMT color
table saturates at the chosen display limits. Reported extrema and saturation
percentages expose out-of-range amplitudes. Differentiating decimated snapshots
can differ from solver-step station velocity, especially at higher frequencies.

Each worker keeps one read-only mapping and an isolated temporary GMT directory
and CPT for its assigned frames. The worker closes its PyGMT session before
deleting that directory and removes the corresponding exit callback, including
when rendering fails. A reused worker starts a new session for its next chunk.
Both GMT and Python temporary files, including GMT download locks, use the
worker's private directory. Remote catalogue refresh is disabled with
`GMT_DATA_UPDATE_INTERVAL = off` before importing PyGMT; these plots use local
displacements and the installed `vik` colour table and require no remote data.
See the [GMT catalogue setting](https://docs.generic-mapping-tools.org/6.4/gmt.conf.html#gmt-data-update-interval).

Older script versions could print `Not available in classic mode` and
`Module 'end' failed with status code 30` at interpreter exit because the
session directory was already deleted. They could also report
`gmt_data_server.txt.download` or `Bad record counter` errors when concurrent
catalogue refreshes collided. Those diagnostics do not change already saved
PNGs; use the updated script to fix session teardown and catalogue isolation.
Rendering failures still propagate as errors rather than being suppressed.

No MPI communication or shared marker file is
required; disjoint frame ranges and numbered PNG names preserve time order.
`--check-only` checks layout and samples three velocity frames without loading
plotting libraries; it is not a full-file or physical validation. For restarted
output, supply `--first-time` for the first snapshot actually stored (default:
0 s), since the binary has no timestamps and the script cannot detect restarts.
Read native float64 files on a matching-byte-order machine.

## Open MPI placement and CPU affinity for animation workers

Check `mpirun --version` first: these examples use **Open MPI** options. Intel
MPI and MPICH have different affinity controls. On Linux, child processes
inherit CPU affinity through process creation and execution. Consequently,
binding a Python rank to one core can constrain all its joblib workers to the
same core. See the [Linux affinity documentation](https://man7.org/linux/man-pages/man2/sched_setaffinity.2.html).

For **two exclusively allocated nodes, each providing at least 100 available
cores**, start one Python rank per node and let its workers use the available
CPU set:

```bash
mpirun -np 2 \
  --map-by ppr:1:node --bind-to none --report-bindings \
  python /path/to/Hercules_GPU/script/animations_mpi.py \
  --plane-file planes_fullDomain/planedisplacements.0 \
  --nx 201 --ny 126 --dx 500 --dy 325 --dt 0.02 \
  --workers 100 --output-dir Animations_velocity_run1
```

`ppr:1:node` places one rank on each node; `--bind-to none` adds no Open MPI CPU
binding. Scheduler-enforced CPU sets still apply. Each rank starts 100 workers,
giving 200 rendering workers across two nodes. Workers are not individually
pinned. Set `-np` to the intended number of nodes for this one-rank-per-node
layout and size `--workers` to the cores actually available to each rank.

For an explicit core group, use
`--map-by ppr:1:node:PE=4 --bind-to core --report-bindings` together with
`--workers 4`, as in the smaller example above. The rank and its workers share
four physical cores under the default core-based PE interpretation. To put two
such ranks on each node, use `ppr:2:node:PE=4` and allocate at least eight cores
per node; the total `-np` must match the intended rank count. Mapping and binding
options control placement within an existing allocation; they do not request
nodes or CPUs from the scheduler. See the [Open MPI mapping and binding reference](https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man1/mpirun.1.html).

`--report-bindings` reports each launched rank's binding. For a Linux check of
the CPU set actually visible to Python, run the same placement policy with:

```bash
mpirun -np 2 --map-by ppr:1:node --bind-to none --report-bindings \
  python -c 'import os, socket; cpus = sorted(os.sched_getaffinity(0)); print(socket.gethostname(), "allowed logical CPUs:", len(cpus), cpus, flush=True)'
```

The reported IDs/count are logical CPUs; simultaneous multithreading can make
this count larger than the physical core count. Repeat with the `PE=4` policy
when checking that alternative. This checks the parent rank's allowed CPU set,
not individual worker placement or measured rendering performance.
