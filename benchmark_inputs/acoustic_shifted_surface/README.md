# Shifted acoustic free-surface smoke test

This synthetic 1 km cube has a 125 m water layer at solver Z=250 m,
with air above and elastic bedrock below. A small Ricker double-couple in
the solid excites the fluid. These material/source values are test inputs,
not an Izmir calibration or an analytical wavefield reference.

Use a fresh copy for each run and a rebuilt GPU executable. From that copy:

```bash
export HERCULES_PLANE_METADATA=1
export HERCULES_PLANE_IO_MODE=mpiio
mkdir -p output/planes output/stations output/srctmp
mpirun -np 4 "$PSOLVE_BIN" parameters.in
python3 check_output.py
```

The checker requires NumPy. It checks positive exposed surface and
gravity/CPML face counts, CPML enrollment of stations and plane points,
complete finite nonzero plane output, and station/plane agreement at an
acoustic CPML face probe with nonzero motion. It also checks the 1024 equal-area
interface fragments and their coordinates in `output/interface.txt`; checking
total area alone misses incorrect subdivision at MPI partition boundaries.
Repeat in fresh copies with `HERCULES_PLANE_IO_MODE=legacy`
and with `acoustic_free_surface_gravity=no`. The pressure-release variant
must still enroll the shifted surface; its wavefield is intentionally different.
MPI-IO and legacy output from otherwise identical cases should agree within
roundoff. Repeat with one MPI rank and compare against four ranks: this
adaptive mesh exercises hanging-node and MPI sharing before acoustic
interface forcing consumes the new solid displacements. Keep the default
acoustic refinement factor; a uniform mesh can hide this ordering defect.
Try both `HERCULES_GPU_MPI_BACKEND=cuda_aware` and `host` on allocated GPUs.
For two completed case copies, compare their planes with:

```bash
python3 "$REPO/script/compare_plane_outputs.py" \
  case-4-ranks/output/planes case-1-rank/output/planes \
  --shapes 11x11 1x1 --expected-snapshots 800 --atol 1e-10 --rtol 1e-5
```

This fixture does not validate the
large Izmir case, a gravity-wave dispersion relation, or CPML reflection levels.
The checker also rejects displacement above 1 m for this deliberately weak
synthetic source. That generous bound catches gross exponential growth; it is
not an amplitude acceptance threshold for other cases. Backend agreement and
finite numbers alone cannot establish stability.

CPU-only regressions of the production surface assembly and output routines:

```bash
python3 -m unittest discover -s script/tests -p 'test_acoustic*.py'
```

The surface test checks an unshifted model, the 2812.5 m Izmir datum,
a building shift, buried water and a submerged face. The output test checks
a compact pulse's analytical integrated CPML gain and residual velocity.
The neighbor test checks final leaf bounds, gaps in the spatial cache,
corrected material classification, and absent remote octants.

## Recorded Horizon validation

Solver commit `8b04824`, tested on October 1, 2026 with NVIDIA 26.5,
CUDA 13.2 and Open MPI 5.0.11. Job `40545` completed all six 8 s runs:
four ranks with MPI-IO, legacy plane output, host communication, full host
displacement/pressure copies, and gravity disabled; plus one rank with
MPI-IO. All runs used `acoustic_formulation=velocity_potential`.
The gravity-disabled run tests a pressure-release boundary, not the legacy
pressure-variable formulation.

- Mesh: 27,392 elements and 2,400 hanging nodes, with different acoustic and
  solid resolutions. All variants have 1,024 correctly located interface
  fragments and 1,024 exposed acoustic top faces. Gravity-enabled runs
  enroll 448 acoustic CPML surface elements.
- Gravity-enabled plane peaks: approximately
  `(1.45948e-8, 1.72186e-8, 8.64858e-6)` m in X/Y/Z.
- Maximum one-rank/four-rank plane difference: `1.362e-17` m.
  The largest difference across the tested I/O and communication modes,
  including full host copies, was `1.189e-14` m. Each comparison passed
  `atol=1e-10 m, rtol=1e-5`, and each run passed station/plane parity.
- Job `40548` extended the four-rank MPI-IO case to 100 s, retaining
  `dt=0.001 s` and saving every 125 steps. Fields remained bounded, and the
  final 10 s vertical plane peak was `2.088e-9` m. Station/plane differences
  stayed below `4.4e-15` m. Its coarser output sampling gives a recorded
  vertical peak of `8.358e-6` m.

This fixture previously exposed exponential growth above `1e225` m when
the acoustic step read displacements before MPI/hanging-node synchronization.
It also catches subdivision based on remote aggregate octants: that path
produced 988 rather than 1,024 fragments on four ranks. The boundedness,
geometry, and rank-count checks address separate defects. These results do
not validate the full Izmir model or replace an analytical gravity-wave or
CPML-reflection benchmark.
