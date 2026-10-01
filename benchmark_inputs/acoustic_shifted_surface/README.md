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
mpirun -np 1 "$PSOLVE_BIN" parameters.in
python3 check_output.py
```

The checker requires NumPy. It checks positive exposed surface and
gravity/CPML face counts, CPML enrollment of stations and plane points,
complete finite nonzero plane output, and station/plane agreement at an
acoustic CPML face probe with nonzero motion. Repeat in fresh copies with `HERCULES_PLANE_IO_MODE=legacy`
and with `acoustic_free_surface_gravity=no`. The pressure-release variant
must still enroll the shifted surface; its wavefield is intentionally different.
MPI-IO and legacy output from otherwise identical cases should agree within
roundoff. Try both supported GPU communication backends on allocated GPUs
when testing multi-rank communication. This fixture does not validate the
large Izmir case, a gravity-wave dispersion relation, or CPML reflection levels.

CPU-only regressions of the production surface assembly and output routines:

```bash
python3 -m unittest discover -s script/tests -p 'test_acoustic*.py'
```

The surface test checks an unshifted model, the 2812.5 m Izmir datum,
a building shift, buried water and a submerged face. The output test checks
a compact pulse's analytical integrated CPML gain and residual velocity.
