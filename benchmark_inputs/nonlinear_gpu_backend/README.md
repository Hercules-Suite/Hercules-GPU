# Nonlinear GPU backend acceptance check

Use a calibrated production-style nonlinear case in which at least 25% of the
local/global solid elements are nonlinear. Keep the MPI placement, GPU count,
input files, output rate, and launch environment identical between runs.

GPU is the production runtime default when `HERCULES_NONLINEAR_BACKEND` is
unset. The commands below set both backends explicitly so acceptance results
remain reproducible if defaults change.

This production CPU/GPU timing recipe applies only to material models that
support both production backends. `ManzariDafalias` production updates are
GPU-only. For that model, use the host/CUDA material-point suite and the
independent OpenSees comparisons in `../nonlinear_material_point/` for
constitutive validation, followed by a GPU production run; do not interpret an
unsupported CPU production launch as a performance reference.

CUDA nonlinear-station output gathers only the quadrature-point fields used by
each local station. It does not download the complete nonlinear history arrays
at every `output_stations_print_rate` step. Keep at least one nonlinear station
when checking CPU/GPU histories, and report the total and nonlinear element
counts so a reused single-layer mesh is not mistaken for a layered case.

Run the established reference backend first:

```bash
export HERCULES_NONLINEAR_BACKEND=cpu
/usr/bin/time -p mpirun -np <ranks> <psolve> parameters.in
```

Move the station outputs to a `cpu-output` directory, reset the case outputs,
and run the CUDA backend:

```bash
export HERCULES_NONLINEAR_BACKEND=gpu
export HERCULES_GPU_MPI_BACKEND=cuda_aware
export HERCULES_DISP_SYNC_MODE=boundary
/usr/bin/time -p mpirun -np <ranks> <psolve> parameters.in
```

Move those station outputs to `gpu-output`, then apply the parity and
end-to-end performance gates:

```bash
python3 compare_backends.py \
  --cpu-dir cpu-output \
  --gpu-dir gpu-output \
  --cpu-seconds <cpu-real-seconds> \
  --gpu-seconds <gpu-real-seconds> \
  --nonlinear-elements <count> \
  --total-elements <count>
```

The default gate requires station values to agree within `rtol=1e-6` or
`atol=1e-9`, a nonlinear fraction of at least 25%, and at least 2x end-to-end
wall-clock speedup. Record the machine, CUDA/MPI stack, material model, element
counts, rank/GPU placement, timings, and command line with the result. These
thresholds validate the execution backend; they do not validate the chosen
geostatic, tolerance, damping, substep, or material-property calibration.
