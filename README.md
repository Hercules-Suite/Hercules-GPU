# Hercules GPU

Hercules GPU is a parallel, octree-based finite-element code for large-scale
three-dimensional seismic-wave simulations. It combines adaptive mesh
generation, MPI domain decomposition, and CUDA-accelerated time integration for
regional wave propagation, nonlinear/plastic response, and coupled-physics
studies.

## Documentation

- [Hercules GPU Documentation](https://hercules-suite.github.io/Hercules-GPU/) — installation, input parameters, source models, running, output, and validation
- [Technical README](doc/Technical%20README.md) — concise build/run guidance and input/runtime variable index
- [Extended README](doc/Extended%20README.md) — comprehensive solver controls, formulations, examples, implementation details, and historical validation notes
- [Input and output utilities](script/README.md) — model conversion, topography preparation, validation, output readers, and plane animations
- [`benchmark_inputs/`](benchmark_inputs/) — reproducible input cases and model-specific examples

## Citing Hercules GPU

If you use Hercules GPU in your research, please cite this published workshop paper.

Wenyang Zhang and Yinzhi Wang (2026). [Modernizing Hercules for GPU-Accelerated
Earthquake Ground Motion Simulation on NVIDIA GPUs](https://doi.ieeecomputersociety.org/10.1109/SCW72768.2026.00233).
*SC26 Workshops*, High Performance Computing for Environmental and Earth
Sciences (HPC4EES). DOI `10.1109/SCW72768.2026.00233`.

## Capabilities

- Adaptive octree mesh generation and partitioning
- Distributed-memory execution with MPI
- CUDA-accelerated finite-element time stepping
- Elastic, viscoelastic, acoustic–elastic, and poroelastic simulations
- Nonlinear solid analysis with elastoplastic/plastic and hysteretic constitutive models
- Multiple point, force, finite-fault, plane-wave, and compliant-base sources
- Stations, output planes, checkpoints, and mesh diagnostics
- Three-dimensional material models and topography
- Dashpot and convolutional perfectly matched layer (C-PML) boundaries

## Requirements

A typical GPU build requires:

- C and C++ compilers
- CUDA toolkit and a supported NVIDIA GPU
- MPI implementation
- GNU Scientific Library (GSL)
- PROJ
- zlib and bzip2
- GNU Make

## Build

From the repository root, copy [`user.mk.example`](user.mk.example) if you do
not already have a `user.mk`. Set the MPI compiler wrappers, MPI and CUDA
paths, target GPU architecture, and GSL and PROJ include and library paths for
your system. Build with a job count suited to the build node.

```bash
cp user.mk.example user.mk
${EDITOR:-vi} user.mk
make -j 8
```

The primary executables are generated under `quake/forward/`, including
`qmesh` for mesh generation and `psolve` for the simulation phase.

## Configure a simulation

Current Hercules runs use one `parameters.in` file. It contains mesh, material,
analysis, damping, topography, source-directory, and output configuration.
Source-specific parameters remain in `source.in` and its associated auxiliary
files.

Begin with a current case under [`benchmark_inputs/`](benchmark_inputs/) and
consult the [input-file documentation](https://hercules-suite.github.io/Hercules-GPU/docs/input-files/)
before adapting it to a production model.

## Run

`run_phaseA_phaseB.sh` launches mesh generation and solving. Run it from the
case directory so relative paths in `parameters.in` resolve correctly:

```bash
export HERCULES_ROOT=/path/to/Hercules_GPU
cd /path/to/case

export PSOLVE_BIN="$HERCULES_ROOT/quake/forward/psolve"
export BASE_IN="$PWD/parameters.in"
export PHASE_A_NP=64
export PHASE_B_NP=8
export PHASE_B_OMP_NUM_THREADS=8

"$HERCULES_ROOT/run_phaseA_phaseB.sh"
```

`BASE_IN` is the complete path to `parameters.in`, not the case directory.
Adjust MPI rank and OpenMP thread counts for the target allocation and GPU
layout.

## Repository layout

| Path | Purpose |
| --- | --- |
| `octor/` | Adaptive octree mesh generation, partitioning, and load balancing |
| `quake/` | Material database tools, mesh generation, and the wave solver |
| `benchmark_inputs/` | Reproducible validation and feature cases |
| `script/` | Input preparation, validation, output readers, and animation scripts |
| `doc/` | Extended technical notes and specialized workflows |
| [`docs-site/`](docs-site/README.md) | Documentation website source, tests, and publishing instructions |
| [`user.mk.example`](user.mk.example) | Starter build configuration for a Linux CUDA/MPI system |
| `run_phaseA_phaseB.sh` | Unified Phase A and Phase B launcher |

## License

Hercules GPU is distributed under the GNU General Public License, version 3
(`GPL-3.0-only`). See [`LICENSE.txt`](LICENSE.txt). Existing copyright and
attribution notices are retained. Third-party components remain subject to
their respective licenses.
