# BrickUP strong drained pressure boundary — validation, 2026-11-17

This change adds opt-in `poroelastic_surface_pressure_enforcement=strong` to
`poroelastic_formulation=brick_up` with `poroelastic_surface_pressure_bc=drained`.
Omitted enforcement (`natural`) and explicit `weak` retain the previous
compliance projection. BrickUP does not acquire full Biot's natural fluid
traction, relative-fluid inertia or consistent mass option.

Strong enforcement fixes independent nodes with positive Q1 dry-trace support.
On conforming flat/FEM faces this is ordinary nodal zero excess pressure. Interior
VT/groundwater cuts can fix an entire brick's active pressure support and suppress
its normal pressure gradient. **Cut support remains experimental, and strong is
not uniformly more accurate than weak.** See the
[implementation and input description](../../doc/Extended%20README.md#strong-brickup-drained-surface).

## Build and provenance

The boundary implementation was built on Horizon using `ml reset`, in an isolated
clean build. Its executable SHA256 was
`1db5444e4942636621139ba7df27c5635a1e1133faf64554ca8f0da6bf639cfb`.
The source snapshot included 107 files. The initial build connection was
interrupted; the same source and object tree was resumed successfully. Library
dependencies were checked. A subsequent build corrects the restart history
initialization exposed by this investigation; the boundary operator is unchanged.
The corrected executable SHA256 is
`1777c3925b8184fd9ecbfd5d6dd6a5cee116f1a734998bcebbdb3062553dc893`.

All listed jobs used `cuda-615-SC`. Unless stated otherwise, they used one node
and four MPI ranks/GPUs. Investigation fixtures, scripts, executable/source hashes,
input hashes, logs, station records and analytical figures are local-only.
No production material calibration is implied by these fixtures.

| Jobs | Coverage |
| --- | --- |
| 97508 | Surface saturation, aligned/buried/varying groundwater, explicit weak/omitted controls, zero penalty, sealed normal |
| 97509 | FEM/VT dry and coastal geometry; acoustic pressure and velocity-potential open-pore coupling |
| 97510 | Before/after elastic, viscoelastic, acoustic, acoustic–elastic and full-Biot controls |
| 97511, 97519 | Nonlinear elastoplastic, Rayleigh, non-CPML, active-permeability controls and parser guards |
| 97518 | Independent aligned-surface diffusion comparison, 4/2/1 m grids |
| 97520 | Five strong cases repeated with two nodes/eight ranks, host MPI |
| 97521 | Same five cases with one rank, CUDA-aware transport selected |
| 97523, 97524 | Groundwater and VT cut-surface diffusion comparisons, 4/2/1 m grids |
| 97525, 97534 | Restart diagnostics and incompatible-boundary rejection |
| 97540 | Corrected strong/weak restarts, fresh-run and full-Biot controls |
| 97542, 97543 | Elastic checkpoint/restart control (pressure output disabled) |

A Slurm failure status on a guard job is expected: the wrapper preserves the
nonzero solver exit. Acceptance checks require the intended diagnostic, not just
any failure. No production job reservation was changed. One initial elastic-only fixture
incorrectly requested pressure output and was rejected during parsing; job 97542
corrects that input. It is retained as a fixture error, not a successful test.

## Independent scalar drainage and fluid balance

An extracted host execution of the production diffusion and lumped pressure
update kernels uses zero solid motion and the one-dimensional problem

```math
\partial_t p=D\partial_{zz}p,\qquad
p(0,t)=0,\quad \partial_zp(1,t)=0,\quad
p(z,0)=\sin(\pi z/2),\quad D=1\ \mathrm{m^2/s}.
```

The independent exact solution is
$p(z,t)=\sin(\pi z/2)\exp(-\pi^2t/4)$, evaluated at $t=0.1$ s.
The time step scales as approximately $0.2h^2/D$. Pressure error decreases by
factors 3.996, 3.999 and 4.000 on 8/16/32/64 elements. At 64 elements, relative
pressure error is **0.000248%** and integrated drainage error is **0.00413%**.
Integrated storage loss balances the computed boundary reaction to **2.45e-15**
relative or better. This checks the scalar operator and reaction accounting;
it is not a GPU coupled-consolidation test. The production solver does not yet
write a strong-boundary reaction/outflow budget.

## Full-solver analytical diffusion comparison

The GPU/MPI control uses a 128 m cube, zero initial pressure, isotropic pressure
diffusivity $D=100$ m²/s, drained top, and sealed sides/bottom. The synthetic
material has `rho=1` kg/m³, `Vp=10` m/s, `Vs=5` m/s, porosity 0.4,
`poroelastic_biot_alpha=0`, fluid bulk modulus 1 Pa, fluid density 1 kg/m³ and
viscosity 1 Pa s. With the finite-grain storage convention,
$M=[0.4-0.4/(100-4\times25/3)]^{-1}$ Pa and $k=D/M$ m².
Zero coupling deliberately isolates the pressure equation; these are not soil
properties and do not validate coupled consolidation.

A `biot_explosion` at solver $(64,64,16)$ m has zero stress moment and pressure
moment 1 Pa m³. Its discrete rate samples a normalized Gaussian centered at
0.12 s with standard deviation 0.025 s; its first sample is zero. Runs last
1 s at `dt=0.0025` s on uniform 4/2/1 m grids, with no CPML. The analytical
reference convolves the actual piecewise-linear rate with the finite-box heat
Green function: even images at sealed boundaries and odd images at the drained
top. There is no fitted scale, time shift or filter. Halving convolution
quadrature spacing from 125 to 62.5 microseconds changes the reference by
2.23e-7 relative L2. Actual station coordinates are taken from
`station_mediums.txt`, including the terrain offset for VT.

For the aligned top at $z=0$, the relative L2 pressure errors at $z=4$ m are:

| Uniform cell size | Strong | Weak |
| --- | ---: | ---: |
| 4 m | 25.919% | 25.648% |
| 2 m | 5.910% | 6.366% |
| 1 m | **1.523%** | **5.295%** |

Strong surface pressure is exactly zero. Weak peak surface pressure is
8.34e-7, 1.95e-6 and 3.88e-6 Pa respectively: this fixed per-update blend does
not guarantee monotone surface-error reduction. At 1 m spacing, strong errors
at the other receivers $z=8,12,24,32$ m are 2.316%, 4.558%, 2.269% and 0.905%.
The coarse receiver next to the point source has 86.2% error, showing substantial
source/mesh error. These are joint spatial/time-discretization results at fixed
`dt`, not an isolated boundary convergence order.

## Cut surfaces: measured limits

The cut controls place a horizontal groundwater plane or planar VT terrain at
solver $z=2.5$ m, retaining the source at $z=16$ m. The VT source input depth is
13.5 m, and its station input depths receive the 2.5 m terrain shift. The two
families therefore have different absolute receiver coordinates and must not be
compared as if they were identical receiver sets.

| Geometry and receiver | h | Strong pressure L2 error | Weak pressure L2 error |
| --- | ---: | ---: | ---: |
| Groundwater, z=6.5 m | 4 m | 11.402% | 73.035% |
| Groundwater, z=6.5 m | 2 m | 19.642% | 9.672% |
| Groundwater, z=6.5 m | 1 m | **5.118%** | **2.023%** |
| VT, z=9 m | 4 m | 40.541% | 40.828% |
| VT, z=9 m | 2 m | 8.171% | 8.750% |
| VT, z=9 m | 1 m | **3.099%** | **5.226%** |

The groundwater result demonstrates that strong can be worse than weak. Its
nonmonotone refinement also changes how the existing center-based saturation
represents the plane. The deepest 1 m VT receivers have 0.716% and 0.375% strong
errors, but the receiver only 1 m from the point source has 60.9% error in both
modes. Those near-source errors are not evidence of accurate pressure everywhere.

A separate 8 m dynamic mesh with a groundwater plane at 67 m fixes two node
planes (8,450 independent nodes) under strong enforcement. Receivers at 67,
67.5 and 71 m all report zero pressure, with nonzero pressure at 75 m. This is a
measured overconstraint of the cut brick, not a desirable drainage accuracy
result. Use an aligned mesh where possible and retain weak as a comparison.

## Coupled and regression coverage

The dynamic surface fixture has exactly zero strong surface pressure, versus
41.178 Pa peak with weak enforcement. Fourteen unchanged-path comparisons cover
legacy weak/default aliases, the inactive strong penalty factor, impermeable
normal behavior, and the listed elastic/acoustic/full-Biot controls. Pressure
records match exactly at output precision; largest displacement/velocity relative
L2 differences are 9.81e-13 and 3.59e-12. Ten rank/backend comparisons across
1/4/8 ranks have relative L2 differences below 3.57e-9 for displacement, velocity
and pressure. These are numerical regression comparisons, not independent
physical references for the coupled interfaces.

Active finite-permeability nonlinear and non-CPML strong runs each fix 1,089
independent nodes. Initial zero-permeability versions correctly fix none and are
not counted as active drainage evidence. FEM/VT acoustic controls exercise dry
traces next to wet contacts with both acoustic scalar formulations; wet contacts
retain their existing interface law. Buried-table tests exercise the porous–dry
elastic interface. Applicable viscoelastic and nonlinear paths are checked;
existing unsupported geostatic, CPML and checkpoint combinations remain excluded.

The short coupled runs establish finite output and transport consistency over
their tested intervals. They do not establish arbitrary long-time stability,
general curved-terrain accuracy, minimal cut constraints, or a calibrated
multi-element nonlinear consolidation response.

## Restart correction

The first round trips exposed a pre-existing initialization defect in both strong
and legacy weak BrickUP. The checkpoint displacement prefix restores pre-loop
host order, but the ordinary upload refreshes only the current GPU displacement;
the old GPU displacement was left at zero on restart. The pressure extension
also supplied current/old host order, whereas the loop expects old/current
before its first swap. Initial state output alone missed the ensuing divergence.

The correction seeds both GPU displacement levels once after reading the
checkpoint and puts host pressure into pre-swap order. The file format and
ordinary timestep remain unchanged. Strong and weak resumed pressure histories
match the uninterrupted records at output precision; displacement/velocity
relative L2 differences are below 3.82e-16 over the final 500 steps of each
8,000-step run. A fresh strong run is also unchanged to numerical roundoff.
The original failed restart comparisons remain in the investigation record.

Changing a weak checkpoint to strong with nonzero saved fixed-node pressure
is rejected explicitly. Users must retain the original geometry and hydraulic
configuration when resuming; the guard does not make arbitrary boundary changes
valid. CPML and full-Biot restart exclusions are not relaxed by this work.

The elastic-only round trip (97542/97543) also matches its uninterrupted motion
to output precision apart from roundoff: maximum displacement/velocity relative
L2 differences are 1.35e-18 / 1.21e-18. All 29 required unchanged-path, transport and
restart comparisons passed. The full audit retains 87 accepted executions, two
diagnosed pre-fix restart failures, four intended guard rejections and the one
corrected elastic-fixture input error.
