# Wolf et al. (2022): validation history and full-Biot follow-up

**Current status (2026-11-14):** The matching `biot_explosion` source is now
implemented and Figure 5/7 runs have been compared directly with the supplied
Wolf author histories. **The tested meshes do not reproduce the paper's
accuracy.** See the historical [author-reference results](#wolf-author-reference-comparisons-2026-11-14)
and [report](../../doc/Hercules_Wolf_Reference_Validation_2026-11-14.pdf), then
the [uniform-mesh boundary control](#uniform-figure-5-boundary-distance-control-2026-11-14)
and [inviscid diagnostics](#inviscid-propagation-and-source-diagnostics-2026-11-14).
Full Biot and CPML are available; material-dependent hydraulics and the layered
401-point source remain missing. Sections below are chronological records:
statements about missing inertia, sources or CPML describe the executable at
that stage, not the current parser. The earlier
[full-Biot follow-up](#full-biot-implementation-follow-up) records its distinct
independent point-force and interface comparisons.

The CPML comparisons below that used the v14.0 maximum stretch of 5 can
be reproduced by setting `cpml_kappa_max=5` explicitly. The current
default is 3. Each comparison records its own profile and layer thickness.

Initial assessment: completed on Horizon, 2026-11-13. The adapted homogeneous reduced-Biot
tests pass the stated accuracy targets at 4 m spacing. This is **limited
verification of the current equations**, not reproduction of the paper's
full-Biot results or validation of every poroelastic boundary/interface path.
Original solver source examined and used for the executable: `2fd4f81` on
`origin/main`. That initial assessment changed documentation only. A subsequent
[sealed-dashpot correction and measured follow-up](#sealed-dashpot-correction-follow-up)
changes the solver boundary coefficients; the original results below remain
historical measurements of the earlier executable.

The [PDF report](../../doc/Hercules_Poroelastic_Validation_Report_2026-11-13.pdf)
includes the numerical results, convergence and waveform figures, boundary
assessment, documentation corrections and run provenance.

Reference: Wolf et al., *An efficient ADER-DG local time stepping scheme for 3D
HPC simulation of seismic waves in poroelastic media*, Journal of Computational
Physics 455, 110886, [DOI](https://doi.org/10.1016/j.jcp.2021.110886).
The authors provide input descriptions, SeisSol output and reference solutions
in [Zenodo record 5236133](https://zenodo.org/records/5236133). The downloaded
`model_and_results.zip` matches published MD5
`3d208188296847d77fdf9cf57b627502`. The downloaded paper and third-party archive
are not redistributed here.

## What can be compared

| Paper target | Paper requirements | Current Hercules limitation and valid scope |
| --- | --- | --- |
| Figure 5, section 6.2: homogeneous full space | Viscous and inviscid Biot materials; explosive source in normal stresses **and pore pressure**; solid and relative fluid velocities; eight receivers | Current `brick_up` evolves displacement and pressure, with Darcy diffusion but no independent fluid inertia/velocity state. No equivalent combined explosive source is exposed. A force-source comparison against a **reduced-Biot** Green tensor is a distinct, supported test. |
| Figure 6, section 6.3: two half-spaces | Inviscid fluid, two hydraulic material sets, source below the interface, reflected and transmitted waves; Gar6More3D reference | Zero viscosity is rejected; porosity, Biot coefficient and directional permeability are global inputs. These inputs cannot represent the exact two-material full-Biot problem. Sealed reduced-Biot/elastic equivalence can check existing interfaces, but cannot reproduce this figure. |
| Figure 7, section 6.4: free surface | Inviscid full Biot; traction-free and `p=0`; receivers at 0.5 and 500 m depth | Current drained surface is a weak pressure projection; zero viscosity and propagating slow waves are unsupported. Finite-viscosity surface refinement tests require a matching reduced-Biot reference before being called independent validation. |
| Figure 8, section 6.5: LOHp geometry | 500 m layer, two hydraulic materials, viscous fluid, 401-point explosive line source, solid and relative fluid velocities | Figure 8 is geometry, with waveform/error results in Figures 9–11. Hydraulic heterogeneity, the combined source, and independent relative-fluid-velocity output are missing. The existing sealed layered elastic-twin benchmark tests a narrower problem. |

The paper uses Biot's low-frequency viscous law **with fluid inertia**. Calling
Hercules a low-frequency formulation does not make its equations identical.
In particular, setting permeability to zero is a sealed limit, not an inviscid
limit. Replacing zero viscosity by a small positive number does not recover
the paper's inviscid propagating slow wave in this implementation.

## Material conversion

Use the paper's grain density and porosity to form the bulk mixture density
`rho = (1-phi) rho_s + phi rho_f`. Then convert the **drained** Lamé moduli to
input velocities using that same density:

```math
K_d=\lambda_d+\frac{2}{3}\mu,\qquad
\alpha=1-\frac{K_d}{K_s},\qquad
\frac{1}{M}=\frac{\alpha-\phi}{K_s}+\frac{\phi}{K_f},
```

```math
V_{P,d}=\sqrt{\frac{\lambda_d+2\mu}{\rho}},\quad
V_S=\sqrt{\frac{\mu}{\rho}},\quad
V_{P,u}=\sqrt{\frac{\lambda_d+2\mu+\alpha^2M}{\rho}}.
```

| Quantity | Table 6 homogeneous | Table 7, `z>0` | Table 7, `z<0` |
| --- | ---: | ---: | ---: |
| Mixture density (kg/m3) | 1664 | 2208 | 2022.8 |
| `alpha` | 0.5226666667 | 0.5 | 0.1315789474 |
| `M` (Pa) | 6.019261637e9 | 1.142857143e10 | 1.659464926e10 |
| Drained `Vp` (m/s) | 3296.268286 | 4037.864265 | 2423.437076 |
| `Vs` (m/s) | 1961.161351 | 2331.262021 | 1399.172048 |
| Reduced undrained `Vp` (m/s) | 3442.901750 | 4195.037984 | 2452.565969 |

These are benchmark-specific conversions, not solver defaults. The last row
is not the inviscid full-Biot fast-wave velocity quoted in the paper. Hercules
already derives finite-grain storage from `alpha` and the drained bulk modulus;
the `M=K_f/phi` simplification applies only to its `alpha=1` limit.

## Outer-boundary limitation

**Historical boundary assessment at `2fd4f81`.** The sealed normal-impedance
mismatch described here is corrected by the follow-up below. The lack of a
poroelastic PML and Darcy diffusion absorber remains current.

Hercules currently has **no poroelastic PML**. CPML leaves porous candidates in
its bands as elastic cells. All investigation cases therefore use
`abc_type=dashpot`, but this is not a validated poroelastic absorbing boundary.

Static review of `solver_init()` and `compute_setboundary()` in
`quake/forward/psolve.cc` shows that ordinary solid Lysmer coefficients use
`edata->Vp`, `edata->Vs` and density. They are assembled before poroelastic
classification and before the undrained fast path adds `alpha^2 M` to the
volumetric stiffness. Thus their normal impedance uses the **drained** input
`Vp`, even when the physical sealed P-wave propagates at `Vp_u`. The Darcy
operator has no outer diffusion absorber; its unconstrained exterior has a
natural no-flow condition. Positive mechanical damping alone does not establish
nonreflecting poroelastic behavior.

For the homogeneous Table 6 sealed continuum, the ideal normal-incidence
impedance-mismatch estimate `abs((Vp_d-Vp_u)/(Vp_d+Vp_u))` is about 2.18%.
This is an analytical warning, **not a measured reflection coefficient** for
the finite-element boundary. The mismatch is much larger for the soft-soil
materials in the existing layered benchmark. In particular, a poroelastic
run and its elastic-undrained twin do not have identical outer dashpots.
Late reflected waveforms cannot establish an interior coupling error.

Full-space comparisons must end before any possible mechanical boundary return,
and finite-permeability runs also need diffusion-distance and larger-domain
checks. The estimate `D=M*kappa/eta` gives an uncoupled diffusion scale
`sqrt(D*T)`; it is not a rigorous error bound for the coupled finite domain.
Matched doubled-domain runs retain material, spacing, source, receiver offsets
and time step. Surface and layered runs remain diagnostics wherever outer
reflections have not been excluded. No absorbing-boundary code is changed by
this validation work.

## Prepared investigation

Investigation-specific generators, inputs, analytical references, plots and
launchers are local-only. Existing benchmark inputs are unchanged.
All cases use SI units and solver X north, Y east, Z positive downward. Source
and receiver coordinates, checksums, source revision, time steps and comparison
targets are recorded in the local manifest and each case's metadata.

1. Fifteen homogeneous cases: 64, 32 and 16 m element edges, each with an elastic
   undrained control, sealed fast path, sealed explicit pressure coupling,
   Table 6 viscosity/permeability, and a larger-permeability diagnostic.
   They retain Table 6 material constants, the eight Table 5 offsets and the
   16 Hz, 0.07 s Ricker pulse. A Cartesian force `(1e9, 4e8, -2e8)` N replaces
   the paper's combined explosive source. A 4096 m cube and a 0.8 s record
   exclude the first possible elastic boundary return; this is shorter than
   the paper's 6000 m cube and 1.5 s record. The 64 m mesh is intentionally an
   underresolved pilot, not an accuracy result. Actual mesh size and timestep
   diagnostics must be inspected after startup.
2. Four surface cases: sealed and drained surfaces at 32 and 16 m, a source
   at 500 m depth and receivers from the surface to 500 m depth. These use a
   positive-viscosity force source and support convergence/pressure-release
   diagnostics only. They have no independent half-space reference yet.
3. Six layered cases: coarse/fine copies of the existing case 3, each with
   fast-path poroelastic, explicit poroelastic and elastic-undrained twin.
   Twin velocities are recomputed from the exact specified moduli. The source
   remains a double couple in elastic bedrock. Shared-code equivalence does
   not establish independent accuracy for hydraulic drainage or the paper's
   two-poroelastic-material interface.

The original 16 Hz runs exposed strong spatial dispersion in the elastic
control as well as the porous cases. Twenty additional **compact** cases
isolate convergence economically: 2048 m cube, Table 5 offsets multiplied by
0.25, 8 Hz force pulse delayed by 0.14 s, 0.5 s record, and 32/16/8/4 m spacing.
These retain Table 6 material constants but are explicitly a different
source/geometry test, not a reproduction of Figure 5. Seven further controls
double domain dimensions at fixed spacing (two homogeneous and two surface
cases), halve the 8 m explicit/diffusive time steps, and compare output-only
pressure reconstruction with every-step updates at common output times.

The increased permeability `6e-9 m2` is a mathematical stress test of the
reduced equations, not a physically justified full-Biot approximation at
8 or 16 Hz. The paper value is `6e-13 m2`. Time steps obey both conservative
undrained-wave and explicit-diffusion estimates. A successful dry-wave CFL
message alone is insufficient for this assessment.

The independent homogeneous reference solves

```math
\rho\ddot{\mathbf u}=\mu\nabla^2\mathbf u+
(\lambda_d+\mu)\nabla\nabla\!\cdot\!\mathbf u
-\alpha\nabla p+\mathbf F(t)\delta(\mathbf x),\qquad
\frac{\dot p}{M}+\alpha\nabla\!\cdot\!\dot{\mathbf u}
-\frac{\kappa}{\eta}\nabla^2p=0.
```

Its longitudinal spatial-Fourier denominator, with Laplace frequency `s`
and wavenumber magnitude `q`, is

```math
D(q^2)=(\rho s^2+(\lambda_d+2\mu)q^2)
\left(\frac{s}{M}+\frac{\kappa}{\eta}q^2\right)
+\alpha^2s q^2.
```

Partial fractions yield two decaying radial Green functions; numerical inverse
Laplace transformation supplies displacement, velocity and pressure. This
reference does not reuse Hercules quadrature or time integration. Local checks
compare its sealed limit with the classical isotropic elastic Green tensor,
the small-permeability limit with the sealed solution, finite-permeability
mass balance away from the source, and doubled transform padding. These check
the reference implementation, not the GPU solver.

Compare raw amplitudes and recorded time axes: no fitted scale, polarity
reversal, time shifting, or filtering to hide a discrepancy. Report vector
velocity/displacement relative L2 errors, pressure errors, and per-component
time-frequency envelope and phase misfits using ObsPy's Kristekova measures
(`w0=6`, 80 logarithmic frequencies, 1–40 Hz for 16 Hz sources and 1–20 Hz
for 8 Hz sources). These explicit
analysis choices are not asserted to match the authors' unpublished plotting
settings. The finest-mesh targets are 5% velocity-vector L2 and 10%
pressure L2 at every receiver. Failure requires refinement/investigation;
it must not be converted into a pass by reporting lag alone. Missing stations,
nonfinite samples, incomplete traces and solver failures invalidate a result.

The controls compare one versus multiple GPU ranks and host versus CUDA-aware
MPI, halve the time step, demonstrate mesh convergence, and compare output-only
with every-step diagnostic pressure updates. The sealed elastic control checks
the underlying elastic operator; it does not validate the explicit
pressure-coupling path by itself.

## Measured results

All 52 distinct inputs and six repeated rank/backend controls completed with
solver exit code zero: **58 successful runs** in ten jobs on reservation
`cuda-615-SC`. Completion alone is not an accuracy pass. All eight homogeneous
receivers were checked for finite, complete traces, correct coordinates, mesh
spacing and defined pressure before comparison.

For the compact 8 Hz test, maximum receiver **velocity-vector relative L2
error**, in percent, is:

| Element edge | Elastic control | Sealed fast path | Sealed explicit | Table 6 permeability | Darcy stress test |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32 m | 83.122 | 83.122 | 88.199 | 88.199 | 88.455 |
| 16 m | 22.666 | 22.666 | 23.165 | 23.165 | 23.200 |
| 8 m | 5.257 | 5.257 | 5.365 | 5.365 | 5.362 |
| 4 m | 1.340 | 1.340 | 1.365 | 1.365 | 1.331 |

The 8 m results remain outside the 5% velocity target; only the 4 m cases meet
both targets at every receiver. The following table gives maximum receiver
errors at 4 m. Each column can attain its maximum at a different receiver.

| Case | Displacement L2 (%) | Velocity L2 (%) | Pressure L2 (%) |
| --- | ---: | ---: | ---: |
| Elastic control | 0.831 | 1.340 | Not applicable |
| Sealed fast path | 0.831 | 1.340 | 0.418 |
| Sealed explicit | 0.848 | 1.365 | 0.437 |
| Table 6 permeability | 0.848 | 1.365 | 0.401 |
| Darcy stress test | 0.849 | 1.331 | 0.337 |

Across the finest-mesh velocity components, the largest envelope misfit is
0.471% and phase misfit is 0.386%; for pressure they are 0.300% and 0.080%.
These are supplementary measures, not replacements for the raw L2 targets.
The independent-reference self-checks gave relative differences of
`1.06e-14` against the elastic Green tensor, `3.54e-13` in the small-permeability
limit, `5.54e-6` for the finite-permeability mass-balance residual and `1.19e-8`
when doubling inversion padding.

The original, longer-offset 16 Hz tests at 64/32/16 m **do not pass**. At 16 m,
maximum receiver velocity errors remain 144.51–154.64%, including 144.51% in
the elastic control. Those failed results are retained. They demonstrate
inadequate accuracy for that setup, not an isolated poroelastic defect; the
compact convergence study is a separately identified test, not a relabeling
of the original results.

### Boundary, timestep and execution controls

The following changes use matched source samples, material and receiver
offsets. Aggregate relative L2 combines all recorded receivers for a field;
it is distinct from the maximum-receiver metric above.

- Doubling the compact 16 m domain changed displacement, velocity and pressure
  by less than `2.65e-8` aggregate relative L2 for both the sealed fast path and
  Darcy stress test. This supports boundary independence for those records;
  it does **not** measure or validate dashpot absorption.
- Doubling the surface domain changed all three fields by less than `4.31e-9`
  over 0–0.8 s. Over the full 1.2 s record, aggregate differences grew to at
  most 0.0265% displacement, 0.0184% velocity and 0.0414% pressure. No matching
  independent half-space solution was used, so these remain diagnostics.
- Halving the 8 m time step changed aggregate velocity by 0.1011% for sealed
  explicit coupling and 0.0956% for the Darcy stress test. Pressure changes
  were 0.000392% and 0.02155%, respectively. Spatial error dominates these
  tested cases; this is not a universal timestep bound.
- Output-only pressure reconstruction agreed with every-step updates at
  common output times within `1.22e-9` aggregate relative L2 across the three
  fields. The monitor confirmed 1,001 reconstruction steps and 2,999 skipped
  steps, so the optimized path was exercised.
- The homogeneous explicit and diffusive controls agreed within `1.52e-10`
  aggregate relative L2 for one versus four ranks and host versus CUDA-aware
  MPI. Comparisons use printed station precision, not full-field bitwise
  identity. The layered explicit case agreed within `4.21e-9` for one versus
  four ranks. Its host/CUDA-aware comparison retained pressure sensitivity:
  0.0116% aggregate pressure L2, up to 0.567% at an individual weak-pressure
  receiver, and a maximum absolute pressure difference of 13.54 Pa. This
  residual backend dependence is recorded rather than treated as exact parity.

The layered explicit/twin velocity difference decreases from 16.44% on the
coarse mesh to 4.46% on the fine mesh. Fast-path differences are 1.38% and
1.41%. These full-record differences include potentially different dashpot
reflections and share the elastic solver; they establish only limited
consistency, not independent interface accuracy. A boundary-controlled
half-space/layer reference is still needed.

### Provenance and retained artifacts

Successful jobs: `92788`, `92793`, `92794`, `92795`, `92796`, `92799`, `92800`,
`92812`, `92813`, `92822`. Job `92817` failed because MPI placed 16 GPU ranks
on a node exposing four GPUs; it was cancelled and its logs retained. Job
`92822` explicitly placed four ranks per node on four nodes and completed all
five finest-mesh cases. The failed placement attempt is not counted among the
58 successful runs.

The executable SHA-256 is
`a19809365c0e46c8dbff41000d138dc0b98bb91bc3a897bf2bcd59403f793b80`.
Jobs used `ml reset` and the executable's required GSL/PROJ library paths;
module listings, linked libraries, executable/source provenance and input
hashes are retained. Inputs and outputs are archived privately.
The local records contain per-job `analysis.json`, receiver comparison plots,
`results/controls.json`, `results/convergence.csv`, `results/convergence.png`,
`results/convergence.pdf` and `results/scheduler-status.txt`. Generated inputs,
diagnostic scripts and downloaded reference data are not included in Git.

## Remaining scope of the original assessment

Acoustic-only, acoustic–elastic and acoustic–poroelastic sealed/open-pore
interfaces, viscoelastic attenuation and nonlinear elastoplastic/plastic and
hysteretic paths are outside these linear lossless paper setups. No new
physical feature is implemented here, and these tests confer no validation
on those additional paths. Drained-surface independent accuracy, hydraulic
heterogeneity, full-Biot fluid motion and the exact Figures 5–11 reproductions
remain separate gaps.

## Sealed-dashpot correction follow-up

Completed on Horizon, 2026-11-13. This follow-up changes `solver_init()` in
`quake/forward/psolve.cc`; it does not add full-Biot physics or a poroelastic PML.
The original 58-run assessment above remains an independent historical record.
The PDF now begins with this correction and its results, followed by the
original nine-page assessment as a clearly identified historical appendix.

### Correction and supported scope

Active sealed poroelastic boundary cells use

```math
Z_n=\rho\sqrt{V_{P,d}^2+\alpha^2M/\rho},\qquad Z_t=\rho V_S.
```

Classification now precedes mechanical dashpot assembly. The element-local
storage helper is shared with the pressure/stiffness operators. The correction
applies to both fast and explicit pressure coupling, including eligible
attenuating and nonlinear paths with their **initial** moduli. No tangent or
frequency-dependent impedance update is introduced. It adds no per-step
kernel, allocation, transfer or field exchange. A setup-only MPI reduction
reports corrected and legacy boundary-element counts. Existing free-top,
side/bottom, edge/corner and CPML exclusions are preserved; shear and scalar
acoustic radiation coefficients are unchanged. No input parameter is added.

Eligibility requires all three permeabilities exactly zero, an `undrained`
pressure surface, and no active scalar-acoustic `open_pore` coupling. Models
with finite permeability, a globally `drained` surface or active `open_pore`
coupling conservatively retain the drained approximation, even for sealed
cells elsewhere. Their absorbing accuracy remains unestablished. The natural
unconstrained exterior pressure condition is no-flow, with no Darcy absorber.

### Reflection experiment and metric

Each before/after executable uses identical inputs. For each executable,
subtract a larger-depth reference at identical material, spacing, source,
receivers, timestep and horizontal dimensions. No time shift, filtering or
amplitude fit is used. This is a finite-domain difference, not an independent
half-space analytical solution or a universal reflection coefficient.

Normal incidence uses Table 6 material and a horizontal vertical-force sheet
at depth 256 m, with a 16 Hz Ricker velocity history of amplitude 0.001 m/s and
delay 0.07 s. The source's `base_force_vs_mps` is deliberately set to `Vp_u`
as the traction scale, with `base_force_scale=2`; that field does not alter
material Vs. Receivers are at depth 320 m, at the horizontal center and ±64 m
along X. Small/reference depths are 512/1024 m. The predicted reflected P peak
is about 0.200 s. The 0.14–0.26 s residual norm is divided by the larger-domain
incident norm over 0.03–0.14 s, separately for vector velocity and pressure.
The common timestep is `h/128000` seconds, with diffusion stability reducing
it where necessary. The record ends at 0.32 s.

The initial 1024 m horizontal aperture showed a residual floor: corrected
fast-path center velocity residuals were 2.316%, 0.866%, 0.729% at 16/8/4 m.
Doubling **both** horizontal dimensions to 2048 m while centering source and
receivers reduces aperture contamination. The table reports that wider
experiment, including the coarse result that worsens after correction.

| h (m) | Path | Velocity before (%) | Velocity after (%) | Pressure before (%) | Pressure after (%) |
| --- | --- | --- | --- | --- | --- |
| 16 | fast | 1.332 | 2.258 | 1.305 | 2.184 |
| 16 | explicit | 1.346 | 2.150 | 1.325 | 2.081 |
| 8 | fast | 1.772 | 0.545 | 1.775 | 0.542 |
| 8 | explicit | 1.804 | 0.512 | 1.807 | 0.509 |
| 4 | fast | 2.072 | 0.154 | 2.072 | 0.154 |
| 4 | explicit | 2.081 | 0.147 | 2.081 | 0.147 |

The old continuum impedance mismatch is 2.18%. The wider fine-grid results
approach that mismatch before correction and decrease after correction.
At 16 m the old coefficient partly cancels numerical boundary error, so a
single coarse run would give the wrong conclusion about the physical fix.
Residuals remain nonzero; no claim of perfect discrete absorption is made.

### Oblique incidence and finite permeability

The oblique test uses a vertical 1e9 N point force at (1024,1024,768) m in a
2048 m square domain. Small/reference depths are 1024/2048 m. Receivers at
depth 640 m have Y offsets `640*tan(theta)` for 0°, 30°, 45° reflected-ray
angles. These are spherical point-source fields, **not exact oblique plane
waves**. The residual is measured over 0.20–0.43 s and normalized by the
larger-domain **whole-record** norm over 0–0.45 s. Its denominator differs
from the normal-incidence table; the percentages cannot be directly compared
as plane-wave reflection coefficients. Both converted and other returned
wave contributions can enter this residual.

| h (m) | Path | Angle (°) | Velocity before (%) | Velocity after (%) | Pressure before (%) | Pressure after (%) |
| --- | --- | --- | --- | --- | --- | --- |
| 16 | fast | 0 | 0.343 | 0.501 | 0.376 | 0.665 |
| 16 | fast | 30 | 0.852 | 1.041 | 4.425 | 5.695 |
| 16 | fast | 45 | 1.827 | 1.840 | 15.040 | 14.887 |
| 16 | explicit | 0 | 0.335 | 0.489 | 0.355 | 0.636 |
| 16 | explicit | 30 | 0.828 | 1.010 | 4.454 | 5.695 |
| 16 | explicit | 45 | 1.773 | 1.756 | 14.504 | 14.074 |
| 8 | fast | 0 | 0.451 | 0.436 | 0.338 | 0.314 |
| 8 | fast | 30 | 1.349 | 1.543 | 3.829 | 4.963 |
| 8 | fast | 45 | 1.619 | 1.702 | 22.941 | 23.762 |
| 8 | explicit | 0 | 0.456 | 0.436 | 0.344 | 0.308 |
| 8 | explicit | 30 | 1.347 | 1.542 | 3.858 | 4.974 |
| 8 | explicit | 45 | 1.618 | 1.699 | 23.033 | 23.842 |

Oblique residuals persist and some increase after correction. A scalar
normal-impedance correction is not an all-angle absorber. These results do
not justify extending the sealed coefficient to finite permeability.

The finite-permeability diagnostic retains the original 1024 m aperture and
normal-incidence metric. Before/after results are unchanged to numerical
precision; the table gives the common residual, including source-aperture
effects. The larger permeability is a reduced-equation mathematical stress
test, not a justified full-Biot approximation. Only this 16 Hz source and
these mesh sizes were tested; frequency-dependent absorption is unvalidated.

| h (m) | Permeability (m2) | Velocity residual (%) | Pressure residual (%) |
| --- | --- | --- | --- |
| 16 | 6e-13 | 1.349 | 1.270 |
| 16 | 6e-9 | 1.661 | 1.738 |
| 8 | 6e-13 | 1.856 | 1.821 |
| 8 | 6e-9 | 2.572 | 2.715 |

### Coupled and execution-path checks

- The fine layered poroelastic–elastic test's maximum full-record velocity
  difference from its elastic twin falls from **1.407% to 0.119%** for the fast
  path. Explicit pressure coupling changes from **4.459% to 4.413%**. Remaining
  discrepancies are recorded, not attributed solely to rounding or declared
  solved. This is shared-code consistency, not independent interface accuracy.
- Pure elastic, Rayleigh-viscoelastic, elastic-regime nonlinear, standalone
  displacement-acoustic, acoustic–elastic (pressure and velocity-potential),
  acoustic–poroelastic open-pore (both scalar formulations), drained-surface
  porous and porous-interior/elastic-CPML controls retain their before/after
  station motion within 4.786e-13 relative L2. Defined pressure columns
  agree at printed precision. Undefined pressure at elastic stations remains
  explicitly undefined and is excluded from pressure metrics only.
- Sealed acoustic–poroelastic cases exercise both scalar formulations; sealed
  porous Rayleigh and nonlinear cases receive the new impedance and show the
  expected boundary-dependent changes. In the pressure-formulation sealed
  acoustic, porous Rayleigh and porous nonlinear controls, early velocity
  differences through 0.15 s remain below 6.9e-14 relative L2; later boundary
  returns account for the changed records. All complete with finite nonzero
  motion and defined pressure. The nonlinear cases use synthetic high-yield
  `vonMises_ep` parameters and have zero recorded plastic hardening `Kappa`;
  they exercise the GPU path in its elastic regime. They do not validate
  yielding, permanent deformation, calibrated hysteresis or evolving impedance.
- Eight one/four-GPU and host/CUDA-aware repeats cover sealed fast, sealed
  explicit, sealed acoustic–poroelastic and porous nonlinear cases. Maximum
  station velocity discrepancy is 5.544e-11 and pressure discrepancy 6.427e-09
  relative L2, below the 1e-7 comparison threshold. This is a station check
  at printed precision, not full-field bitwise identity or a multinode test.
- An independent face-geometry check compiles the actual boundary assembler
  locally and verifies all 27 location flags × 8 nodes × 3 components at both
  drained/undrained speeds, with and without `HALFSPACE`. Tangential terms,
  corners, zero interior contribution and top-face exclusion pass.

No new topographic FEM/VT, building, acoustic gravity, BKT attenuation,
yielding/hysteretic, heterogeneous-permeability or multi-node validation is
claimed. Their applicable setup path inherits the initial-modulus correction;
nonlinear tangent changes and frequency-dependent constitutive impedance
remain outside this absorber. Scalar standalone waves were not separately
excited; scalar propagation and output were exercised through fluid–solid
interfaces, while the standalone test uses displacement acoustics. Full-Biot
inertia, poroelastic PML and a general hydraulic absorber remain unimplemented.

### Follow-up provenance

There are **142 completed, checked runs** from 67 valid distinct inputs plus
eight repeated execution controls. Fourteen additional attempts (seven inputs
with each executable) in job `92922` failed startup because scalar acoustic
settings were supplied with `include_acoustic=no`. Their inputs/logs are
retained and excluded; corrected `_v2` inputs completed in job `92930`.
All intended final cases pass finite/complete-output, coordinate, mesh-spacing
and pressure-definition checks. Successful execution is not an accuracy pass.

Jobs: `92914`, `92918`, `92919`, `92920`, `92921`, `92922`, `92927`, `92928`, `92929`, `92930`, `92931`, `92932`, `92936`, `92937`, `92938`, `92939`. Runs use reservation `cuda-615-SC`, `ml reset`, required GSL/PROJ
paths and the supported repository launcher. The before executable is the
original `2fd4f81` binary; the after executable was rebuilt from an isolated
copy with this `psolve.cc` patch. Production files were not overwritten.

Before SHA-256:
`a19809365c0e46c8dbff41000d138dc0b98bb91bc3a897bf2bcd59403f793b80`.

After SHA-256:
`1168bfdbe54969298fbbf1e7f82f06dbe068d23642774e2368ca181f6b400962`.

The local validation records include
inputs, hashes, source patch, station traces, analysis, figures and helper
usage notes.
These investigation helpers/fixtures are ignored and unavailable in a clean
checkout. Solver build succeeded on Horizon. Local boundary geometry checks,
documentation lint, TypeScript, production/rendered-HTML and GitHub Pages
export audits passed. Existing benchmark inputs are unchanged.




## Full-Biot implementation follow-up

Completed on Horizon on 2026-11-13. `poroelastic_formulation=full_biot`
adds independent volume-weighted relative fluid displacement, coupled solid/fluid
inertia and constant low-frequency Darcy drag. `brick_up` remains the default.
The [full implementation/report PDF](../../doc/Hercules_Full_Biot_Implementation_Report_2026-11-13.pdf)
contains equations, raw waveforms, convergence, the failures that led to fixes,
and the supported scope. See the [full-Biot operational reference](../../doc/Extended%20README.md#full-biot-u-w-p-formulation)
for parameter defaults, units, restrictions and station columns.

### Independent physics comparisons

The source-matched full-space reference is an independent continuum Biot Green
tensor in Laplace space, inverted with the sampled source. It uses no Hercules
finite-element or time-update operators. Its finite-difference PDE residuals
are below 1.46e-6 and its small-permeability sealed limit agrees within 1e-7.
The reference includes independent fluid inertia, constant Darcy drag and the
same saturated material conversion as the paper. These are adapted point-force
tests, not the paper's combined normal-stress/pore-pressure explosive source.

A 512 m cube has a centered force `(1e9,4e8,-2e8)` N, a 25 Hz Ricker pulse
delayed 0.04 s, and six receiver offsets `(40,0,0)`, `(60,0,0)`, `(0,40,0)`,
`(0,60,0)`, `(23,23,23)`, `(35,35,35)` m. The initial source sample is zero in
both solver and reference. Records end at 0.12 s, before the earliest outer
return at 0.12589 s. Material: saturated density 1664 kg/m3, porosity 0.4,
tortuosity 2, drained lambda/mu 5.28/6.4 GPa, alpha 0.5226666667,
M 6.019261637 GPa. Inviscid fast/slow P speeds are 3590.329/1055.976 m/s.

The table reports the **maximum receiver relative L2 error** for raw vector
solid velocity, pressure and vector Darcy flux over the full prescribed record.
Velocity timestamps account for the known half-step stagger; there is no fitted
time shift, amplitude scaling or output filter. Targets are 5% for velocity and
10% each for pressure and flux. Sealed flux must be exactly zero. A completed
run is not automatically an accuracy pass.

| Case | h (m) | v error (%) | p error (%) | q error (%) | Job / binary |
| --- | ---: | ---: | ---: | ---: | --- |
| inviscid | 8 | 41.453 | 100.687 | 60.320 | 93081 / v09 |
| inviscid | 4 | 8.663 | 100.393 | 61.236 | 93025 / v05 |
| inviscid | 2 | 2.382 | 44.444 | 27.313 | 93025 / v05 |
| inviscid | 1 | 2.198 | 10.952 | 6.651 | 93033 / v06 |
| inviscid | 0.5 | 3.470 | 3.167 | 3.338 | 93063 / v08 |
| viscous | 4 | 8.606 | 75.599 | 30.319 | 93025 / v05 |
| viscous | 2 | 2.330 | 33.598 | 13.430 | 93025 / v05 |
| viscous | 1 | 2.075 | 8.757 | 3.668 | 93033 / v06 |
| transition | 2 | 2.127 | 1.387 | 2.628 | 93036 / v06 |
| sealed | 2 | 2.633 | 1.623 | 0.000 | 93036 / v06 |

Inviscid means viscosity zero with positive permeability. Viscous uses
viscosity 0.001 Pa s and permeability 6e-9 m2; transition uses 6e-10 m2.
Sealed material has zero permeability and locked relative fluid motion.
The mesh sequence uses dt=h/40000 s. The finest inviscid, 1 m viscous, 2 m
transition and 2 m sealed cases meet the stated targets. Pressure and flux
converge strongly; maximum solid-velocity error is not monotone from 1 to 0.5 m,
so no asymptotic solid-velocity order is claimed. The finest cube has
1,073,741,824 cells on 32 GPUs across eight nodes.

At fixed h=2 m, halving dt from 50 to 25 microseconds (including source sampling)
gives v/p/q errors of **2.546/44.471/27.336%** for inviscid material and
**2.497/33.618/13.449%** for viscous material (job 93087, v09). The large coarse
pressure error persists under time refinement; spatial refinement is necessary.

### Normal-incidence interface accuracy

An independent inviscid modal-scattering reference enforces mechanical
continuity and either zero normal Darcy flux or open-pore pressure/total-flux
continuity. Unit-power balance closes within 2e-12. A 2048 by 2048 m force sheet
uses a 12.5 Hz Ricker pulse delayed 0.08 s; records end at 0.24 s before relevant
outer returns. Water spans depths 0–256 m above porous material to 1024 m;
water density/Vp are 1000 kg/m3 and 1500 m/s. The elastic test has porous
material to 768 m and elastic material to 1536 m, with density/Vp/Vs of
2400 kg/m3, 4000 m/s and 2300 m/s. Source depth is 512 m. Four receivers sample
incident, reflected and transmitted waves. Detailed offsets and force scaling
are retained in case JSON and the PDF.

| Case | h (m) | v error (%) | p error (%) | q error (%) | Job / binary |
| --- | ---: | ---: | ---: | ---: | --- |
| acoustic sealed | 8 | 3.403 | 70.632 | 41.918 | 93029 / v05 |
| acoustic sealed | 4 | 0.841 | 23.144 | 16.969 | 93030 / v05 |
| acoustic sealed | 2 | 0.209 | 5.699 | 4.446 | 93064 / v08 |
| acoustic open | 8 | 6.988 | 56.624 | 65.695 | 93029 / v05 |
| acoustic open | 4 | 1.131 | 23.980 | 21.632 | 93030 / v05 |
| acoustic open | 2 | 0.287 | 6.350 | 5.352 | 93064 / v08 |
| elastic sealed | 8 | 2.162 | 61.825 | 41.556 | 93031 / v05 |
| elastic sealed | 4 | 0.641 | 23.209 | 18.929 | 93031 / v05 |
| elastic sealed | 2 | 0.158 | 5.702 | 4.772 | 93065 / v08 |

Here h is the porous cell edge; water cells are one refinement level smaller.
All three 2 m cases meet the predefined targets. The fine acoustic geometry
has 1.476 billion cells on 32 GPUs; the elastic geometry has 805 million on
16 GPUs. Acoustic pressure and velocity-potential formulations at 4 m agree
within 5.70e-11 relative L2. Initial 25 Hz/coarse attempts remain recorded;
the table is the paired 12.5 Hz refinement sequence. No independent oblique
interface, free-surface or general heterogeneous hydraulic validation is claimed.

### Numerical and execution corrections

The equal-order collocated fluid/pressure discretization has a spurious grid
branch. Unstabilized inviscid h=2 m gave 174.4% pressure and 120.5% flux errors
even though solid motion was much closer. Full Biot therefore includes
fourth-order numerical pressure dissipation `-dt K_a C^-1 K_a p`, with
`poroelastic_pressure_stabilization=0.05` by default. It preserves constant
pressure and assembled total fluid content, dissipates pressure energy, scales
as h cubed in the inviscid limit and decreases with drag; impermeable directions
have zero coefficient. It is not physical diffusion or filtering. Raising the
coefficient to 0.1/0.2 at h=2 m only changes pressure error to 42.859/40.943%
and flux error to 26.341/25.211%; it does not substitute for spatial resolution.

The production algebra passes 200 independent generalized-eigenvalue and 1000
matrix-exponential/energy/momentum checks. A separately assembled nonuniform
operator passes 1000 energy/content tests, and a 2001-mode Fourier check passes
the combined wave/dissipation bound to numerical roundoff. This is not a proof
for every nonlinear, heterogeneous or boundary configuration.

Two concurrency bugs were exposed by the new fields and adaptive/multinode
checks:

- Different solid/fluid message tags reuse the same GPU receive buffers.
  Fences now identify the actual schedule/direction storage before MPI can
  overwrite a buffer still read by CUDA unpack. Before the fix, the fine
  inviscid v/p/q errors were 1381.9/442.5/1004.4% (93044), and fine elastic
  interface velocity error was 264.4% (93046). Corrected jobs 93063/93065 give
  the small errors in the tables. The failed candidates are retained and
  excluded from the accepted convergence sequence.
- A pre-existing host OpenMP assignment loop shared its dangling-node index.
  Making that index private to each iteration removes nondeterministic
  adaptive host/GPU differences as large as 1.88e-5. Turning assignment OpenMP
  off independently isolated the cause before the fix. Final adaptive and
  small-reproducer comparisons agree at approximately 1e-12.

The stabilization exchange also allocates its scalar buffers when batching is
active. Binary station setup now creates/reads the fluid companion directory.
Full Biot uses assembled solid force directly in the coupled update rather than
recovering force by subtracting large displacement values.

### Standalone, coupled and output checks

Final v09 host/CUDA-aware checks cover elastic, Rayleigh-viscoelastic,
elastic-regime von Mises, displacement-acoustic, acoustic–elastic,
reduced-Biot acoustic sealed/open (both scalar open formulations), porous
Rayleigh/nonlinear and drained cases against the pre-feature executable.
All 22 final comparisons agree within 3.03e-12 relative L2. Ten final full-Biot
transport comparisons, including a 56.6-million-cell adaptive mesh with 197,120
hanging nodes, agree within 4.69e-12. One/four GPU, scalar batching and GPU overlap
controls are separately retained. These are station comparisons, not bitwise
full-field identity or performance benchmarks.

Zero X permeability produces exactly zero wx/qx/vrx. Binary displacement agrees
with ASCII within its known six-decimal scientific formatting (relative L2
1.558e-7 before rounding); fluid companions agree within 3e-13 and include a
coordinate-frame sidecar. Defined columns are finite/complete; undefined elastic
pressure remains explicitly undefined and is excluded only from pressure metrics.
Invalid tortuosity, nonpositive mass matrix, checkpointing and porous CPML inputs
are rejected with the expected diagnostic.

Nonlinear tests have zero recorded plastic hardening Kappa: they exercise the
elastic regime, not yielding, permanent deformation or hysteresis. Rayleigh is
tested; BKT and each other attenuation/constitutive family remain separate gaps.
Drained upper-surface behavior has execution checks but no independent waveform
reference here. Building, acoustic-gravity and topographic cases are not validated.

### Remaining scope and provenance

There is still **no poroelastic PML**. Lateral/bottom fluid boundaries are sealed,
and mechanical dashpots are approximate. In this initial implementation,
topographic FEM/VT, buried/spatial
water-table cuts, checkpoint/restart, geostatic ramps, fixed-base and DRM histories
are unsupported and rejected. Hydraulic constants remain global. The paper's
exact multilayer/source configurations still need layer-dependent hydraulic tables
and a simultaneous normal-stress/pore-pressure source; an independent free-surface
reference is also needed. This follow-up does not claim literal reproduction of
all paper figures.

There are **161 completed runs with audited finite/complete outputs**,
including known failed-accuracy candidates, plus 17 non-guard failed attempts
and 8 expected guard rejections. Job 93047 additionally failed in harness
setup before its first case. These counts are not counts of physics passes.

Fine GPU physics uses v08 in jobs 93063–93065. v09 changes only the host
OpenMP assignment index; transport/regression/timestep jobs are 93081–93087.
Final v10 also uses the common acoustic material classification for remote
fluid-boundary neighbors. Accepted near-zero acoustic Vs (1e-7 m/s) matches
exact zero on one/four GPUs within 2.35e-11 relative L2 in jobs 93092–93093.
Earlier accepted coarse physics uses v05/v06. Each metric carries its binary/job;
all old failures are retained. The baseline is the completed sealed-dashpot
implementation, source commit `345e53db5f28b677cf21087487b2b58d2eb0b2a6`.

- Final v10 executable SHA-256: `3a54110ae89287fc648463ea1e812793b453fde1fe82ca2dd1224581f38b1770`.
- Fine-GPU v08 SHA-256: `12d639faca9e36ec03012e9b3a5a278e1dd15638509008a6be9ebe34a89b0a12`.

The isolated Horizon build uses `ml reset`, default NVIDIA 26.9/CUDA 13.3/
OpenMPI 5.0.11 and the required GSL/PROJ paths. Runs use `cuda-615-SC` and the
supported launcher. Immutable binaries, inputs, references, raw station traces,
source/input SHA manifests, job metadata, build logs, figures and PDF generator
are retained in ignored local diagnostics; they are unavailable in a clean
checkout. Existing benchmark inputs and production binaries were not overwritten.
Solver builds, independent checks, documentation lint/TypeScript, both production
and GitHub Pages rendered audits, and desktop/narrow-screen guide inspection pass.


## FEM topography follow-up (2026-11-13)

The maintainer requested **FEM** topography for full Biot. This follow-up extends
source baseline `f06bd68dac5450c84d1105ce17adea6e386a2c0c`; the earlier section and
PDF describe their original flat-only scope. The new
[PDF addendum](../../doc/Hercules_Full_Biot_FEM_Topography_Report_2026-11-13.pdf)
contains the implementation summary, waveforms and limits.

### Geometry and hydraulic changes

`include_topography=yes`, `computation_method=fem` now works with `full_biot`.
FEM retains whole axis-aligned Hex8 cells, so regular coupled mass,
pressure-force/continuity quadrature, stabilization, GPU/MPI stepping and
receiver interpolation remain applicable. No per-step operator was changed.

- Local/remote air and owner-confirmed carved voids are distinguished from
  material neighbors. Acoustic open pores remain an interface condition.
- `undrained` constrains relative-fluid motion along each exposed mesh-face
  normal. `drained` permits flow with natural zero excess fluid-pressure
  traction on **all** exposed faces, including vertical stair risers. Domain
  lateral/bottom truncations remain sealed. Startup reports face areas by axis.
- Full-Biot saturation follows physical terrain for both `full` and `sqd`
  material sampling. Retained surface cubes whose centers lie above terrain
  are fully saturated and use depth zero for the porous-depth limit. This
  matches FEM's discrete solid domain and avoids an unintended elastic cap.
- `vt` still aborts with a specific cut-cell-fluid-operator diagnostic.
  Buried/spatial groundwater grids, checkpoint/restart, geostatic ramps,
  fixed-base and DRM histories remain unsupported. No poroelastic PML was added.

The hill's 8 m drained run reports exposed X/Y/Z face areas
`20480 / 20480 / 262144 m²`; the vertical risers therefore participate in
hydraulic drainage. Faces that meet at stair corners constrain their respective
normal components in the sealed case. The normal is the **discrete face normal**,
not the normal of the continuous interpolated terrain.

### Independent free-surface comparison

An independent inviscid modal solution imposes zero total traction and either
zero normal Darcy flux (undrained) or zero pore pressure (drained). Its modal
reflection matrix satisfies unit energy-flux balance. No Hercules element,
time-stepping or boundary matrix enters the reference; no amplitude fitting or
time shifting is used. The material is the same full-Biot material documented
above: density 1664 kg/m³, fluid density 1040 kg/m³, porosity 0.4, tortuosity 2,
positive permeability `6e-9 m²`, zero viscosity, stabilization 0.05.

The FEM surface is at solver Z=64 m in a 2048×2048×1024 m domain. A vertical
solid-force sheet 128 m below it uses `2*rho*Vp_u*0.001 N/m²` times a 12 Hz
Ricker pulse delayed 0.08 s. Stations are 32/64/192 m below the surface. Histories
run to 0.32 s, with `dt=h/40000`. Velocities/flux are compared at `t-dt/2`,
pressure at `t`. The reported relative L2 maximum is over the three receivers.
Targets set before evaluation: 5% velocity and 10% pressure/flux.

| Surface | h (m) | Velocity % | Pressure % | Darcy flux % | Job |
|---|---:|---:|---:|---:|---:|
| undrained | 8 | 5.980 | 57.162 | 63.323 | 93180 |
| undrained | 4 | 1.814 | 30.193 | 34.014 | 93184 |
| undrained | 2 | 0.441 | 8.130 | 9.143 | 93193 |
| drained | 8 | 2.070 | 49.044 | 58.276 | 93180 |
| drained | 4 | 0.758 | 23.765 | 19.851 | 93185 |
| drained | 2 | 0.208 | 5.638 | 4.981 | 93194 |

The 2 m cases pass these chosen tolerances; coarse cases do not. The translated
8 m topography-off twins agree within `5.5e-11` relative L2 across solid motion,
pressure and relative-fluid output. This is a separate shared-code consistency
check, not the independent reference.

### Translated interfaces and combined coast

The earlier normal-incidence acoustic/poroelastic and elastic/poroelastic
geometries are shifted down by 64 m and run with FEM topography. Profiles and
source sheets are shifted in solver Z; station input depths remain surface
relative. Hydraulic properties, source pulses and references are unchanged.
These cases retain a flat outer surface and do not establish arbitrary coastal
scattering accuracy. Maximum raw receiver errors at 2 m porous spacing are:

| Contact | Velocity % | Pressure % | Darcy flux % | Job |
|---|---:|---:|---:|---:|
| Acoustic no-flow | 0.209 | 5.689 | 4.446 | 93212 |
| Acoustic open-pore | 0.287 | 6.350 | 5.352 | 93224 |
| Elastic no-flow | 0.158 | 5.702 | 4.772 | 93196 |

A separate 512 m coastal integration model combines a dry FEM hill, flat water
surface, sloping explicit seabed, porous solids and a deeper elastic region.
Jobs 93207–93208 execute acoustic `pressure` and `velocity_potential`, `no_flow`
and `open_pore`, and shared/field-specific constraints. Its mesh contains
nonmatching acoustic/solid traces and active pore-pressure hanging constraints.
Both scalar formulations and host/CUDA-aware MPI agree within `3.4e-13` across
applicable outputs. Shared versus field-specific maps are execution checks,
not an assertion that the two constraint discretizations are identical.
No independent coastal scattering reference was evaluated.

### Hill convergence and applicability limits

The Gaussian hill has elevation `64*exp(-r²/100²)` m, centered in a 512 m cube,
with `topographybase_zcoord=128 m`; the common topography input grid spacing is
16 m. A force `(1e8,4e7,-2e7) N` 96 m below the summit excites non-normal paths.
Five receivers are 16–64 m below the local surface. A 6 Hz Ricker pulse delayed
0.2 s runs to 0.4 s, with `dt=h/40000`. These are **finite-domain mesh differences**,
including the fixed approximate truncations, not errors against an exact hill
solution or a demonstration of absorber accuracy. Differences use the coarser
physical timestamps, with the half-step correction for velocity/flux.

| Surface | Mesh pair (m) | Velocity % | Pressure % | Darcy flux % |
|---|---|---:|---:|---:|
| undrained | 4 → 2 | 5.520 | 13.576 | 12.822 |
| undrained | 2 → 1 | 2.671 | 6.159 | 5.314 |
| drained | 4 → 2 | 4.797 | 12.879 | 6.903 |
| drained | 2 → 1 | 2.300 | 5.614 | 3.046 |

Jobs 93202–93204 cover the 6 Hz refinement series. Its 2→1 m differences are
below 2.7% velocity, 6.2% pressure and 5.4% flux. The earlier **25 Hz** hill
(jobs 93179/93183) remains under-resolved: 4→2 m differences reach
48.85% velocity, 81.11% pressure and 66.77% flux. Both series are retained;
the frequency change is not hidden or presented as refinement of one problem.
A new production model still needs its own geometry and slow-wave resolution
study. FEM's staircase is not an exact curved-surface boundary.

### Transport, other physics and regressions

- Adaptive hill: 138125 elements, 36744 hanging nodes, and 17904 kept
  pore-pressure constraints. Host/CUDA-aware and four/eight-rank runs agree
  within `2.5e-13` (93198–93200), for both hydraulic surface conditions.
- Binary station mode retains relative-fluid output (93179); fluid companion
  histories agree with ASCII mode within `2.0e-13`.
- Viscous and impermeable hill cases execute with finite fields (93179).
  The coarse impermeable elastic-undrained twin differs by 2.91% in velocity;
  this consistency check is not exact algebraic equivalence and does not
  validate permeable flow or drainage.
- The nonlinear effective-skeleton path matches the linear case in its elastic
  regime to `1.1e-12` across output quantities (93179). It does **not** validate
  yielding, permanent deformation or hysteresis. Material/geostatic/tolerance
  values remain case-specific. Rayleigh is exercised; BKT and every other
  attenuation/constitutive family remain separate gaps.
- Baseline/new-executable comparisons include elastic, Rayleigh, acoustic,
  acoustic–elastic, reduced poroelastic, acoustic–poroelastic, nonlinear
  elastic-regime, and full-Biot paths (93187–93188, 93197–93198). Applicable
  differences are within `2.4e-12`. Existing reduced-Biot FEM hill and
  acoustic–elastic coastal cases agree at the written precision (93210–93211).
- Source-extracted C++ checks pass 168 face-classification cases: all six normals,
  local/remote material and air, carved voids, refined neighbors, hydraulic
  contacts and domain truncations. Eighteen production classification checks
  cover full/sqd, retained surface cubes, finite depth, groundwater, air,
  reduced and flat controls. These local-only harnesses are not CUDA tests.

### Follow-up provenance

There are **82 completed runs with audited applicable finite/complete
outputs**, including coarse failed-accuracy cases, and 6 failed or rejected attempts.
One rejection is the expected VT guard; four are stale regression fixtures
(`reg_elastic` and `reg_poro_drained`, with both binaries), replaced by their
existing corrected `*_v2` fixtures. Undefined pressure in nonporous elastic
receivers is NaN by contract and is excluded from finite-pressure checks.
The four-node fine acoustic attempt in 93195 exhausted GPU memory during
setup, before timestepping, and was stopped before another expensive case.
The sealed/open cases are rerun separately on eight nodes (93212/93224). Completed-run
counts are not counts of physics passes.

- Baseline executable SHA-256: `3a54110ae89287fc648463ea1e812793b453fde1fe82ca2dd1224581f38b1770`.
- New v01 executable SHA-256: `9bbf94dffbfc77e5ec9021c88d863a1dcd815919062900fcdce8c2ed976a5e2f`.
- `psolve.cc` SHA-256: `fc99b14fbcc36ba2717e10af8f24c7bb891ec67d70e6f77d4a7cc23e6e0dcf80`.
- `poroelastic_biot_setup.h` SHA-256: `c67a217f6988ebac4a17e0cb5bda6b28c0bf77d27c88ba1278f955091cca02a9`.

Builds/runs use `ml reset`, the default Horizon modules, required GSL/PROJ paths,
`cuda-615-SC`, and the supported launcher. Libraries are checked before jobs.
Inputs, hashes, references, raw traces, job metadata and diagnostic/PDF helpers
remain ignored/local-only and are unavailable in a clean checkout. Existing
benchmark inputs and production binaries are preserved. The Horizon build,
selected GPU/MPI benchmarks, source checks, frozen-lockfile install, lint,
TypeScript, production documentation audits and GitHub Pages export audits pass.

## VT topography follow-up (2026-11-13)

This extends the full-Biot implementation from source baseline `77c85c6` to
`computation_method=vt`. It supersedes the earlier statement that full-Biot VT
is rejected. The earlier flat/FEM validation retains its historical scope.
The [VT PDF addendum](../../doc/Hercules_Full_Biot_VT_Topography_Report_2026-11-13.pdf)
contains implementation details, waveforms and explicit accuracy limitations.

### Discretization and supported scope

Native VT divides each cut cube into five P1 tetrahedra, using its existing
partition, retained volume fractions and minimum-volume cutoff. For each
retained tetrahedron with volume Vt and fraction ft, each vertex receives
`ft*Vt/4` for fluid inertia/drag and pressure capacity. The pressure coupling
uses `Q[a,b,d] = ft*Vt*grad(Na)[d]/4`; continuity is its negative transpose.
The tetrahedral pressure stabilization is symmetric, preserves constants and
has the same retained-volume geometry as the solid stiffness. Full-cell
operators remain Hex8. Physical saturation includes retained VT solid volumes
whose cube centers lie above the continuous terrain.

A drained surface has natural zero excess fluid pressure. An undrained surface
restricts relative-fluid motion to terrain tangent modes at supporting cut-cell
nodes. Normals use the bilinear height-grid derivatives. Zero-permeability axes
and sealed interface/domain locks are intersected with those modes. Setup
diagonalizes anisotropic drag in the allowed subspace; the compact CUDA update
includes the corresponding fluid-inertia reaction in solid momentum.
**This is a boundary approximation in an O(h)-thick support band, not exact
integration on reconstructed cut faces.** Shared geometric hanging constraints
also approximate varying normals. Near-surface hydraulic accuracy must be
measured separately from solid-motion accuracy.

VT retains shared geometric pressure/vector constraints. Full-Biot VT acoustic
coupling requires `acoustic_interface_projection=face_mass` for both hydraulic
options. The solid/fluid use reciprocal tetrahedral traces at the existing
four face quadrature points, including remote neighbors; acoustic fields retain
Hex8 traces. Partial-face sampling and discarded material need convergence.
Existing flat-water/complete-face requirements still apply. Field-specific
constraints and local matching remain unavailable with VT.

Pressure and relative-fluid receivers use tetrahedral barycentric weights in
cut cells. Full-Biot station setup preserves the owning cube's eight vertices
separately from the native four-vertex kinematic record; a receiver in a
discarded tetrahedron aborts with an actionable diagnostic. Binary station mode
continues to write the fluid companion. Sorted surface-cell property lookup
replaces the old per-cell scan without changing the selected geometry.

The setup enforces `dt <= 0.1*h/c_fast`, including when stabilization is zero,
and `dt <= h^4/(576*max(nu4))` for nonzero stabilization. The native normalized
gradient operator has row-sum bound `24*sqrt(max(nu4))/h^2`; matching positive
volume fractions cancel against nodal capacity. These initial-modulus bounds
do not certify all nonlinear tangent stiffnesses or resolve slow waves.

There is still no poroelastic PML. Porous CPML candidates, buried/spatial
water-table cuts, checkpoint/restart, geostatic ramps, fixed-base and DRM
histories remain unsupported. Mechanical dashpots remain approximate. Hydraulic
coefficients are global; this is not direct replication of all paper figures.

### Independent off-grid planar free surface

The surface is at solver Z=67.2 m, so the 8, 4 and 2 m meshes contain genuinely
partial tetrahedra. The domain is 2048 × 2048 × 1024 m; the solid force sheet
is at Z=192 m (surface-relative depth 124.8 m). A 12 Hz Ricker pulse is delayed
0.08 s; histories end at 0.32 s, with `dt=h/40000`. The force amplitude is
`2*rho*Vp_undrained*0.001` N/m². Large-domain pre-return windows exclude boundary
returns. The independent modal reference imposes zero total traction and
sealed normal flux or drained excess pressure. It uses no fitted amplitude or
time shift. Velocities are compared at `t-dt/2`, pressure at `t`.

Material: saturated density 1664 kg/m³, fluid density 1040 kg/m³, porosity 0.4,
tortuosity 2, drained lambda=5.28 GPa and mu=6.4 GPa, alpha=0.5226667,
M=6.0192616 GPa, permeability=6e-9 m², viscosity=0 (inviscid), pressure
stabilization=0.05. Fast/slow P speeds are 3590.329/1055.976 m/s.

The original bulk receivers are 32, 64 and 192 m below the surface. Maximum
raw relative L2 errors follow; targets are solid velocity ≤5% and pressure/flux
≤10%. Both finest bulk cases pass; coarse pressure/flux cases do not.

| Surface | h (m) | Velocity % | Pressure % | Darcy flux % | Job |
|---|---:|---:|---:|---:|---|
| undrained | 8 | 3.947 | 56.417 | 59.979 | 93290 |
| undrained | 4 | 1.643 | 30.267 | 32.535 | 93286 |
| undrained | 2 | 0.406 | 7.920 | 8.501 | 93280 |
| drained | 8 | 2.565 | 48.874 | 59.901 | 93290 |
| drained | 4 | 0.641 | 22.416 | 18.457 | 93286 |
| drained | 2 | 0.242 | 5.188 | 4.471 | 93281 |

An additional probe at depth 1 m tests the surface layer explicitly. Its sealed
normal flux remains under-resolved: **64.0% error at h=2 m fails the 10% target**.
Do not use bulk agreement to infer shallow hydraulic accuracy.

| Surface probe | h (m) | Velocity % | Pressure % | Darcy flux % |
|---|---:|---:|---:|---:|
| undrained | 8 | 2.527 | 85.387 | 100.000 |
| undrained | 4 | 0.528 | 18.369 | 71.842 |
| undrained | 2 | 0.155 | 5.375 | 63.996 |
| drained | 8 | 1.588 | 180.626 | 42.603 |
| drained | 4 | 0.371 | 38.438 | 13.494 |
| drained | 2 | 0.124 | 6.791 | 3.355 |

### Curved terrain, coupling and regression

The 512 m hill case uses elevation `64*exp(-r²/100²)` m, base Z=128 m and a
16 m height grid. A vector point force 96 m below the summit emits a 6 Hz pulse
delayed 0.2 s; five receivers at depths 16–64 m record through 0.4 s.
The table gives maximum self-convergence differences after physical-time
interpolation; there is no exact curved-terrain reference. Sealed pressure/flux
are not resolved to 10%, and sealed flux differences are not monotone.

| Surface | Meshes (m) | Velocity % | Pressure % | Darcy flux % |
|---|---|---:|---:|---:|
| undrained | 4 → 2 | 0.673 | 18.948 | 11.293 |
| undrained | 2 → 1 | 0.827 | 14.492 | 11.919 |
| drained | 4 → 2 | 0.528 | 5.187 | 1.493 |
| drained | 2 → 1 | 0.206 | 1.123 | 0.359 |

Coastal tests combine VT land, a flat water surface, sloping seabed and deeper
elastic material, with both acoustic scalar formulations and both hydraulic
conditions. These establish execution/transport consistency, not coastal
scattering accuracy. Coarse normal-incidence interface comparisons retain large
pressure/flux errors (up to about 71%/66%) and are not counted as accuracy passes.

- Host/CUDA-aware MPI and four/eight-rank comparisons, including adaptive and
  coastal cases: maximum nonzero-field relative L2 difference
  `1.4975e-12`. Roundoff-level zero fields are checked
  absolutely rather than divided by near-zero reference norms.
- Before/after elastic, acoustic, acoustic–elastic, reduced poroelastic,
  acoustic–poroelastic, Rayleigh and nonlinear elastic-regime regressions:
  maximum difference `4.4351e-12`. Selected native
  VT elastic and reduced cases are unchanged at written precision.
- Indexed-versus-original VT lookup runs: maximum difference
  `3.3128e-13`. Lookup affects setup cost only.
- Halving dt in the coarse vector-source hill changes velocity/pressure/flux by
  0.189%/0.159%/0.207%.
- Binary fluid output agrees with ASCII companions. Viscous/anisotropic drag,
  one impermeable coordinate axis, Rayleigh damping and nonlinear effective
  skeletons execute. The nonlinear test remains in the elastic regime;
  yielding, hysteresis and every BKT attenuation law are not validated.

### Build, tests and complete execution record

Horizon compilation used only `ml reset`, its default CUDA/MPI modules and the
required GSL/PROJ library paths. Production executables were untouched.
Jobs use account `ASC26075`, partition `debug`, reservation `cuda-615-SC`.
Investigation helpers and inputs are local-only; they are unavailable in a clean
checkout. No benchmark inputs from earlier work were changed.

Local non-CUDA checks passed:

- 20,000 affine-reproduction/partition-of-unity points, both partitions.
- 8,000 anisotropic constrained-mode cases.
- 1,000 cases extracted from the actual CUDA pressure/stabilization functions:
  adjoint work balance, symmetric positive-semidefinite stabilization, constant
  nullspace and the `24/h²` normalized row bound.
- 8,000 cases extracted from the boundary kernel: tangent/axis constraints and
  the reaction in coupled solid momentum.
- 1,501,500 indexed-versus-linear VT lookup comparisons.

The final documentation frozen install, lint, TypeScript, Sites production
build/rendered audits (73 pass, one export-only skip), and Pages export/audits
(74 pass) passed.

There are **105 completed runs**, plus 15 failed/rejected case
executions and six cancelled/replaced fine jobs. Completed runs are not all
accuracy passes. Earlier cut-cell receiver results from v01–v03 are retained
but superseded by the v04 output correction. v03 fine planar receivers are
outside cut tetrahedra, so those results are unaffected. The final v05 adds the
discarded-receiver guard and passes surface, adaptive, coastal, binary and valid
legacy regression checks. Final executable SHA-256:

`f1d6da37d10e92465474227fcd9d34352292beb73ecaf9612ba6e6bbed14ba22`

| Build | Changes / selected jobs |
|---|---|
| v01 | Initial VT operators; hill h4 93259, hill h2 93266; initial transport and regression checks |
| v02 | Remote solid face bounds corrected (off-face node IDs are dummies); coastal and backend jobs 93274–93277 |
| v03 | Indexed surface-cell lookup; planar h2 93280/93281, hill h1 93282/93283, identity checks 93284 |
| v04 | Correct cut-cell station cube metadata; surface/output checks 93286–93288 |
| v05 | Reject discarded-tetrahedron receivers; final checks 93290/93291/93293, baseline 93294 |

All failed and cancelled attempts remain visible:

| Job | Case | Reason |
|---|---|---|
| 93260 | `coast_open_pore_velocity_potential_shared` | Timestep exceeds conservative VT bound |
| 93269 | `coast_open_pore_pressure_dt` | Remote trace dummy-node geometry defect; fixed in v02 |
| 93269 | `coast_no_flow_pressure_dt` | Remote trace dummy-node geometry defect; fixed in v02 |
| 93269 | `coast_open_pore_velocity_potential_dt` | Remote trace dummy-node geometry defect; fixed in v02 |
| 93269 | `coast_no_flow_velocity_potential_dt` | Remote trace dummy-node geometry defect; fixed in v02 |
| 93272 | `vt_reduced` | Stale fixture rejected by parser; repaired before comparison |
| 93272 | `reg_elastic` | Stale fixture rejected by parser; repaired before comparison |
| 93272 | `reg_poro_drained` | Stale fixture rejected by parser; repaired before comparison |
| 93273 | `vt_reduced` | Stale fixture rejected by parser; repaired before comparison |
| 93273 | `reg_elastic` | Stale fixture rejected by parser; repaired before comparison |
| 93273 | `reg_poro_drained` | Stale fixture rejected by parser; repaired before comparison |
| 93284 | `reg_elastic_v2` | Stale fixture rejected by parser; repaired before comparison |
| 93284 | `reg_poro_drained_v2` | Stale fixture rejected by parser; repaired before comparison |
| 93285 | `reg_elastic_v2` | Stale fixture rejected by parser; repaired before comparison |
| 93285 | `reg_poro_drained_v2` | Stale fixture rejected by parser; repaired before comparison |
| 93262 | `offset_undrained_h4` | Cancelled to replace slow setup with equivalent indexed lookup |
| 93263 | `offset_drained_h4` | Cancelled to replace slow setup with equivalent indexed lookup |
| 93264 | `offset_undrained_h2` | Cancelled to replace slow setup with equivalent indexed lookup |
| 93265 | `offset_drained_h2` | Cancelled to replace slow setup with equivalent indexed lookup |
| 93267 | `hill6_undrained_h1` | Cancelled to replace slow setup with equivalent indexed lookup |
| 93268 | `hill6_drained_h1` | Cancelled to replace slow setup with equivalent indexed lookup |

Final solver source SHA-256 values:

| File | SHA-256 |
|---|---|
| `quake/forward/kernel.cu` | `b499f47a11b4d3ec6870bab9a17cfbd2e82dc71eeb9b459ce19eceec7f76c4fc` |
| `quake/forward/kernel.h` | `f03c90345f44b8d4b29eafae8124d861622fc747468666c90716d0b915db185d` |
| `quake/forward/poroelastic_biot.h` | `63f2f4185f90c697f56d6fa52268be8c9e8bbb1a2734a4c9e34c6a505e29d4f7` |
| `quake/forward/poroelastic_biot_runtime.h` | `2205c4f02c96e6b89367a7298e3c2a77f7e1a6023a99190e1cd1b7e00c749d26` |
| `quake/forward/poroelastic_biot_setup.h` | `3f11ca1de8eb10f91472301d51210acd7987bcee195b419817831c0edcbcc639` |
| `quake/forward/psolve.cc` | `33efe241f73a0cf77dc5ee3a0e205a8213d80892f1434d696b284fdb2dbbe8ec` |
| `quake/forward/psolve.h` | `6e354ae2d2a4b4c99fdfc27a8f1f0a62b826bd11b7dc7154dbd79b2031c513fb` |
| `quake/forward/psolve_gpu_launch.cu` | `985555ff76836815e79afaafbb6196a9029c8e6803e477ad55079ed8ad206ccc` |
| `quake/forward/psolve_gpu_launch.h` | `d1133ce43874fefc1dc50a1bc8635fb3ee727edf60113a1d7ee41264dd37b865` |
| `quake/forward/topography.cc` | `e150c2e2b73849f727e7f31a18acc1025647a76292b533cc253570e6f25f388e` |
| `quake/forward/topography.h` | `4fca3bf8179bffbb2ea3f608ef15a9fd87dd3448d36383c961202acba43b9668` |

## Full-Biot open-pore point-projection follow-up

Completed on Horizon, 2026-11-13, from baseline `7f90f69`. Full-Biot open pores
now accept `acoustic_interface_projection=point` (default) or `face_mass` on flat
and FEM-topography Hex8 meshes. The pre-existing CUDA kernels already exchange
water pressure into the solid/fluid momenta and return solid velocity plus Darcy
flux with matching spatial weights. The change removes the startup restriction
and prints the selected quadrature. No time-stepping kernels, field exchanges,
constraint maps or station interpolation changed. `acoustic_interface_coupling=yes`
remains mandatory. Full-Biot VT still requires `face_mass` for both hydraulic
conditions because only that path uses its tetrahedral trace weights.

[PDF report](../../doc/Hercules_Full_Biot_Open_Pore_Point_Report_2026-11-13.pdf).
This follow-up supersedes earlier statements requiring `face_mass` for all
full-Biot open-pore contacts. Earlier results retain their original scope.

### Analytical normal-incidence comparison

The earlier 12.5 Hz inviscid force-sheet setup and continuum modal reference are
unchanged: 2,048 m square sheet, source at solver Z=512 m, water to Z=256 m and
porous material to Z=1,024 m; delay 0.08 s, recording through 0.24 s. Receiver
coordinates are (1024,1024,Z) for Z=192,320,384,640 m. Water density is 1,000
kg/m3 and Vp is 1,500 m/s. The porous properties are the Table 6-derived
inviscid material above, with permeability 6e-9 m2 and viscosity zero; tortuosity
is 2 and stabilization 0.05. The force-sheet amplitude is
`2*rho*Vp_u*0.001 N/m2`; 0.001 m/s is the prescribed source-history scale.
Porous h=8/4/2 m uses dt=0.0002/0.0001/0.00005 s. Water cells are one level
finer. Raw comparisons honor the documented staggered velocity times, with no
fitted amplitude, fitted timing or output filtering. The source is adapted,
not the paper's combined explosive source.

Maximum receiver relative L2 errors for `point`:

| Porous h (m) | v (%) | p (%) | q (%) | Job |
| --- | ---: | ---: | ---: | --- |
| 8 | 6.9879 | 56.6241 | 65.6952 | 93336 |
| 4 | 1.1313 | 23.9797 | 21.6324 | 93337 |
| 2 | 0.2873 | 6.3502 | 5.3524 | 93353 |

`face_mass` gives the same errors to the reported precision; the largest
normal-incidence projection discrepancy is 6.747e-12.
The 2 m case meets the predefined 5% velocity and 10% pressure/flux targets;
8 m and 4 m fail at least one target. A nearly uniform face field makes midpoint
and four-point integration agree, so this comparison alone cannot distinguish
their behavior on spatially varying fields.

The reference script was derived and coded during this investigation from
continuum Biot modes, with boundary residuals and modal power balance checked.
It is independent of Hercules's discretization but not independently authored,
and has not been cross-validated with Gar6More3D. Author-provided Wolf reference
histories are separate and were not used as the reference for this table.

### Localized source, FEM terrain and longer records

A 512 m cube has water to Z=128 m, with the same porous material below.
A vector force (1e9,4e8,-2e8) N at (256,256,256) m has a 25 Hz Ricker history
delayed 0.04 s. Seven receivers are at (256,256,64), (304,288,64), (336,256,96),
(256,256,160), (304,288,160), (336,256,192), and (288,280,320) m. The 0.12 s
window precedes the shortest outer return (>0.12 s even using the porous fast
speed everywhere). These spatially varying fields exercise oblique arrivals;
there is no independent oblique-interface Green-function reference here.

Maximum receiver relative L2 differences between `point` and `face_mass`:

| Porous h (m) | v (%) | p (%) | q (%) |
| --- | ---: | ---: | ---: |
| 8 | 0.1282 | 0.1265 | 0.0319 |
| 4 | 0.0149 | 0.0216 | 0.0112 |
| 2 | 0.0064 | 0.0069 | 0.0031 |

The two projections converge toward each other, but the `point` 4 m versus 2 m
mesh differences remain 12.2844/71.7078/15.4127% in
v/p/q. The corresponding `face_mass` mesh differences are similar. In particular,
the localized pore-pressure field remains under-resolved. Small projection
differences do not demonstrate absolute scattering accuracy.

All three constraint methods (`field-specific-constraints`,
`shared-geometric-constraints`, `local-matching`) execute for the flat model.
The FEM coastal case combines dry topography, flat water, sloping seabed and a
deeper elastic region; shared and field-specific constraints are tested with
both acoustic formulations. Over 0.12 s the largest projection differences
are 0.295%/0.281%/0.420% for v/p/q. These are execution/transport and projection
comparisons, not validation against an exact coastal solution.

One-second records include boundary returns and repeated interface interaction:

| Model | v difference (%) | p difference (%) | q difference (%) |
| --- | ---: | ---: | ---: |
| Localized flat model | 2.4517 | 1.9540 | 0.8785 |
| FEM coastal model | 3.0681 | 1.9284 | 2.8791 |

Both projections remain finite. At sampled receivers, every late-window
(0.5–1 s) v/p/q peak is below 44% of its corresponding early (0–0.5 s) peak.
This establishes no observed runaway growth over these records, not a general
stability theorem or a validation of the approximate mechanical dashpots.
There is still no poroelastic PML.

### Algebra, transport, output and regression checks

- 4,800 tests execute the extracted production interface kernels on the host,
  checking paired spatial work for both projections/formulations, six faces,
  matching/nonmatching fragments and local/remote solid vertices. Maximum
  normalized residual is 2.541e-16. CPML is excluded. This is spatial adjoint
  consistency, not a fully discrete energy-conservation proof.
- Twenty-one one/four/eight-GPU and host/CUDA-aware MPI comparisons agree within
  6.664e-13 relative L2 for nonzero fields. Near-zero quantities use
  absolute checks; no signal filtering is applied.
- Eleven before/after regressions agree within 5.967e-12. They cover
  elastic and acoustic standalone paths, acoustic–elastic, reduced/full
  poroelastic paths, acoustic–poroelastic contacts, attenuation and an eligible
  nonlinear skeleton in its elastic regime. The coastal setup also includes
  poroelastic–elastic contact. Yielding, hysteresis, every attenuation family,
  CPML absorption and general heterogeneous hydraulic interfaces remain gaps.
  All 406 comparable station columns, including strain/stress/nonlinear state,
  pass relative 1e-9 or absolute 1e-12 tolerances, with identical NaN masks.
- Acoustic pressure and velocity-potential runs agree within
  3.431e-12. Halving the localized 4 m timestep changes v/p/q by
  0.499%/0.054%/0.025%. A viscous open-pore integration case runs; this is not an
  analytical viscous-interface accuracy test.
- Binary fluid output agrees within 4.94e-13; binary solid displacement agrees
  with ASCII within its printed precision (1.36e-7 relative L2).
- Coupling disabled and VT point cases abort with their intended diagnostics.

### Run and source provenance

86 runs completed with finite applicable station quantities, plus
two intentional guard rejections. Submission 93352 failed before solver startup
because inputs were missing; its two planned cases were resubmitted as 93353.
No failed attempt is counted as a completed validation.

Jobs: 93336–93346 (excluding 93347), 93348–93350, 93353 and 93362.
All use the `cuda-615-SC` reservation and account `ASC26075`. Default modules
were selected with `ml reset`; required GSL/PROJ paths were retained and binary
library dependencies checked. The isolated executable is `psolve-v01`:

```text
eb1726b49d3efe4fa39dc45198398accb09f0d280d4316688d69d6316014484b
```

All 94 tracked solver-source hashes match the local files. The
previous VT executable is retained separately as `psolve-base`. Production
executables were untouched.
Local-only inputs, scripts, source/input hashes, job metadata and raw histories
are in the private validation archive and are unavailable in a
clean checkout. The local README describes preparation, submission, analysis,
outputs and overwrite behavior; no investigation helper is published in `script/`.

### Exact Wolf-paper reproduction: remaining work at the point-projection stage

Full fluid inertia, zero-viscosity flow and natural drained conditions are now
available. Reproduction is feasible in stages; the current tests are adapted
verification problems and do not constitute reproduction of Figures 5–8.

| Target | Existing reference | Required next work |
| --- | --- | --- |
| Figure 5, homogeneous full space | Author-supplied viscous/inviscid Green-function histories | Add the simultaneous isotropic normal-stress and pore-pressure source with exact normalization, sign and temporal convention; match all eight receivers, material values and record length. Start here. |
| Figure 7, flat free surface | Author-supplied Gar6More3D histories at 0.5 and 500 m depth | Use the matching combined source and zero total traction/pore pressure; refine shallow response independently of bulk motion. |
| Figure 6, two porous half-spaces | Author-supplied inviscid Gar6More3D solid velocities | Add material-dependent hydraulic properties and verify porous–porous mechanical/hydraulic continuity. Current hydraulic controls are global. |
| Figure 8 geometry, Figures 9–11 layered results | Author-supplied finite-difference histories, including 0.625 m reference | Add the hydraulic material map and combined 401-point line source with 50 m spacing, prescribed Gabor pulse and 5 s duration; verify near-surface slow-wave resolution. |

The paper's source file has diagonal stress entries -1 and pore pressure +1,
scaled by 1e10. Converting it into u–w–p forcing must preserve units and the
source time integration convention; a mechanical point force is not equivalent.
The paper defines relative fluid velocity as porosity times the intrinsic
fluid-minus-solid velocity, matching Hercules Darcy flux `q`, not `q/porosity`.
The archive conversion scripts also apply reference amplitude scales and, in
the free-surface case, a prescribed 0.001 s time-origin adjustment. Those
normalization conventions must be traced explicitly before comparison; they
must not be replaced by fitted waveform amplitudes or fitted time shifts.
For the viscous Table 6 test the permeability is 6e-13 m2, not the adapted
6e-9 m2 used above. Inviscid propagation with positive permeability and zero
viscosity is independent of that drag coefficient, but the viscous test is not.
A poroelastic PML is not necessary if the domain excludes boundary returns
throughout the comparison window. Lower-order Hercules meshes may need much
finer resolution than SeisSol's high-order discretization; equal element size
does not imply equal accuracy.

The [author archive](https://zenodo.org/records/5236133) is already downloaded and
checksum-verified. No additional download is needed to begin this work. The
source and geometry contracts are described in
[Wolf et al.](https://doi.org/10.1016/j.jcp.2021.110886), sections 6.2–6.5.

## Poroelastic CPML follow-up (2026-11-13–14)

Implemented a mixed Cartesian complex-frequency-shifted CPML for both `brick_up`
and `full_biot`. Saturated cells now retain their porous degrees of freedom in
the absorbing layer. The implementation and Horizon validation are described in
[the PDF report](../../doc/Hercules_Poroelastic_CPML_Report_2026-11-14.pdf) and
[the operational reference](../../doc/Extended%20README.md#poroelastic-cpml).
This is an adapted absorption validation, not reproduction of Wolf Figures 5–8.

### Operator and supported scope

For stretches `S_j=kappa_j+d_j/(s+alpha_j)` and `J=S_x*S_y*S_z`, the weak
operators use `J` for coupled inertia and pressure storage, `J/(S_j*S_l)` for
stiffness and reduced Darcy diffusion, and `J/S_j` for pressure/volume coupling.
Full Biot retains solid/fluid cross inertia and `J*s` Darcy drag. Exact local
leading drag relaxation is retained. Pressure storage advances its convolution
history while sharing one nodal pressure across element profiles. These are
mixed weak-form adaptations, not an elastic absorber applied only to the solid.

The existing elastic CPML supplies drained stiffness and solid diagonal inertia.
Additional GPU histories number 372 scalars per full-Biot layer cell and 156 per
reduced cell. Existing scalar/vector MPI exchanges assemble contributions;
there is no extra per-step communication round or host staging in CUDA-aware mode.
Automatic reference speed includes porous fast/undrained speed. Reduced sealed
fast-path folding is disabled with CPML so pressure advances explicitly.
Both solid and relative-fluid displacements respect the outer clamp.

Open-pore water pressure now applies the same stretched traction to solid and
fluid, and acoustic return flux includes both motions. Both acoustic scalar
formulations and flat/FEM `point` and `face_mass` projections are exercised.
VT full-Biot acoustic coupling retains its `face_mass` requirement. Poroelastic–
elastic no-flow contacts and the existing effective-skeleton forces remain.
Nonlinear constitutive and BKT histories remain excluded from the CPML itself;
eligible nonlinear elastoplastic/plastic and hysteretic interiors retain their
existing behavior. The regression fixtures do not certify yielding/hysteresis.

Native field-specific VT maps use compatible P1 solid/full-Biot pressure traces,
Q1 acoustic traces and Q1 regular-buffer traces. Force distribution and field
assignment use identical positive weights. A face midpoint can select two
coarse diagonal vertices from the existing geometric halo; incompatible incident
traces or missing parents abort. Reduced BrickUP retains its existing geometric
VT pressure constraints and pressure interpolation. This is not a migration of
its cut-cell pressure discretization.

Required limits:

- `cpml_stiffness_mode=native`; reduced stiffness aborts with active porous CPML.
- Full-Biot `open_pore` with any porous CPML requires
  `acoustic_interface_method=field-specific-constraints`. Shared geometric maps
  showed late growth in FEM and VT; `local-matching` is not qualified for this
  combination. The guard also applies when the water interface is outside CPML.
- Porous VT cut cells cannot intersect the layer. Use a regular Cartesian buffer
  with flat mesh-aligned terrain through CPML. This is checked at setup.
- Buried/spatial groundwater, geostatic ramps and checkpoint/restart are
  unsupported for porous CPML. Existing full-Biot fixed-base/DRM exclusions remain.
- The full-Biot pressure-stabilization spatial operator is retained, with
  stretched storage normalization. This finite-mesh dissipation is not an exactly
  matched continuum term. Refine and check timestep/stabilization sensitivity.
- Layer pressure and fluid outputs are diagnostics, not physical receiver
  predictions; `Ru` is undefined there. The partition weight is not calibrated
  for the added porous-layer memory/work.

### Same-grid boundary comparison

Use a 512 m cube and an enlarged reference, with identical near-source mesh,
source, receivers, material and top surface. Core references are 2,048 m cubes;
source and receivers are translated by 768 m in X/Y, with unchanged Z. The
25 Hz Ricker vector force is `(1e9,4e8,-2e8)` N at `(256,256,256)` m, delayed
0.04 s and recorded for 0.5 s. Velocity samples use their prescribed `t-dt/2`
times. No fitted shifts, amplitude correction or output filtering is applied.
Errors below are worst receiver relative L2 differences against the enlarged
numerical reference, not absolute continuum accuracy or a reflection coefficient.

Synthetic material: mixture density 1,664 kg/m3; fluid density 1,040 kg/m3;
porosity 0.4; tortuosity 2; drained lambda 5.28 GPa, mu 6.4 GPa; alpha
0.5226666667; M 6.019261637 GPa. Fast/slow speeds are 3,590.329/1,055.976 m/s.
The inviscid test has positive permeability and zero viscosity; the viscous test
uses viscosity 0.001 Pa s and permeability 6e-13 m2. Sealed permeability is zero.
The reduced comparison uses viscosity 0.001 Pa s and permeability 6e-9 m2.

CPML profile: 64 m thickness, reflection target 1e-4, power 3, kappa 1, alpha
6.25 Hz, no alpha floor, multiaxial cross damping or acoustic sponge.
These are case-specific settings, not new defaults. Timesteps are `h/80000` s.

| Case | h (m) | Solid velocity % | Pressure % | Darcy flux % |
| --- | ---: | ---: | ---: | ---: |
| Inviscid CPML | 16 | 9.929 | 7.281 | 10.061 |
| Inviscid CPML | 8 | 0.835 | 0.581 | 0.844 |
| Inviscid CPML | 4 | 0.191 | 0.315 | 0.252 |
| Inviscid dashpot | 4 | 13.763 | 37.247 | 24.988 |
| Viscous CPML | 8 | 0.877 | 0.604 | 1.159 |
| Sealed CPML | 8 | 0.842 | 0.576 | identically zero |
| Reduced CPML | 8 | 0.849 | 0.582 | no independent fluid state |

Normal incidence is tested separately with a horizontal force sheet in a
2,048 m wide domain, small/reference depths 512/2,048 m, source Z=384 m,
receivers Z=352,368,400,416 m and 0.25 s records. Fast/slow modal projections
use the independent continuum generalized eigenvectors. Maximum modal boundary
differences are 0.681% fast P and 1.638% slow P; S-wave velocity/flux differences
are 0.565%. S-wave pressure is analytically zero; its observed maximum below
4.8e-11 Pa is roundoff and is not assigned a relative error. A near-grazing
point-force geometry gives 1.920% velocity, 2.483% pressure and 1.916% flux.
This is a finite set of incidence tests, not an all-angle guarantee.

At 8 m, increasing thickness to 128 m gives 0.109%/0.569%/0.224% differences.
Kappa 2 gives 1.001%/0.565%/1.008%; multiaxial cross fraction 0.1 gives
7.585%/5.894%/7.275%. The default profile gives 14.041%/1.396%/14.020% on this
particular pulse. A requested reflection target is not a measured error guarantee.
Halving timestep changes the 8 m histories by 0.122%/0.162%/0.122%.
Removing stabilization changes pressure/flux substantially; that comparison
changes the volume discretization too and is not a pure boundary-error measure.
Strong drag, anisotropic permeability and drained-surface controls are also run.

### Water, coast and the stability correction

At 8 m porous spacing, a flat 128 m water layer is compared against the
corresponding enlarged reference. Open-pore `face_mass` maximum boundary
differences are 1.082% velocity, 1.210% pressure and 0.267% Darcy flux;
`point` gives 1.077%/1.210%/0.267%. No-flow `face_mass` gives
1.060%/1.142%/0.258%. Both acoustic formulations agree within 3e-9 relative L2
in the recorded controls. This does not establish exact oblique interface accuracy.

The coastal case combines a dry curved hill, a flat water surface, sloping
seabed, porous layer and deeper elastic region. Terrain cuts end inside the
physical domain, leaving the required flat buffer. A 1,024 m enlarged reference
preserves the local mesh and geometry by translating source/receivers 256 m in
X/Y. Through 0.25 s, VT differences are 0.179%/0.310%/0.150%; FEM differences
are 0.186%/0.499%/0.155%. These isolate boundary sensitivity, not terrain error.

**The shared-map open-pore trials fail long-time stability.** Initial 0.5 s
records were finite, but the VT case grew by roughly 1e7 relative to initial
receiver peaks over 3 s; a FEM shared-map control also grew. Half/quarter timesteps,
removing stabilization and alternate profiles did not resolve it. M-PML 0.1
suppressed short-record growth but regrew by 10 s. Stronger cross damping could
be worse. Those histories remain archived and are not successful validations.

With native field-specific VT constraints and ordinary CPML (cross fraction
zero), the same coastal model remains bounded through 10 s. Across all receivers,
last-0.5-s / first-0.5-s peak ratios are at most 0.1734% solid velocity,
0.01182% pressure and 0.09042% flux. At receiver 0 they are
0.03220%/0.004921%/0.03105%. A drained VT field-specific case remains bounded
through 3 s; corresponding all-receiver ratios are 1.509%/0.1365%/1.617%.
Flat full-Biot inviscid/viscous, flat water and field-specific FEM cases also
have decaying late peaks over their 3 s records. These are finite-duration
observations, not a general eigenvalue or discrete energy stability proof.

### Algebra, transport, regression and provenance

Host algebra checks cover 300 random profiles with 200 steps each: exact cube
pressure-gradient integrals agree with independent Gauss quadrature within
8.33e-17; rational inertia/drag identities within 1.58e-14; discrete stretched
storage/fluid-content identities within 9.79e-18. Extracted production coupling
kernels pass 4,800 spatial adjoint-work checks with nonunit cofactors, both
projections/formulations, six faces and local/remote matching/nonmatching traces;
maximum residual is 2.33e-16. These are spatial/filter checks, not a proof of the
complete time integrator. Native VT traces pass 108 support/affine probes.

The v04 final matrix compares native full-Biot and elastic VT on one/four/eight
GPUs and host/CUDA-aware MPI: maximum aggregate relative L2 1.19e-11.
Fifteen baseline regressions differ by at most 4.29e-13, covering standalone
elastic, acoustic, reduced/full porous, acoustic–elastic, porous–water,
viscoelastic and nonlinear elastic-regime paths; four include nonporous CPML.
All 782 comparable station-column checks pass relative 1e-9 or absolute 1e-12
criteria, including strain/stress/state columns and applicable pressure/flux.
Near-zero reference columns are judged by the absolute criterion. Selected
final v02/v03 repeat checks preserve the earlier CPML results. The v05 change
preserves reduced VT pressure interpolation while retaining native full-Biot
maps. Unchanged-full-Biot/elastic/reduced-flat repeats and reduced-VT
eight-rank/host controls agree within 2.43e-10 aggregate relative L2. The reduced
VT drained 3 s case has last/first 0.5 s peak ratios no greater than 3.34%
velocity and 0.188% pressure; these are bounded-execution checks, not independent
reduced coastal absorption accuracy.
The final station diagnostic checks receivers at X=32,60,68 m around the
64 m layer boundary: pressure and fluid histories are finite and all three
are correctly labeled porous; the first two are absorbing-layer diagnostics.
No fixture calibrates yielded, liquefying or arbitrary hysteretic response.

Builds use only `ml reset`, the existing GSL/PROJ paths and library-resolution
checks. Final isolated executable SHA-256 (`psolve-v05`):

```text
14acee59b5b7100ec6138d142fbf7aa34674fa96d41e478f96258f4ff567f422
```

Earlier core/profile measurements use v01/v02; native VT long-time/reference
measurements use v03; guarded regression/MPI measurements use v04. Their
provenance, input snapshots and job lists are retained with each run. v04 SHA:
`8a9b4cf8f8f51afe55627027d47a6c7a5491d7c1150529b6dcc119f435fd08bd`;
v03 SHA: `98525e047a884ad33cd03e310d572af32a173c3641d7d39335cf7ca5b2cf4710`.
All 95 final solver-source/build-file hashes match the local sources.
Production executables were untouched.

Jobs and exact case/binary/rank/backend associations are in the private `jobs.json`
manifest, together with input
snapshots, raw station histories, source hashes, metrics, failed attempts and
helper usage. They are unavailable in a clean checkout. The reservation is
`cuda-615-SC`, account `ASC26075`. Failed fixture preparation is distinguished
from physics: eight original coastal inputs were rejected before propagation
because terrain crossed water free-surface cells; four elastic-only inputs
retained an inapplicable `open_pore` setting and were corrected; job 93471 lacked
prepared reference inputs. One linear output fixture requested nonlinear-only
`Ru` and was rejected before propagation; the corrected pressure-only fixture
is recorded separately. Six recorded guard rejections include two repeated
checks; the four final checks reject reduced CPML stiffness, porous VT cuts,
shared open-pore constraints and unqualified local matching. Finite-output audits
include retained unstable profile trials and must not be reported as universal
stability passes.

The formulation follows coordinate-stretch principles from
[Martin, Komatitsch and Ezziani (2008)](https://doi.org/10.1190/1.2939484), with
second-order finite-element CPML guidance from
[Xie et al. (2014)](https://doi.org/10.1093/gji/ggu219) and
[Matzen (2011)](https://doi.org/10.1002/nme.3205). The mixed pressure-storage and
coupling adaptation above is derived for Hercules and tested separately; it is
not a copied implementation from those papers. CPML now permits smaller
validation domains. At this stage exact Wolf reproduction still needed the
combined stress/pore-pressure source; the author-reference follow-up below
adds that source. Heterogeneous cases still need hydraulic material fields.
The author archive is already available; no additional download is needed now.

## Controlled 128 m CPML / 10 s follow-up (2026-11-14)

The requested homogeneous 128 m case remains bounded for 10 s. However, the
shared-constraint coastal open-pore case still grows with 128 m CPML, so the
normal solver's field-specific-constraint requirement is retained. The
[follow-up PDF](../../doc/Hercules_CPML_Thickness_Followup_2026-11-14.pdf)
contains the envelopes and comparison protocol. No production solver code or
input interface changed in this follow-up.

### Results and decision

Each entry is the largest, over applicable receivers, of the maximum absolute
component during 9.5–10 s divided by that receiver's peak during 0–0.5 s.
Pressure is scalar; velocity and Darcy flux use their three components.
Different receivers can control different columns. These ratios describe
late growth/decay, not waveform accuracy or reflection coefficients.

| Case | Velocity late/early | Pressure late/early | Darcy flux late/early |
| --- | ---: | ---: | ---: |
| Homogeneous, 128 m | 0.001873 (0.1873%) | 0.00008363 (0.008363%) | 0.001873 (0.1873%) |
| Coastal shared constraints, 64 m | 2.109e38 | 2.476e38 | 2.061e38 |
| Coastal shared constraints, 128 m | 2.777e9 | 1.976e9 | 1.652e9 |

All three completed their 100,000 time steps with finite applicable station
quantities. The two enormous coastal responses are unstable failures despite
being finite; they are not physically meaningful predictions. The homogeneous
case has no acoustic interface and cannot by itself qualify open-pore coupling.

The 128 m coastal layer delays and slows the growth. For the receivers
controlling the table, 0.1 s window peaks first exceed their initial peaks
around 4.15–4.25 s, versus about 1.45 s for 64 m. Descriptive logarithmic fits
to 7–10 s envelopes give growth rates about 3.71/s versus 10.31/s; these are
measured signal slopes, not independently computed operator eigenvalues.
A three-second run of the thicker case would have missed the later failure.
The result does not establish that every still-thicker layer is unstable, but
it rules out using this 128 m / 10 s experiment to remove the current guard.

### Controlled inputs and mesh checks

The homogeneous case extends `control_thick` from 0.5 to 10 s, with identical
512 m domain, 8 m cells, material, source and 0.0001 s timestep. Its first
0.5 s agrees with the previous run within 3.93e-13 relative L2.

The coastal comparison retains the original dry hill, flat water surface,
sloping seabed, porous–elastic depth interface, full Biot, inviscid fluid,
`shared-geometric-constraints`, `face_mass`, velocity-potential acoustics and
0.0001 s timestep. CPML profile controls are identical: reflection target 1e-4,
power 3, kappa 1, alpha 6.25 Hz, zero alpha floor, zero multiaxial cross fraction
and zero acoustic sponge. The profile's damping scale is recomputed for each
thickness by the existing algorithm. Both layers use 8 m porous/solid cells
and 4 m water cells: 8 versus 16 porous-cell widths and 16 versus 32 water-cell
widths. No smoothing, fitted amplitudes or time shifts are applied.

The original 64 m case occupies 512 x 512 x 512 m. For 128 m, extend each lateral
exterior and the bottom by 64 m, giving a physical 640 x 640 x 576 m model;
translate terrain/source/receivers by +64 m X/Y. Thus the inner CPML faces stay
at the same physical locations relative to the source and receivers, including
original bottom-interface Z=448 m. Terrain cuts remain outside CPML.

The initial 640 x 640 x 576 m fixture failed octree extraction before propagation
(job 93562). The corrected fixture translates the entire physical model +448 m
Z inside a 640 x 640 x 1024 m octree box. `topographybase_zcoord` changes from
128 to 576 m; the added upper space is air and surface-relative source/receiver
depths remain unchanged. This supplies an enclosing 1024 m root with
`maximum_octant_level=7`, preserving 8 m cells. No octree code was modified.

After undoing translation, receiver coordinates agree within 1e-6 m and all
receiver medium, cell-size and material metadata agree. The two coastal early
records through 0.08 s agree within 6.41e-8 relative L2. The 64 m rerun agrees
with the previous three-second run within 1.99e-13. These checks support a
matched physical comparison rather than moving terrain into the thicker layer.

### Executables and archive

Horizon jobs 93560 (homogeneous), 93561 (64 m coast) and 93564 (128 m corrected
coast) each used four GPUs, CUDA-aware MPI and the `cuda-615-SC` reservation,
account `ASC26075`. Only `ml reset` was used; library resolution and required
GSL/PROJ paths were retained. Job 93562 is a failed mesh preparation attempt,
not a propagation result.

The homogeneous run uses the normal v05 executable, SHA-256
`14acee59b5b7100ec6138d142fbf7aa34674fa96d41e478f96258f4ff567f422`.
The coastal runs use an isolated copy of the current source whose sole code
difference disables the shared-constraint startup rejection. Every physical
operator, including CPML and VT coupling, is unchanged. All 95 source/build-file
hashes match the current source plus that recorded one-guard diagnostic diff.
Diagnostic executable SHA-256:

```text
3d120a83acee8acbcb0b0e6d1e5fca29d121a115359b568f4a4f56e7a75a0edd
```

The bypass is local-only and is not committed or installed as the normal solver.
Local-only fixtures, hashes, jobs, raw output, comparison metrics and scripts are
in the private validation archive. These helpers are
unavailable in a clean checkout; their README records dependencies, inputs,
commands, outputs, overwrite behavior and limitations.

## Wolf author-reference comparisons (2026-11-14)

The matching source is implemented and the original author histories are now
used directly for Figure 5 (inviscid and viscous full space) and Figure 7
(drained flat free surface). **These tests do not yet reproduce the paper's
accuracy.** They are distinct from the earlier adapted point-force tests.
The [nine-page report](../../doc/Hercules_Wolf_Reference_Validation_2026-11-14.pdf)
contains waveform overlays, all controls, source derivation and limitations.
No further reference download is needed for these cases.

### Source contract and reference provenance

The archive README attributes full-space Green functions to Martin Galis using
Karpfinger's code, free-surface/interface solutions to Gar6More3D, and LOHp
histories to the Galis/Gregor/Moczo finite-difference code. The source's diagonal
stress entries −1 and pressure +1, with factor 1e10, prescribe **first-order
rates**, not instantaneous displacement-form forces. The new source selects
`type_of_source=biot_explosion`, with required finite signed coefficients
`biot_stress_moment_nm` and `biot_pressure_moment_nm`, both 1e10 N m here.
For tension-positive total stress and compression-positive pore pressure:

```math
\dot{\boldsymbol\sigma}^{src}=-M_s S(t)\mathbf I\delta(\mathbf x-\mathbf x_s),
\qquad \dot p^{src}=M_p S(t)\delta(\mathbf x-\mathbf x_s),
\qquad H(t)=\int_0^t S(\tau)\,d\tau.
```

The existing discrete `forcefunction.in` layout supplies rate samples S in
1/s; the source integrates their piecewise-linear interpolant exactly once,
starting at zero and holding the final integral. Both moments use the same
containing regular Hex8 element:

```math
\mathbf f_a^{src}=(M_s-\alpha M_p)H(t)\nabla N_a(\mathbf x_s),\qquad
\Delta r_a^p={M_p\over M}N_a(\mathbf x_s)[H(t_{n+1})-H(t_n)].
```

Subtracting alpha times the injected pressure moment prevents counting it
twice in total stress. The pressure load enters the actual continuity RHS
before existing constraint/MPI assembly; it is not a receiver correction.
Solid loading uses the existing CPU/GPU source path and pressure adds a small
GPU kernel on its owning rank, without another communication round.

The Figure 5 source is the original 16 Hz Ricker rate sampled every 0.001 s,
with peak at 0.07 s. All eight receivers at ±600/±1000 m on X and ±345/±575 m
on each diagonal coordinate are compared through 1.5 s; the singular source
receiver is excluded. The author converter's scale 7.025e12 is fixed. Table 6
uses rho_s=2080, rho_f=1040 kg/m³, porosity 0.4, tortuosity 2,
lambda_d=5.28e9, mu=6.4e9, Ks=2e10, Kf=2.5e9 Pa, permeability 6e-13 m²,
and viscosity 0 or 0.001 Pa s. Input saturated density is 1664 kg/m³,
Vp_d=3296.2682863 and Vs=1961.1613514 m/s; alpha=0.5226666667.
Relative fluid velocity is the volume-weighted Darcy flux q, not q/porosity.

Figure 7 uses the original 5 Hz rate, peak 0.25 s, source depth 500 m,
zero total traction and zero pore pressure at the top, and a 2 s record.
Receivers (500,400,0.5) and (500,400,500) m are compared. Table 7 uses
rho_s=2500, rho_f=1040 kg/m³, porosity 0.2, tortuosity 2,
lambda_d=mu=1.2e10, Ks=4e10, Kf=2.5e9 Pa, permeability 6e-13 m²,
zero viscosity, saturated rho=2208 kg/m³ and alpha=0.5.

Gar6More3D Ux/Uy/Uz files contain **displacements**. Their fixed author scale is
1e10/(2 pi² 5²), and the prescribed time adjustment is +0.001 s.
The author converter imports a helper absent from the archive. Here velocity
is a centered NumPy gradient of displacement; an identical author derivative
convention cannot be claimed. Direct displacement and archived SeisSol
velocity comparisons provide additional checks. The small reference artifact
around 1.5 s noted in the paper is retained in full-window errors.

### Error measures and controlled runs

Vector relative L2 is the time-history error norm divided by the reference
norm, without fitted amplitude, time shift or filtering. Tables report the
largest receiver value. Displacements use output t; velocities and q use
t−dt/2 and are interpolated to author sample times. Diagnostic EM/PM in the
PDF use a Morlet wavelet, w0=6, 80 logarithmic frequencies (1–40 Hz for
Figure 5; 1–12.5 Hz for Figure 7); their processing is not claimed identical
to the paper. Undefined relative errors for zero components are not used.

| Figure 5 case | h, m | Stabilization | Solid v L2, % | Darcy q L2, % |
| --- | ---: | ---: | ---: | ---: |
| Inviscid, graded | 32 | 0.05 | 304.84 | 100.86 |
| Inviscid, graded | 16 | 0.05 | 113.60 | 100.52 |
| Inviscid, graded | 8 | 0.05 | 111.87 | 109.38 |
| Inviscid, graded | 8 | 0 | 130.83 | 361.14 |
| Inviscid, graded | 8 | 0.005 | 116.39 | 128.40 |
| Inviscid, CPML | 8 | 0.005 | 97.03 | 128.39 |
| Inviscid, CPML | 4 | 0.005 | 116.41 | 149.70 |
| Viscous, graded | 32 | 0.05 | 551.92 | 857.44 |
| Viscous, graded | 16 | 0.05 | 161.40 | 199.68 |
| Viscous, graded | 8 | 0.05 | 111.94 | 131.91 |
| Viscous, CPML | 8 | 0.005 | 38.91 | 52.91 |

The graded full-space cube is 8192 m with a central ±1280 m fine region and
coarsening to 256 m. The external boundary is remote enough, but mesh
transitions return reflected energy during the comparison. The CPML control
uses a 4096 m cube with uniform interior/CPML spacing except for the source
patch; layer thickness is 256 m, target 1e-6, kappa 1, zero cross damping,
alpha 16 Hz (5 Hz for the free surface), and stabilization 0.005. Its full-space
source near Z=3073 m puts the return from the unabsorbed top after 1.5 s.
The source neighborhood is refined to h/4. Figure 5 dt is 0.0002/0.0001/
0.00005/0.000025 s for h=32/16/8/4 m; Figure 7 dt is 0.00016/0.00008/
0.00004 s. The 4 m inviscid control has 1,076,947,390 elements on 128 GPUs;
the prepared 4 m viscous case was not run.

CPML reduces late reflected fast-wave contamination, especially in the viscous
case. Reducing stabilization alone does not recover the inviscid slow wave;
viscous results are almost unchanged because stabilization decreases with
drag. These are still failures to match the reference closely.

At 4 m the inviscid full-window maxima are **116.41% solid velocity and
149.70% Darcy flux**, larger than at 8 m. The slow wave gains amplitude but
remains delayed and rings. In a separate early-fast-arrival comparison against
archived SeisSol output, worst velocity error decreases from 32.61% at 8 m to
8.18% at 4 m. That window ends at receiver distance/c_fast + 0.16 s. This
arrival improves, but neither complete waveform is resolved. The 4 m pressure
difference from archived SeisSol is 149.57% relative L2. All 18 author-case
coordinate/material/source audits and configured output-count checks pass.

| Figure 7 CPML h, m | Receiver depth, m | Displacement vector L2, % | Velocity vector L2, % |
| --- | ---: | ---: | ---: |
| 16 | 0.5 | 94.11 | 109.52 |
| 16 | 500 | 81.64 | 100.42 |
| 8 | 0.5 | 33.74 | 43.57 |
| 8 | 500 | 30.93 | 46.47 |

Secondary comparisons against archived SeisSol pressure give maximum relative
L2 differences of 128.30% (inviscid full space) and 38.27% (viscous full space)
for the 8 m CPML runs. The 8 m free-surface pressure differences are 89.92% at
0.5 m and 64.56% at 500 m depth. These numerical reference histories are
additional evidence, not independent pressure Green functions.

### Resolution diagnosis and implementation checks

An independent continuous radial modal Green function checks absolute inviscid
source normalization against the author files within 0.593% solid and 0.613%
fluid relative L2, using a smooth Ricker rather than the sampled interpolant.
It uses no FE operators or fitted gain. This supports the source convention,
not the accuracy of the simulated propagation.

The continuum fast/slow speeds are 3590.329/1055.976 m/s, giving a 66 m slow
wavelength at 16 Hz. On a uniform axial inviscid unstabilized grid, the lumped
Hex8 elastic symbol is 4 sin²(kh/2)/h² and the collocated pressure-gradient
symbol is sin(kh)/h. The resulting slow branch has no 16 Hz propagating mode
at 32 or 16 m. At 8 m its phase/group speeds are 933.27/695.92 m/s; at 4 m
they are 1030.14/978.80 m/s, giving about 2.39 radians of phase error over
1 km even at 4 m. This semidiscrete axial calculation excludes stabilization,
time stepping, source patches, CPML and obliquity. It explains a resolution
mechanism without attributing every measured error to it.

The 23 original-source regressions agree with preceding binaries within
2.86e-10 relative field L2. The production cumulative-history function passes
2460 evaluations against independent adaptive quadrature, maximum absolute
error 1.20e-15. Separate stress/pressure source runs superpose to the combined
open-pore result within 3.57e-7 (ASCII output precision). New-source tests cover
host/CUDA-aware MPI, 1/4/8 ranks, pressure/velocity-potential acoustics, FEM/VT,
reduced/full Biot, attenuation and a below-yield nonlinear case. This is
integration evidence, not independent validation of nonlinear yielding or
all acoustic/porous scattering combinations.

CPU VT testing exposed stale interior displacement data: boundary-only host
visibility gave about 48% field error with the new source and 90–105% with
the old force source. The correction selects full host visibility whenever
CPU VT is active, even if boundary mode was requested; default GPU execution
remains device-resident. Sixteen final CPU/GPU coupled comparisons agree within
4.21e-10. A separate exclusion makes CPU VT skip cells already assigned to
Cartesian CPML, matching the GPU mapper and preventing duplicate stiffness.

An isolated elastic VT fixture (jobs 93781/93782/93783) contains 1155 cells
reassigned to Cartesian CPML. Keeping full host copies in both CPU controls,
removing duplicate stiffness reduces CPU/GPU field differences from 4.60% to
2.65e-11 without damping and from 2.62% to 1.55e-11 with Rayleigh damping.
This isolates the stiffness correction from the host-copy correction.
Porous VT cuts within CPML remain rejected; no geometry guard is relaxed.

The source requires explicit pore pressure (`brick_up` or `full_biot`), exactly
one owning regular porous Hex8 cell outside CPML/VT topographic tetrahedra,
discrete unfiltered finite input, and no checkpoint/restart/geostatic delay.
Pure acoustic/elastic source cells and the folded sealed-pressure path are
rejected. Existing effective-skeleton attenuation/nonlinear dispatch remains;
these checks do not qualify every material law.

### Reproducibility and remaining paper targets

Jobs 93629–93631 are graded baselines; 93643/93644 stabilization controls;
93655/93657/93658 CPML controls; 93671 the 4 m inviscid control.
Jobs 93648/93649/93656/93660/93661/93669/93673/93676 cover regressions, new
sources, MPI, superposition, FEM and guards. Jobs 93677–93680 and 93720–93723
isolate and verify the CPU VT fix. Invalid old FEM fixtures fail the existing
sloping-water-surface guard; the checkpoint fixture is rejected by the existing
CPML checkpoint guard. Rejected setup runs are retained and excluded from
accuracy results.

Runs use reservation `cuda-615-SC`, account `ASC26075`, `ml reset` and the
existing non-module GSL/PROJ paths. Imported meshes use direct MPI with one
CPU thread per rank; generated-mesh integration cases use the supported
launcher. Binary v02 SHA-256:

```text
f8e4404893e753c26690073117aa03e32d7eb3b1ccda91e993ff0a181bddccb5
```

Final v04 (CPU VT corrections) SHA-256:

```text
67e83b1a8675c1597dd15feb41493df4b78df925e864d549536a158dcc8034c0
```

All 94 tracked final solver-source hashes match the remote build. Its embedded
Git revision is inherited from the isolated copied tree, so exact manifests
and binary hashes are authoritative. Local-only helpers/results are in the private
validation archive; their README documents commands,
dependencies, input units, output, overwrite behavior and limitations. They
are unavailable in a clean checkout. Author archives and large meshes are
not redistributed.

Next work should address full-space/free-surface spatial, source-patch and
timestep convergence. Figure 6 requires distinct hydraulic materials across
the interface. LOHp additionally requires the prescribed 401-point line,
50 m spacing, 5 s Gabor-pulse record, and near-surface slow-wave resolution
against the supplied 0.625 m FD reference. Neither is reported as reproduced.
The shared open-pore CPML restriction remains: homogeneous reference tests
cannot establish that the acoustic hydraulic interface guard can be removed.

## Uniform Figure 5 boundary-distance control (2026-11-14)

The [boundary-distance report](../../doc/Hercules_Wolf_Figure5_Uniform_Boundary_Test_2026-11-14.pdf)
compares four completed inviscid runs: 4096 m and 8192 m cubes, each with
`poroelastic_pressure_stabilization=0.005` and `0`. Every cell is 8 m,
including source and CPML cells; there is no source patch or grading. These
are different meshes from the source-refined cases in the preceding section.
The element counts are 134,217,728 and 1,073,741,824, respectively, with one
leaf level and zero hanging nodes. All eight receiver cells are also 8 m.

The 4096 m cube includes a 256 m CPML on the sides and bottom. Its source is
at `(2049,2049,3073)` m. The 8192 m cube translates the source and receivers
by `(2048,2048,0)` m in solver coordinates, preserving physical offsets and
subcell phases. Its lateral CPML entrances move 2048 m outward and its bottom
4096 m downward; the top is unchanged. The far-domain earliest continuum
fast-wave return from a CPML entrance is 1.860 s, and the fixed unabsorbed
top's image-path return is 1.568 s, before adding source delay. Both exceed
the 1.5 s record. These geometric estimates do not prove numerical accuracy.

All four runs use the same 50 microsecond timestep, 0.2 ms station cadence,
original 700 author rate samples at 1 ms spacing, 16 Hz dominant source
frequency, 0.07 s delay, and stress/pressure moments of 1e10 N m each.
The meshing control `simulation_wave_max_freq_hz=24.271798903240597` with
`simulation_node_per_wavelength=10` selects uniform 8 m cells in this
homogeneous model; it is not a source cutoff or filter. CPML controls are:

```text
abc_type = cpml
cpml_thickness_m = 256
cpml_reflection_target = 1e-6
cpml_poly_order = 3
cpml_kappa_max = 1
cpml_alpha_max_hz = 16
cpml_mpml_cross_fraction = 0
cpml_alpha_floor_fraction = 0
cpml_pml_material_damping = none
```

The reflection target is a profile-design input, not a measured coefficient.
Viscosity is zero in every run; changing gamma changes numerical pressure
dissipation, not physical Darcy drag. Worst receiver full-record differences,
normalized by the corresponding far-domain numerical trace, are:

| Gamma | Solid velocity difference, % | Darcy flux difference, % | Pressure difference, % |
| ---: | ---: | ---: | ---: |
| 0.005 | 0.003676 | 0.000744 | 0.000755 |
| 0 | 0.003544 | 0.000669 | 0.000727 |

The late-window velocity differences are 0.003093% and 0.002805%,
respectively. Relative to the full author-reference norm, full-record
solid-velocity differences are at most 0.002755% and Darcy-flux differences
at most 0.000454%. **Moving these boundaries does not remove the late
ringing. Their reflected contribution is not the dominant discrepancy in
this uniform 8 m, 1.5 s experiment.** This does not establish accuracy of
the slow wave, the unchanged top boundary, or general CPML configurations.

| Cube edge, m | Gamma | Worst solid velocity L2 versus author, % | Worst Darcy flux L2 versus author, % |
| ---: | ---: | ---: | ---: |
| 4096 | 0.005 | 109.9279 | 201.8727 |
| 8192 | 0.005 | 109.9276 | 201.8727 |
| 4096 | 0 | 369.4692 | 2783.3000 |
| 8192 | 0 | 369.4692 | 2783.3000 |

Velocity comparisons respect the `t-dt/2` staggering and use the unchanged
author scaling 7025000000000, without fitted gains or shifts. The PDF gives
waveforms, full-record FAS against the reference, late-window FAS and the
boundary geometry. Pressure has a same-grid numerical comparison here, not
an independent analytical pressure reference. The source-history check is
separate: piecewise-linear interpolation changes its 16 Hz rate amplitude
by 0.0841%, and its 1–80 Hz FAS differs from the smooth Ricker by 0.1246% L2.
Correct source history does not ensure accurate numerical radiation.

Jobs 94065 (both near cases), 94066 (far gamma 0.005) and 94074 (far gamma 0)
completed with exit code zero using the unchanged v04 executable and hash
recorded above. Job 94074 replaces the cancelled, unstarted 94067. It was
initially submitted without a reservation; its later externally assigned
`polaris` reservation was explicitly confirmed by the maintainer. Run
snapshots retain inputs, hashes, modules, libraries, job records, monitors
and station files in the same Horizon archive. Helpers and raw data remain
local-only and unavailable in a clean checkout. No solver correction or
boundary/interface restriction change is part of this experiment.

## Inviscid propagation and source diagnostics (2026-11-14)

The [diagnostic report](../../doc/Hercules_Wolf_Figure5_Inviscid_Diagnostics_2026-11-14.pdf)
contains all three completed investigations, waveform/reference spectra,
energy and dispersion figures, input controls and provenance.

This follow-up isolates volume propagation, source/receiver cell placement,
and mesh/timestep sensitivity. All Figure 5 cases keep the 4096 m cube,
material, 256 m CPML profile and original source-rate samples above.
The additional source-free tests use separate periodic operator fixtures;
they are not Figure 5 domains or CPML tests. No production solver correction
is introduced.

### Source-free production GPU operators

Job 94083 completed 30 periodic tests, covering continuum and discrete
fast/slow modes, axial and body-diagonal propagation, 8/4 m spacing,
50/25 microsecond timesteps, and a short-wavelength numerical branch. The
fixture invokes the unchanged production stiffness, pressure traction,
coupled mass, continuity, pressure update and stabilization kernels. Its
independent 3D Fourier evolution is checked against eight-point Hex8
quadrature. The maximum field relative L2 difference between GPU output and
that discrete prediction is **3.36e-8**. Agreement with a discrete operator
does not establish agreement with the continuum equations.
The frequency table below uses a 50 microsecond timestep.

| Continuum mode | Spacing, m | Continuum frequency, Hz | Discrete frequency, Hz | Phase-speed error, % |
| --- | ---: | ---: | ---: | ---: |
| Fast, axial | 8 | 14.0247 | 13.9995 | -0.180 |
| Fast, axial | 4 | 14.0247 | 14.0184 | -0.045 |
| Fast, diagonal | 8 | 24.2915 | 24.0704 | -0.910 |
| Fast, diagonal | 4 | 24.2915 | 24.2361 | -0.228 |
| Slow, axial | 8 | 16.4996 | 14.8982 | -9.706 |
| Slow, axial | 4 | 16.4996 | 16.0911 | -2.476 |
| Slow, diagonal | 8 | 14.2891 | 13.2484 | -7.283 |
| Slow, diagonal | 4 | 14.2891 | 14.0221 | -1.869 |

The axial slow mode has a fixed 64 m wavelength. Halving the timestep changes
its predicted frequency by approximately 1e-5 Hz, whereas halving the mesh
spacing substantially reduces phase error. These are fixed-wavenumber
comparisons; the earlier fixed-16-Hz dispersion calculation is a different
cross-section of the same dispersion relation.

The collocated pressure-gradient symbol decreases to zero at axial Nyquist.
Consequently the unstabilized slow dispersion curve folds and has two axial
wavelengths at 16 Hz:

| Spacing, m | Longer wavelength, m (cells) | Shorter wavelength, m (cells) | Short-branch q/v |
| ---: | ---: | ---: | ---: |
| 8 | 58.3294 (7.2912) | 21.9124 (2.7391) | -54.6995 |
| 4 | 64.3836 (16.0959) | 9.1117 (2.2779) | -256.212 |

The continuum slow wavelength is 65.9985 m with q/v = -8.99569. The extra
short branch is a numerical solution, not a third physical compressional
Biot wave. Periodic GPU tests initialize compatible integer-wavenumber
short modes and recover their predicted frequencies. This establishes that
the production volume operators support the extra mode without a point
source or CPML. It does not by itself assign every part of a point-source
tail to that mode. Negative group velocity on the folded branch alone is
not evidence of growing energy.
For the longer-wavelength 16 Hz root, the unstabilized semidiscrete group
speeds are 695.9 m/s at 8 m and 978.8 m/s at 4 m, compared with 1056.0 m/s
in the continuum. This explains a mechanism for the delayed, dispersed
slow pulse; it is not a fitted arrival-time correction or a prediction of
the complete stabilized finite-domain waveform.

For zero viscosity and gamma zero, the periodic fixture measures the
staggered central-difference energy invariant:

```math
E_{n+1/2}=\tfrac12 V\sum_a\bigl(\rho|\mathbf v_a|^2+
2\rho_f\mathbf v_a\cdot\mathbf q_a+m|\mathbf q_a|^2\bigr)
+\tfrac12\mathbf u_{n+1}^{T}K_d\mathbf u_n
+\tfrac12\mathbf p_{n+1}^{T}C\mathbf p_n,
\qquad C=\operatorname{diag}(V/M),\quad V=h^3.
```

Here velocities are displacement increments divided by dt, and K_d is the
assembled drained stiffness. The maximum relative energy drift is
**1.44e-12**, including ten-second physical slow and short-grid modes.
With gamma 0.005, their final/initial energy fractions are 0.112064 and
1.604e-7, respectively: stabilization damps the physical slow mode as well
as the short branch. This intentional numerical dissipation is not viscosity.
The energy test excludes CPML, MPI partition/constraint exchange, sources,
interfaces, topography, attenuation and nonlinear constitutive evolution.

The private diagnostic executable adds only a local include and pre-MPI
dispatch to an isolated copy of the verified v04 source manifest. No
production kernel body is edited, and no public input option is introduced.
Its SHA-256 is:

```text
e28cf8168a63891168d7b33fbcb62b9b5749ffb5cca9da7163aaa94ec20eb2c8
```

### Cell placement and separated source loads

Job 94077 translates the source and every receiver by `(3,3,3)` m, keeping
all offsets fixed and placing the source at the center of its 8 m cell.
Both injection weights and receiver interpolation phases change, as do
boundary distances by 3 m. It is not an isolated source-weight experiment.

| Source/receiver placement | Gamma | Worst solid velocity L2 versus author, % | Worst Darcy flux L2 versus author, % |
| --- | ---: | ---: | ---: |
| Original | 0 | 369.469 | 2783.300 |
| Centered | 0 | 156.629 | 640.442 |
| Original | 0.005 | 109.928 | 201.873 |
| Centered | 0.005 | 103.038 | 137.486 |

At x3, 600 m along the positive X axis, centering reduces gamma-zero late
velocity RMS by about 63.4% and flux RMS by 77.2%, over 0.758195–1.499 s.
Neither centered result closely matches the continuum waveform. Diagonal
receiver d4 changes much less; the placement sensitivity is directional.

Job 94078 runs stress-only and pressure-only sources at the original
placement with gamma zero. The nonzero moment remains 1e10 N m. Their sum
recovers the combined trace within 2.21e-7 solid velocity, 1.26e-12 Darcy flux
and 2.84e-7 pressure relative L2, consistent with station-output precision.
This verifies superposition, not continuum accuracy. At x3, pressure-only
late RMS is 1.8224e-6 m/s for solid velocity and 1.0182e-4 m/s for flux,
compared with 8.8658e-8 and 5.2757e-6 m/s for stress-only. The pressure load
produces most of this receiver's late ringing, but both separate responses
still disagree with their continuum references.

Separate-load references use independently diagonalized radial outgoing
Biot Green functions with a smooth 16 Hz Ricker. Their combined response
differs from the sampled author histories by approximately 0.56–0.61% L2;
that baseline mismatch is retained, without fitted gain or shift. Worst
receiver solid-velocity/Darcy-flux errors are 45.49%/1180.18% for stress-only
and 422.65%/2549.95% for pressure-only. These are absolute waveform accuracy
errors; they are not the much smaller superposition errors above.

### Uniform refinement and timestep halving

Jobs 94080/94081 halve every cell from 8 m to 4 m at the same 50 microsecond
timestep, with gamma 0/0.005 respectively. The 4096 m cube contains
1,073,741,824 uniform cells, including 4 m source cells and 64 CPML cells
across the same 256 m thickness. Runtime audits confirm a single leaf level,
zero hanging nodes, unchanged material and source history, the expected
receiver positions/cell sizes, and complete finite 7500-row station records.
Physical source and receiver coordinates are fixed. The source's fractional
cell position consequently changes from 1/8 to 1/4; the periodic tests above
separately isolate the volume propagation error.

| Uniform spacing, m | Gamma | Full-record solid velocity L2, % | Full-record Darcy flux L2, % | Early-fast solid velocity L2, % |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 0 | 369.469 | 2783.300 | 40.182 |
| 4 | 0 | 289.612 | 3052.520 | 10.311 |
| 8 | 0.005 | 109.928 | 201.873 | 40.178 |
| 4 | 0.005 | 122.880 | 157.050 | 10.310 |

All entries are worst-receiver vector errors against the author histories,
with no fitted shift/gain or output filtering. The early-fast window ends
at `0.07+r/c_fast+0.16` seconds, before the slow pulse. **Refinement improves
the fast arrival, but the complete inviscid waveform remains inaccurate.**
With gamma zero, worst Darcy-flux error increases at 4 m. With gamma 0.005,
the flux error decreases, but the solid-velocity full-record error increases.
The slower, dispersed arrival and subsequent ringing remain visible in
the waveform and FAS overlays. Holding dimensionless gamma fixed also
decreases its dimensional fourth-order dissipation coefficient
`gamma*c_slow*h^3` by a factor of eight; no physical viscosity is added.

Jobs 94079/94082 independently halve the gamma-zero timestep to 25
microseconds at 8/4 m, keeping 0.2 ms output cadence. Maximum differences
from the corresponding 50 microsecond run, normalized by that run, are:

| Uniform spacing, m | Solid velocity difference, % | Darcy flux difference, % | Pressure difference, % |
| ---: | ---: | ---: | ---: |
| 8 | 0.011709 | 0.012745 | 0.010351 |
| 4 | 0.052028 | 0.050126 | 0.048814 |

These comparisons use common author times and the actual half-step velocity
timestamps. Reference errors are essentially unchanged: at 4 m and 25
microseconds they are 289.611% solid velocity and 3052.528% Darcy flux.
Thus halving the timestep does not resolve the discrepancy at either spacing.

The completed evidence identifies slow-wave dispersion and a spurious
short-wavelength volume branch, strong source/receiver cell-phase
sensitivity, and little sensitivity to timestep halving. Moving the
boundaries has negligible effect in the tested **8 m** control; no separate
enlarged-domain 4 m test was run, so its quantitative boundary bound is not
transferred to the 4 m result. Periodic energy conservation does not establish
CPML/interface stability or accurate point-source radiation. These results
motivate further investigation of the fluid/pressure spatial discretization
and its point-source coupling before claiming paper-level accuracy.

### Completed-run provenance and scope

All twelve distinct production cases across the boundary and diagnostic
matrices completed with exit code zero, as did all thirty periodic fixtures.
Every production snapshot's binary hash and four staged input/metadata files
were verified. The diagnostic matrix consists of the two 94065 baselines
and eight cases in jobs 94077–94082; the far controls add 94066 and 94074.
The runtime reservation snapshots record `cuda-615-SC` for 94065/94066/
94077/94078/94083, `polaris` for 94074/94079/94080, and `IO_testing` for
94081/94082. Observed external assignments were preserved after the
maintainer's explicit reservation-control instruction. No further jobs or
reservation changes are required for these completed comparisons.

The existing immutable v04 executable is unchanged. The separate periodic
driver is local-only, with its distinct hash and narrow scope above. Inputs,
provenance snapshots, station results, comparison tables and plots remain in
the private validation records;
they are unavailable in a clean checkout. The two reviewed PDFs and this
record preserve the measured results without redistributing author data.

This follow-up changes documentation only. It does not implement a new
discretization, alter a default or relax a guard. Homogeneous inviscid
Figure 5 and the periodic fixtures do not qualify poroelastic–acoustic or
poroelastic–elastic interfaces, FEM/VT topography, attenuation, nonlinear
elastoplastic/plastic or hysteretic response, heterogeneous hydraulics,
or the other Wolf paper figures. Their existing validation gaps remain.

## Uniform 4 m default-CPML-profile comparison (2026-11-14)

Job **94139** completed successfully on **32 nodes / 128 GPUs** under
`IO_testing`, from 14:11:12 to 14:30:56 America/Chicago, elapsed 19:44.
The maintainer explicitly selected gamma 0.005, retained 256 m CPML thickness,
corrected the allocation to 32 nodes, and subsequently authorized moving this
job from `cuda-615-SC` to `IO_testing`. The superseded eight-node job 94137
was cancelled before starting. These instructions do not authorize reservation
changes for other jobs.

This comparison uses the same 4096 m cube, uniform 4 m mesh, 16 Hz dominant
source-rate history, source/receiver coordinates, homogeneous inviscid
`full_biot` material, gamma 0.005, 50 microsecond timestep, 1.5 s duration and
0.2 ms station cadence as the earlier custom-profile job 94081. Meshing
`simulation_wave_max_freq_hz=48.543597806481195` and ppw 10 remain unchanged;
the meshing frequency does not filter the source. The physical CPML thickness
remains 256 m, or 64 cells. Only the profile overrides are removed to exercise
the verified v04 defaults:

| CPML control | Previous job 94081 | Default-profile job 94139 |
| --- | ---: | ---: |
| Reflection target | 1e-6 | 0.01 |
| Polynomial order | 3 | 2 |
| Kappa maximum | 1 | 5 |
| Alpha maximum, Hz | 16 | 0.5 |
| Alpha floor fraction | 0 | 0.1 |

Cross-axis damping remains zero and additional material damping remains
`none`. The executable and both launchers are unchanged. The completed
snapshot matches all four staged input/metadata hashes and the established
v04 binary hash. Runtime output confirms the default profile, full Biot,
128 GPU ranks, 1,073,741,824 elements at level 10 only, zero hanging nodes,
303,038,464 porous CPML cells, and all eight receivers at their expected
coordinates in 4 m porous cells. Every solid/fluid station file contains
7,500 finite samples through 1.4998 s. Both Slurm and the solver wrapper
report successful exit.

The full-record vector L2 errors below use the fixed author normalization,
common author times 0.001–1.499 s and the actual half-step v/q timestamps,
with no fitted gain, time shift or output filter. Each entry is the worst
of eight receivers:

| CPML profile | Solid velocity L2, % | Darcy flux L2, % | Early-fast solid velocity L2, % |
| --- | ---: | ---: | ---: |
| Previous custom profile | 122.880063 | 157.050316 | 10.310234 |
| Default profile | 122.880542 | 157.050324 | 10.310234 |

Maximum full-record differences between the two numerical runs, normalized
by the corresponding previous-profile trace, are **0.288496%** for solid
velocity, **0.007641%** for Darcy flux, and **0.006665%** for pressure.
Pressure here is a comparison between numerical runs, not an analytical
pressure-error estimate. Late spectra retain the same dominant ringing
near 34 Hz at axial receiver x3 and 27 Hz at diagonal receiver d4. Their
windows start at `0.19+r/c_slow` and end at 1.499 s; the spectra use a Hann
window and weighted-mean removal. Zero padding only interpolates frequency
samples and does not improve spectral resolution.

**Using the default CPML profile does not remove the delayed slow pulse or
late ringing in this 4 m experiment. Figure 5 remains unvalidated.** These
differences are not measured reflection coefficients, nor do they establish
that either profile is generally superior. The execution also changes from
64 to 128 GPU ranks, so the comparison is not strictly an isolated CPML
parameter change. The top boundary is unchanged and no farther-domain 4 m
run has been performed. No conclusion is extended to other durations,
viscosities, interfaces, topography or nonlinear behavior.

The private archive retains the raw output for job 94139. The local-only
`analyze_default_cpml_fig5.py` produces per-receiver audits, reference errors,
profile differences, late metrics and waveform/full-FAS/late-FAS PNG/PDF
figures. It is unavailable in a clean
checkout and does not modify raw data. No solver implementation or default
was changed for this comparison.

## 5 Hz uniform 4 m diagnostic and authorized submissions (2026-11-14)

The maintainer initially requested a lower-frequency diagnostic with submission
withheld, then explicitly authorized the two runs recorded below. The local-only case
`fig5_diag_5hz_h4_gamma005_cpml_defaults_256m` is prepared from completed job
94139. Every parsed `parameters.in` value and the complete `source.in` file
are unchanged. The sampled rate changes from a 16 Hz Ricker with 0.07 s delay
to a **5 Hz Ricker with 0.224 s delay**, preserving 700 samples at 1 ms, peak
rate, stress/pressure moments, source/receiver coordinates and no filtering.
This is a modified-frequency diagnostic, not the original Figure 5 experiment.

The geometry remains a **4096 m cube with uniform 4 m cells** and 256 m CPML
at the sides and bottom, using the same default profile. The inviscid full-Biot
material, gamma 0.005, 50 microsecond timestep, 0.2 ms output and 1.5 s record
are retained. Meshing `simulation_wave_max_freq_hz=48.543597806481195` and
ppw 10 remain fixed; setting the meshing frequency to 5 Hz would change the
mesh. The predicted last slow arrival is near 1.17 s; the main pulse fits
before the unabsorbed top's conservative 1.552 s return bound. The short
far-receiver quiet tail limits late-window spectral resolution. The fixed
CPML is thinner relative to wavelength at 5 Hz, so absorption still needs
measurement.

Preparation checks verify unchanged v04 solver-source hashes, baseline input
preservation, source-history integration and pressure increments against the
compiled source-history function. The prescribed-rate FAS peaks at 4.997 Hz;
its relative L2 difference from the ideal spectrum over 0.5–25 Hz is 0.0123%.
Independent continuum full-Biot Green references were regenerated for all
eight receivers, for solid velocity and Darcy flux. The same formula checks
against the original author histories within 0.62% L2 at 16 Hz, retaining the
fixed normalization without fitted gains or shifts. Both smooth and sampled
rate references are saved: their 5 Hz waveform discrepancy is approximately
1.17–1.23% because the radiated velocity also depends on the rate derivative.
The regular point-reference files explicitly exclude small causal-front
windows containing ideal distributional terms. No analytical pressure error
is asserted.

Inputs, references, audits, source-FAS and reference-only figures remain in the
private validation archive. The local helper
`prepare_5hz_fig5.py` and its usage/validity notes are unavailable in a clean
checkout. No solver, public input interface or default was modified.

After preparation, the maintainer explicitly authorized the 5 Hz case on
32 nodes under `cuda-615-SC` and the 16 Hz / 2 m case on 128 nodes without
a reservation. The previous 2 m / 40-node attempt, job 94169, failed during
GPU initialization after 18:08 with CUDA out-of-memory errors. Its inputs,
outputs and reservation history are preserved. Separate immutable submission
copies retain every scientific-input byte from the corresponding prepared
cases; only resource/provenance metadata changes.

| Job | Source central frequency | Uniform spacing | Nodes / GPUs | Reservation | State verified at 15:25 CST |
| --- | ---: | ---: | ---: | --- | --- |
| 94190 | 5 Hz | 4 m | 32 / 128 | `cuda-615-SC` | Running since 15:24:38 |
| 94192 | 16 Hz | 2 m | 128 / 512 | None | Pending, priority; start unknown |

Job 94190 uses case `fig5_diag_5hz_h4_gamma005_cpml_defaults_256m_n32`
and a 45-minute wall-time limit. Job 94192 uses
`fig5_diag_h2_gamma005_cpml_defaults_256m_n128` and a three-hour limit.
Both retain gamma 0.005, 256 m default-profile CPML, the 4096 m cube,
50 microsecond timestep, 1.5 s duration and the immutable v04 executable.
The 2 m source remains the original 16 Hz author history. Each node runs
four GPU ranks. Staged input hashes, executable/launcher hashes and libraries
after `ml reset` were verified. Slurm confirms both node counts and GPU-rank
exports, and explicitly no reservation field or submission option for 94192.
No other job's reservation was changed. **Submission and scheduler checks do
not establish runtime memory fit, successful completion or waveform accuracy.**

## Completed 5 Hz uniform 4 m comparison: job 94190 (2026-11-14)

Job **94190 completed successfully at 15:44:14 CST**, elapsed **19:36**,
on 32 nodes / 128 GPUs under `cuda-615-SC`. Slurm reports `0:0` and the
solver wrapper reports zero. Downloaded input/metadata hashes and the v04
executable hash match the prepared case. Runtime output confirms full Biot,
gamma 0.005, the requested default CPML profile and 256 m thickness,
1,073,741,824 uniform 4 m elements at level 10, no hanging nodes and
303,038,464 active porous CPML cells. The maximum reported setup GPU use is
34.105 GiB. All eight solid and fluid station records contain 7,500 finite
samples through 1.4998 s, at the expected coordinates and in the expected
material. Fluid-output identities `vr=q/porosity` and `vf=vs+vr` also check
within the respective text-output precision. Pressure records are complete
and finite; no analytical pressure error is reported.

The comparison uses the independently generated continuum Green functions
at each source frequency. The primary reference uses the actual piecewise-
linear prescribed source rate; the ideal smooth Ricker reference is retained
as a second check. Velocity and flux times are the station timestamps minus
25 microseconds, without fitted time shifts, gains or filters. Vector L2
errors exclude the documented small causal-front windows; unmasked errors
are also saved. The 16 Hz baseline below is **recomputed with the same
independent-reference convention**, so it differs slightly from the earlier
author-history error table. Each entry is the worst of eight receivers:

| Metric, relative L2 % | 16 Hz / 4 m, job 94139 | 5 Hz / 4 m, job 94190 |
| --- | ---: | ---: |
| Full-record solid velocity | 122.554779 | 18.767962 |
| Full-record Darcy flux | 156.600485 | 23.405234 |
| Early-fast solid velocity | 10.016883 | 1.730874 |
| Early-fast Darcy flux | 35.262278 | 2.470686 |

The early-fast window ends 1 ms before the continuum slow-wave causal front,
at `r/c_slow - 0.001 s`. Across individual 5 Hz receivers, full-record errors
are 10.69–18.77% for solid velocity and 13.09–23.41% for Darcy flux. Using the
smooth reference instead gives worst errors of 18.55% and 23.15%, respectively;
source-table interpolation does not explain the remaining discrepancy.

The slow-P packet retains phase and shape errors. Its main lobe is delayed
by approximately **3.4–5.8 ms**, measured against the same-sign lobe of the
smooth reference; no trace alignment is applied. At x3 the delay is 3.4 ms
for both fields, and at d4 it is 5.6 ms for solid velocity and 5.8 ms for
Darcy flux. The corresponding selected-lobe amplitudes are about 3–7.5% below
the reference, while the following opposite lobe can be larger. Comparing
the largest absolute peaks without matching their signs would incorrectly
identify a lobe swap as a much larger delay.

Full-record radial FAS over 0.5–15 Hz agree more closely than the waveforms:
the amplitude-spectrum L2 errors range from 5.23–9.74% for solid velocity and
0.69–1.57% for Darcy flux. FAS agreement alone does not establish phase accuracy.
Spectra use Tukey alpha 0.1; no waveform filter or amplitude normalization is
applied. Eightfold zero padding only interpolates frequency bins.

Late oscillations are greatly reduced. Define the late window from
`t0 + 4/(pi*f0) + r/c_slow` to the final physical velocity time, 1.499775 s.
At 5 Hz, across all receivers the late RMS is **0.033–0.098% of the smooth
reference's full-record peak for solid velocity**, and **0.0013–0.0024% for
Darcy flux**. Maximum late absolute amplitudes are below 0.160% and 0.0049%,
respectively. Under the same dimensionless window rule, the 16 Hz baseline
has late RMS ranges of 12.5–21.7% and 13.8–23.1% of the corresponding peaks.
Small residuals remain: x3 has a weak Darcy-flux feature near 38 Hz. The
far-receiver 5 Hz late windows are only 0.074–0.078 s long (about 13 Hz nominal
frequency resolution), so precise spectral peaks there are not meaningful.
Late spectra use a Hann window and weighted-mean removal on contiguous data.

**Lowering the source frequency strongly improves the 4 m result, but does
not establish quantitative convergence of the slow wave.** The remaining
packet delay and shape errors are consistent with accumulated numerical
dispersion; this run does not isolate their cause. The same physical CPML
thickness is smaller relative to wavelength at 5 Hz, and no farther-domain
5 Hz control has been run. The top cannot return within 1.5 s, but earlier
CPML-interface reflections are not excluded by that fact. No conclusion is
extended to longer records, finite viscosity, coupled interfaces, topography,
other formulations or nonlinear materials.

The ignored local helper `analyze_5hz_fig5.py` writes input/runtime audits,
per-receiver errors, same-sign slow-peak metrics, full/late spectra, reference
discrepancies and seven PNG/PDF comparison figures.
It refuses unfinished/failed inputs, preserves raw data, and performs no
submission. The private archive retains the raw output for each case in job 94190.
This comparison changes no solver code,
input defaults or other job reservations.

## Completed 16 Hz uniform 2 m comparison: job 94192 (2026-11-14)

Job **94192 completed successfully at 16:19:02 CST**, elapsed **42:51**,
on **128 nodes / 512 GPUs with no reservation**. Both Slurm and the solver
wrapper report successful exit. Unlike the earlier 40-node attempt, this
allocation fits: maximum reported setup GPU use is **67.472 GiB per GPU**.
Runtime and downloaded-input audits confirm the immutable v04 executable,
8,589,934,592 uniform 2 m elements at level 11, no hanging nodes,
2,424,307,712 active porous CPML cells, original physical source/receiver
coordinates, unchanged 16 Hz source history, full Biot, gamma 0.005,
50 microsecond timestep and 256 m CPML with its default profile. All eight
solid/fluid station pairs have 7,500 complete finite samples. Fluid-output
velocity/flux identities also pass. No pressure-reference error is claimed.

The following comparison retains the original Wolf 16 Hz reference and its
fixed normalization. Actual author samples are **3.125 ms apart**, spanning
0.003125–1.496875 s after the established selection inside 0.001–1.499 s.
Hercules v/q output is sampled every 0.2 ms and interpolated from its physical
time `station timestamp - dt/2` to those author times. The prescribed source
rate's 1 ms sampling is a separate quantity. No fitted gains, shifts or
waveform filters are used. Each entry is the worst vector L2 error across
eight receivers:

| Metric, relative L2 % | 4 m, job 94139 | 2 m, job 94192 |
| --- | ---: | ---: |
| Full-record solid velocity | 122.880542 | 100.010951 |
| Full-record Darcy flux | 157.050324 | 124.971321 |
| Early-fast solid velocity | 10.310234 | 2.318893 |
| Early-fast Darcy flux | 10.320026 | 2.322726 |

The early-fast window uses the previous author-comparison definition,
ending at `0.07+r/c_fast+0.16 s`. These values are directly comparable with
the previous default-CPML author-reference table, rather than the independent
continuum-reference convention used in the 5 Hz comparison. The fast-motion
error decreases by about 4.45 times. The full-record solid and flux errors
decrease by only about 19% and 20%, respectively, and remain large.

The 2 m slow packet is sharper and its extended oscillatory tail is much
smaller, but the main pulse still has substantial timing and shape errors.
For late windows starting at `0.19+r/c_slow`, radial RMS values are:

| Receiver / quantity | 4 m late RMS, m/s | 2 m late RMS, m/s | Reduction factor |
| --- | ---: | ---: | ---: |
| x3 solid velocity | 1.94080e-7 | 1.32182e-9 | 146.8 |
| x3 Darcy flux | 1.97696e-6 | 1.22779e-8 | 161.0 |
| d4 solid velocity | 2.57124e-7 | 1.25183e-8 | 20.5 |
| d4 Darcy flux | 2.48418e-6 | 1.18153e-7 | 21.0 |

Across receivers, the 2 m late RMS is 0.060–0.953% of the reference's full
peak for solid velocity and 0.062–1.000% for Darcy flux. Residual oscillations
remain; the larger d4 residual has a spectral feature near 48 Hz in this
window. Full-record radial Darcy-flux FAS errors over 1–48 Hz are only
2.85–3.92%, despite large waveform errors. This demonstrates why agreement
of amplitude spectra alone cannot establish waveform/phase accuracy. Solid
velocity FAS includes fast/slow interference and has larger errors,
37.18–48.47%. Full spectra use Tukey alpha 0.1; late spectra use Hann with
weighted-mean removal. Zero padding does not improve spectral resolution.

**Refinement improves the fast wave and strongly reduces late ringing, but
the 16 Hz / 2 m case still does not reproduce the slow waveform accurately.**
This is not evidence of an asymptotic full-waveform convergence rate. The
physical source stays fixed while its fractional cell position changes from
one quarter to one half; receiver interpolation also changes. Fixed gamma
means the dimensional fourth-order stabilization coefficient decreases
eightfold with `h^3`. CPML thickness/profile are fixed but its cell count
doubles, and GPU rank count/hardware placement change. No farther-domain
2 m control isolates boundary effects. No conclusion is extended to other
materials, interfaces, topography, formulations or nonlinear behavior.

The ignored local-only `analyze_h2_fig5.py` preserves raw data and writes
audits, per-receiver author-reference errors, full/late spectral metrics and
five PNG/PDF comparison figures. Its optional same-sign
extrema diagnostics are explicitly not group-delay estimates: with distorted
multi-lobe packets the largest same-sign peak may select a different lobe.
The private archive retains the raw output for job 94192.
This analysis changes no solver code or job reservations and submits no jobs.

## Mixed-element and mass-treatment investigation (2026-11-14)

This follow-up investigates whether the regular-brick full-Biot field pairing
contributes to the Figure 5 slow-wave errors. It changes **no production solver
code, supported input options, defaults, or benchmark inputs**. All alternative
discretizations below are independent local prototypes. No new Horizon jobs
were submitted and no reservations changed.

The current regular-brick operators use continuous nodal Q1 solid displacement
`u`, relative-fluid displacement `w = phi*(U_f-u)`, and pressure `p`, with lumped
coupled inertia and pressure storage. Pressure traction and continuity use the
same eight-node shapes in `quake/forward/kernel.cu`; the fourth-order pressure
stabilization is separate numerical dissipation. Native VT tetrahedra are not
Q1 bricks and are outside this investigation.

### Method and independent checks

The tests retain the original Figure 5 material, tortuosity 2, finite storage,
and both solid and fluid inertias. They set viscosity and numerical gamma to
zero and exclude time-discretization, sources, boundaries, interfaces and CPML.
The continuum slow and fast speeds are 1055.976329 and 3590.329282 m/s.

Eight axial discretizations were assembled by exact element quadrature and
Bloch-periodic vertex/midpoint identification. For nonzero-frequency modes,
pressure elimination gives the positive-semidefinite displacement stiffness
`K_d + B* C^-1 B`, where `B = [alpha B_u, B_w]`, `C` is pressure storage and
`*` denotes conjugate transpose. The coupled mass includes the actual cross
integral between the solid and fluid bases. This spectral elimination is not
a proposed production update; arbitrary independent static pressure offsets
require the separate pressure-coupling test.

Additional checks comprise seven 3-D corner-Nyquist pressure patterns, 324
pressure-coupling samples along four directions, and 224 full vector 3-D Bloch
systems: seven schemes, 4/2 m cells, continuum frequencies 5/16/32/48 Hz, and
axial, face-diagonal, body-diagonal and general oblique directions. The 3-D
candidate RT0 field is linear normal to each face and constant tangentially;
the corresponding axial reduction is a continuous P1 fluid field. P0 pressure
is constant in each cell. No 3-D RT0 mass-lumping rule is inferred from the
separate axial lumped comparison.

Numerical verification:

- Current Q1 axial eigenvalues agree with the previously GPU-verified analytic
  production symbol to 1.55e-16 relative.
- Separate six-element periodic assemblies agree with the complete Bloch
  spectra of all eight schemes to 7.32e-16 relative.
- Five- versus seven-point axial quadrature differs by at most 4.13e-16 relative
  in stiffness; generalized-eigenproblem residuals are below 9.31e-16 relative.
- Independent analytic Q1/Q1 and Q1/P0 3-D divergence formulas agree with
  quadrature to 2.23e-16 absolute on unit bricks.
- Full 3-D eigenproblems reproduce the axial reductions within 2.58e-11
  relative, including the richer spaces and their neutral modes.

These checks validate the independent calculations against each other and the
current volume symbol. They do not validate a candidate Hercules implementation.

### Exact pressure-coupling defect in the current pairing

For a regular uniform periodic grid, the current pressure-gradient symbol is

```math
g_i = \frac{\sin(k_i h)}{h}
      \prod_{j\ne i}\frac{2+\cos(k_jh)}{3}.
```

Each nonconstant corner-Nyquist pattern, with every `k_i h` equal to zero or
pi and at least one equal to pi, has `g = 0`. All **seven** such alternating
pressure patterns give normalized Q1 pressure-to-vector coupling below
1.13e-16 in independent 3-D quadrature. Both the solid and fluid Q1 fields are
blind to these pressure patterns. This is an exact spatial coupling nullspace,
not merely an inference from the finite-domain waveforms.

Finite storage and fluid inertia mean this finding alone is **not proof of an
exponentially growing dynamic instability**. It explains the zero-gradient
endpoint underlying the already verified folded slow branch. Gamma damps
pressure modes; it does not change the approximation spaces or establish
parameter-uniform mixed stability. The present test does not assign every
part of the observed source-waveform discrepancy to this nullspace.

### Mass and pressure storage are distinct contributors

Keeping the same Q1 fields and varying only inertia/storage integration gives
the following **fixed-frequency, axial** results at 4 m and 16 Hz. Percentages
are signed errors relative to continuum phase/group speed, not waveform L2:

| Coupled inertia | Pressure storage | Slow phase-speed error, % | Slow group-speed error, % |
| --- | --- | ---: | ---: |
| Lumped, current | Lumped, current | -2.446942 | -7.308379 |
| Lumped | Consistent | -1.218247 | -3.664182 |
| Consistent | Lumped | -1.169323 | -3.521288 |
| Consistent | Consistent | +0.012652 | +0.014328 |

Thus the current physical slow-wave dispersion is not attributable to linear
interpolation alone. Integrating both masses consistently greatly improves
this physical branch for this material/frequency. However, **all four Q1
variants retain the pressure nullspace and folded slow branch**. Consistent
mass changes the required inversions and cannot simply replace the current
local nodal update. Small speed errors at one frequency do not establish
broadband waveform accuracy, general material robustness or GPU efficiency.

### Comparison of candidate field spaces

The following table uses the same 4 m macro-element width and fixed 16 Hz
axial physical root. Except for the explicitly lumped rows, inertia and storage
are consistent. Field order is `u / w / p`.

| Discretization | Phase-speed error, % | Group-speed error, % | Additional observation |
| --- | ---: | ---: | --- |
| Current Q1/Q1/Q1, lumped | -2.446942 | -7.308379 | Two slow-branch roots at 16 Hz |
| Q1/Q1/Q1, consistent | +0.012652 | +0.014328 | Fold remains |
| Q2/Q1/Q1 | -0.012185 | -0.061350 | Fold remains through Q1 fluid field |
| Q1/Q2/Q1 | +0.048138 | +0.190239 | One additional axial neutral fluid mode |
| Q2/Q2/Q1 | +0.023683 | +0.116424 | One additional axial neutral fluid mode |
| Q1/RT0/P0 | +0.598085 | +1.796196 | No axial fold; solid-pressure issue in 3-D |
| Q2/RT0/P0 | +0.592714 | +1.780009 | Passes the tested corner-coupling screens |
| Q1/RT0/P0 axial analogue, lumped | -0.610552 | -1.829415 | Axial comparison only |

Increasing only the solid-displacement order does not fix the fluid-pressure
pairing. Increasing the fluid order to Q2 while retaining Q1 pressure removes
the tested positive-frequency axial fold, but introduces an axial fluid mode
with zero stiffness and zero projected pressure coupling at nonzero wavenumber.
It is a neutral inviscid mode, not evidence of exponential growth. The one-
dimensional test has no physical transverse motion to explain that extra mode.
It must be assessed before adopting the richer nodal space.

RT0/P0 retains fluid-pressure coupling for all seven tested 3-D pressure patterns.
However, **Q1 solid displacement with P0 pressure** loses solid-pressure coupling
for the four patterns alternating in two or three directions. This matters
when relative flow is suppressed: improving the fluid space alone does not
establish stability in the sealed/undrained limit. Q2 solid displacement removes
these particular zeros. Sampled positive coupling is not a proof on arbitrary
meshes or across all parameter limits.

The full 3-D systems support the same qualitative physical-wave improvement.
At fixed continuum wavenumber for 16 Hz and 4 m cells, the slow phase-speed
errors across the four directions range from -2.3381% to -2.3294% for the
current scheme, +0.0033% to +0.0127% for consistent Q1/Q1/Q1, and +0.1794% to
+0.5998% for Q2/RT0/P0. These fixed-wavenumber errors are a different cross-
section of dispersion than the fixed-frequency table above.

Equal macro-element size is not equal cost. On a periodic regular brick mesh,
the current scheme has approximately 7 scalar unknowns per cell; Q2/Q2/Q1 has
49, and Q2/RT0/P0 has 28, before boundaries and implementation-specific storage.
These counts are not GPU memory estimates. Face-bubble enrichment of the solid
space may be a less expensive alternative to full Q2, but was not tested here.

### Recommended next implementation study and required coverage

The most useful mixed candidate from these screens is **enriched solid
displacement, face-normal fluid flux, and cell pressure**, with Q2/RT0/P0 as
the concrete tested starting point. Consistent-mass Q1 is also a valuable
control for isolating physical dispersion, but it leaves the gradient nullspace.
Neither is approved here as a production replacement. Before a large-domain
Wolf rerun, a candidate needs finite-domain pulse, energy and fluid-content
conservation tests, viscosity/permeability and sealed-limit sweeps, and an
explicit assessment of accurate, affordable mass inversion and time integration.

Any later implementation must account for all applicable paths:

| Area | Required assessment before claiming support |
| --- | --- |
| Elastic, viscoelastic, nonlinear solid | Richer solid shapes change stiffness, mass, stress/strain sampling and constitutive integration; preserve plastic/hysteretic state and standalone limits. |
| Reduced poroelastic formulation | Its different unknowns and diffusion operator need their own compatibility assessment; full-Biot results do not validate `brick_up`. |
| Acoustic standalone and poroelastic-acoustic | Preserve the acoustic solver; derive pressure traces and normal flux/traction exchange for both open pores and no flow, including projection and energy balance. |
| Acoustic-elastic and poroelastic-elastic | Reassess interface interpolation, continuity and conservation wherever enriched solid or face-fluid unknowns meet existing solid nodes. |
| FEM and VT topography | Whole-brick enrichment, cut-tetrahedral spaces, face orientation, exposed faces and hanging constraints require separate treatment. |
| CPML and dashpots | Re-derive discrete derivative/memory and coupled-inertia placement for the chosen spaces; existing absorption tests do not transfer automatically. |
| GPU/MPI and outputs | Define shared-face ownership/orientation and both communication backends; preserve device-resident exchanges, source normalization and station pressure/flux conventions. |

Relevant primary literature is [Lotfian and Sivaselvan's fully dynamic mixed
method](https://arxiv.org/html/1506.06785v1), for incompressible constituents,
and [Rodrigo et al.'s stabilization analysis](https://doi.org/10.1016/j.cma.2018.07.003),
for consolidation and parameter robustness. The latter cautions against treating
unenriched P1/RT0/P0 as automatically robust. Their results motivate the present
tests; they do not validate these compressible brick prototypes by analogy.

The ignored local-only helper is
`investigate.py`, with usage and limitations
in the adjacent README. Its `results/` contains JSON provenance/audits, seven
CSV tables, two PNG/PDF figures, and `mixed-biot-investigation.pdf`. These helper
files and generated artifacts are unavailable in a clean checkout. The existing
Figure 5 reference results and their validation gaps remain unchanged.

## Consistent Q1 inertia/storage GPU trial (2026-11-14)

The maintainer requested trying consistent inertia/storage while retaining Q1
fields, then asked whether Q2 would be easier once the nonlocal mass inversion
was considered. This trial keeps all production code and input options unchanged.
It adds a private diagnostic entry in an isolated v04 build, reuses the existing
Q1 stiffness, pressure traction, coupled update, continuity and pressure-update
kernels, and applies consistent mass inverses only in that periodic fixture.
It does not add a supported `parameters.in` option.

### Completed job and numerical method

Job **94305** ran on one Horizon node with four independent GPU workers under
`cuda-615-SC`. It completed at **17:41:45 CST**, elapsed **1:16**, Slurm exit
`0:0`, with **all 38 cases successful**. No existing job reservations changed.
The isolated executable SHA-256 is
`654f3c49306e4afe2ad81ed5c84bb7d9af948a816ce2ea9b075a0dcc334e42c1`.
The build used `ml reset`, default modules and the required existing GSL/PROJ
library paths. Inputs, source/binary identities and job metadata were archived.

For homogeneous regular periodic Q1 bricks, the normalized consistent scalar
mass factors as `S = T_x T_y T_z`, with each one-dimensional factor having
weights `(1,4,1)/6`. The experiment applies an analytical periodic inverse
along each axis. Axes up to 49 nodes use all inverse coefficients; longer
axes truncate only the exponentially decaying inverse tail at radius 24.
This specialized inverse is **not a general mesh or MPI mass solver**.

For consistent inertia, `S^-1` acts on the solid and fluid force vectors before
the unchanged local density-block update. For consistent storage, it acts on
the assembled pressure increment. Constant Darcy drag uses the same spatial
mass as inertia, so its existing exact local relaxation commutes with this
homogeneous inverse. That property does not extend to arbitrary materials.
For nonzero gamma, the middle inverse in `K C^-1 K` and the outer pressure
update both use consistent storage. Therefore the gamma-zero comparison is
the cleanest isolation of inertia/storage effects; fixed nonzero gamma does
not imply unchanged numerical dissipation after changing `C`.

The cases include axial/diagonal fast and slow modes at 8/4 m, shear modes,
viscous checks with `eta = 0.001 Pa s`, separate inertia/storage controls,
four ten-second physical/short-wave energy tests, and eight 5/16 Hz pulses.
Every case retains a 50 microsecond timestep. Verification gives:

| Check | Maximum error |
| --- | ---: |
| Broadband vector mass-inverse residual against independent 27-point mass application | 1.897e-14 relative |
| GPU field histories/profiles versus independent discrete evolution | 3.435e-7 relative L2 |
| Lossless staggered energy drift across all gamma-zero, inviscid cases | 6.988e-13 relative |

The independent mode reference evolves all 13 states, including initial
pressure, drag and stabilization. Pulse profiles use a separate five-state
axial Fourier evolution at every represented wavenumber. These checks avoid
assuming pressure is determined solely by displacement when initial content
is nonzero. Viscous results are implementation checks against that discrete
reference, not a new viscous continuum validation.

### Controlled pulse accuracy

The pulse tests prescribe smooth **initial slow-wave data**, with solid
velocity having a 5 Hz or 16 Hz Ricker shape. They are not point-source tests
and do not use the Wolf author station histories. The periodic axial length
is 4096 m, h=4 m, record length 1.5 s; transverse fields are constant. The
continuum slow polarization and pressure are initialized together. Velocities
are compared at `t-dt/2`; displacement and pressure at `t`. No fitted gain,
time shift or waveform filtering is used.

At the receiver 1000 m from the initial packet center, full-record relative
L2 errors are:

| Frequency, Hz | Gamma | Solid velocity, lumped Q1, % | Solid velocity, consistent Q1, % | Darcy flux, lumped Q1, % | Darcy flux, consistent Q1, % |
| --- | --- | ---: | ---: | ---: | ---: |
| 5 | 0 | 15.149215 | 0.404862 | 15.214417 | 0.137867 |
| 16 | 0 | 140.784459 | 5.257166 | 144.504006 | 3.660674 |
| 5 | 0.005 | 15.141825 | 0.425481 | 15.207424 | 0.143892 |
| 16 | 0.005 | 138.649050 | 7.679462 | 142.058305 | 5.554501 |

Consistent inertia/storage greatly improves these traveling pulses with the
same Q1 fields. Nonzero gamma adds attenuation and raises the consistent 16 Hz
pulse error relative to gamma zero. This does not justify removing stabilization
for point sources: **the Q1 pressure nullspace and folded short branch remain**.
The ten-second short-mode test confirms that the mode still exists and remains
energy-conserving when gamma is zero. With gamma 0.005 its final energy fraction
is 2.994e-7; the tested physical slow mode retains 0.733941 of its initial energy.
No energy growth is observed in these cases. These conclusions do not remove
any CPML or interface restrictions.

### Effort and cost compared with Q2

Consistent Q1 is the narrower first implementation study because it preserves
the existing unknowns, connectivity and interpolation. A general implementation
still needs matrix-free or sparse mass inversion, convergence controls,
distributed exchanges and constrained updates at boundaries/interfaces. Wave
time stepping can remain central difference; the mass inversion becomes nonlocal.

Q2 requires additional DOFs, quadrature, hanging constraints, GPU kernels,
source/station interpolation and interface handling, even if higher order is
confined to poroelastic regions. The preceding Q2 prototypes **also used
consistent mass**, so their accuracy results do not demonstrate an inexpensive
local update. A higher-order diagonal-mass method is possible with suitable
quadrature, but must be designed and assessed for the complete Biot pairing;
see [Teukolsky's distinction between exact and quadrature-lumped mass](https://doi.org/10.1016/j.jcp.2014.12.012).

In this small pulse fixture, consistent-Q1 loop time was about **1.8 times**
the lumped time at gamma zero and **1.6 times** at gamma 0.005. These timings
include diagnostic reductions, copies and output, use a specialized tensor
inverse, and run each case on one GPU. They are not production throughput or
MPI scaling predictions and cannot establish whether a general consistent-Q1
implementation would be cheaper than a mass-lumped higher-order method.

The results support continuing the Q1 route before a Q2 migration. A production
decision still needs a representative partitioned mass solve, constrained
interface/boundary derivation and an actual point-source comparison. No full
Figure 5 reproduction, general material/topography support, CPML validation,
coupled acoustic/elastic validation or nonlinear validation is claimed here.

The investigation archive is local-only. Its README documents preparation, the
isolated build, private input format,
batch execution and `analyze.py JOB` . The analysis for job 94305 contains audited
CSV/JSON results, pulse/energy PNG/PDF figures
and a three-page report, `consistent-q1-comparison.pdf`.

## Production consistent full-Biot mass and storage (2026-11-14)

The production solver now has an opt-in `poroelastic_mass_matrix=consistent`
path for `full_biot`. `lumped` remains the default. This is an implementation
and operator/integration validation record; the full Figure 5 waveform
comparison with this path is still pending at the time of this entry.

### Discretization and supported scope

- Physical porous Hex8 cells use exact Q1 inertia and pressure-storage mass
  integrals. Native VT cells use P1 tetrahedral integrals with the existing
  retained-volume fractions. No Q2 fields or new mesh nodes are introduced.
- The coupled inertia includes solid/fluid cross mass. Matrix-free projected
  PCG uses the actual field-specific prolongation and its transpose, owned
  independent unknowns, sealed fluid modes, VT tangent modes, and outer clamps.
  CPU-staged and CUDA-aware MPI use the established exchange paths; ranks
  without local porous elements still participate in the collective solves.
- The consistent storage inverse is used both for the final pressure update
  and inside `K_a C^-1 K_a` stabilization. Source quadrature, source rate and
  source normalization are unchanged.
- **CPML remains hybrid:** porous CPML cells retain their lumped inertia,
  storage and paired history discretization. Nonporous contributions remain
  lumped. Physical/CPML shared nodes receive both contributions. This is not
  a conversion of CPML histories to consistent integration.
- Darcy drag and physical porous-cell solid mass damping/dashpots retain
  lumped spatial operators and enter a joint two-stage L-stable SDIRK2 solve.
  Its fixed coefficient is `1-1/sqrt(2)`, unrelated to pressure-stabilization
  gamma. With neither damping term active, the update is central difference.
  CPML polynomial-inertia damping and nonporous solid damping retain their
  existing integration; leading Darcy drag also uses SDIRK2 inside CPML. Physical
  stiffness-proportional damping remains explicit with an additional timestep
  check. This prevents applying the old local explicit physical mass-damping
  update to a consistent inertia with smaller eigenvalues.
- The earlier homogeneous periodic consistent-Q1 experiment used matching
  consistent Darcy weighting in its simplification. Its viscous modal values
  therefore are not direct references for this production lumped-drag path.
- Flat/FEM geometry, native VT, porous–elastic no-flow and porous–acoustic
  no-flow/open-pore paths are wired. Existing full-Biot/CPML restrictions remain,
  including native VT acoustic `face_mass`, field-specific open-pore CPML,
  Cartesian CPML buffers, and exclusions for checkpoint, groundwater cuts,
  geostatic, fixed-base and DRM histories. `brick_up` rejects the new option.

Controls are optional, with case-insensitive `lumped | consistent` selection.
`poroelastic_mass_tolerance=1e-10` is the double-precision default (`1e-5` for
single precision), measuring relative preconditioned residual.
`poroelastic_mass_max_iterations=200` is the default per solve. Every accepted
solve recomputes its true residual. Nonconvergence aborts; it never silently
accepts an approximate update beyond tolerance. Startup and final diagnostics
report the matrix scope, timestep ceiling, iteration counts and residuals.

### Production build and independent checks

The final production executable is `psolve-consistent-v04`, SHA-256:

```text
5f3f6b296735c83fe309bff8502cd904063897161f700ac63292bcb748cad476
```

It was built on Horizon with `ml reset`, the default compiler/CUDA modules and
the existing GSL/PROJ paths. Linked libraries were checked. All 99 relevant
solver source/header/make files match the local source manifest. Private audit
entry points exist only in separately labelled executables; they are absent
from this production binary. Inputs, binary hashes, module/library records,
logs and station outputs are archived privately by job.

Independent CPU integration uses 3x3x3 Gauss quadrature for Hex8 mass and four
exact degree-two points for tetrahedral mass. It assembles through the host
constraints/MPI, compares with the production GPU operator, then supplies that
independently assembled RHS to the actual PCG solve. Tests cover `M`,
`M+gamma*dt*D`, and pressure capacity `C`, including heterogeneous material,
hanging constraints, native VT, CPML remainder and projected boundary modes.
Final audit jobs **94373/94374**, on one/two nodes respectively, gave maximum
matrix-product relative L2 error **3.142e-16** and manufactured-inverse relative
L2 error **7.917e-10**. These norms differ from the stopping norm; accepted true
preconditioned residuals remain at or below the requested `1e-10`.

The time-integration audit invokes the production displacement update with
manufactured velocity `v(t)=sin(600*t)*shape` and independently integrated
forcing `M*shape*600*cos(600*t)+D*shape*sin(600*t)`. The isolated Darcy ODE
excludes mechanical dashpots. Over 0.01 s, errors at 100/200/400 steps were:

| Steps | Velocity relative L2 error | Reduction on halving dt |
|---:|---:|---:|
| 100 | 1.51574135e-4 | — |
| 200 | 3.94958050e-5 | 3.8377 |
| 400 | 1.00901770e-5 | 3.9143 |

This establishes second-order convergence for that manufactured update. It is
not a continuum poroelastic-wave benchmark or a nonlinear constitutive test.

### Integration, MPI and unchanged-default regressions

Final production jobs **94370/94371** completed 12 short cases: inviscid,
viscous, sealed and anisotropic flow; full-Biot viscoelastic/nonlinear solid
paths; FEM/VT open-pore and no-flow coupling; porous–elastic no-flow; Hex8
`point` projection; and acoustic velocity-potential coupling. Most records are
0.02 s at 25 microseconds; VT records use 0.004 s at 5 microseconds. These are
execution/finite-output checks, not complete source-pulse comparisons.
The original VT timestep was rejected by the new ceiling
`5.38179736429e-6 s`; reducing it to 5 microseconds passed.

Job **94375** repeated inviscid, viscous and point open-pore cases using host
MPI transport. Against CUDA-aware transport, the largest station-vector
relative L2 difference was **2.00e-10**. Comparisons use displacement, velocity,
pressure and fluid-motion groups separately, excluding the time column and
preserving matching undefined quantities. Two-node audits in **94374** also
exercise MPI partitioning and heterogeneous/hanging geometry.

Job **94376** repeated seven default-path cases against the pre-change Wolf `psolve-v04`
baseline binary (job **94330**): elastic, viscoelastic, acoustic,
acoustic–elastic, reduced porous nonlinear/drained, and lumped full Biot.
The largest physical-output relative L2 difference was **9.05e-14**, consistent
with atomic/MPI roundoff. The invalid-formulation and one-iteration-limit
cases correctly aborted. Documentation lint, TypeScript, Sites rendering and
GitHub Pages export audits pass; neither website nor `community` was published.

These checks do not certify single-precision builds, every attenuation or
nonlinear elastoplastic/hysteretic law, yielding/hardening accuracy, arbitrary
cut-cell conditioning, or long-time hybrid-CPML/interface stability. They do
not remove the Q1 pressure nullspace or the spurious short-wave branch.

### Full-size pilot and pending Figure 5 run

Successful 32-node pilots **94356/94367** used the same 4096 m cube, uniform
4 m mesh, 16 Hz central-frequency source, 50-microsecond timestep, gamma 0.005,
and default CPML profile with 256 m thickness. Each ran 50 steps (0.0025 s).
The mesh has 1,073,741,824 elements; 303,038,464 lie in porous CPML. The largest
reported GPU occupancy before stepping was approximately 40 GiB; this is a
startup snapshot, not a claim about a measured future peak.

The pilot used about 55–56 iterations per solve (maximum 58). Removing
unnecessary scalar-workspace clearing and regular-cell temporary matrix writes
reduced average solver compute time from **72.18 s to 58.46 s** for 50 steps.
Initialization is excluded from these timings. Extrapolation gives roughly
**10 hours** for 30,000 steps on 32 nodes, subject to iteration counts and
machine load. The subsequent physical-damping correction is inactive in this
undamped physical-domain case; its CPML contributions retain the existing scheme.

Full Figure 5 job **94379** uses final production v04, **32 nodes / 128 GPUs**,
reservation **`cuda-615-SC`**, and a 16-hour time limit. It retains the pilot's
physical inputs and runs to 1.5 s. Its waveform/FAS errors, source-generated
short-wave content, and late CPML behavior must be measured after completion.
No improvement over the author reference is claimed from the pilots alone.

## Consistent-mass performance controls (2026-11-14)

This follow-up implements the four requested performance investigations while
preserving Q1 fields, the consistent mass/storage equations, the hybrid CPML
mass treatment, and the existing damping integration. It adds four independently
selectable controls; their defaults retain the previous algorithm:

```text
poroelastic_mass_preconditioner = element_inverse
poroelastic_mass_initial_guess = previous
poroelastic_mass_krylov = single_reduction
poroelastic_mass_operator = auto
poroelastic_mass_tolerance = 1e-10
```

The alternative/default values are `lumped`, `zero`, `standard`, and `element`,
respectively. These controls are inactive with `poroelastic_mass_matrix=lumped`.
The tolerance default remains `1e-10` in double precision. The initial performance/tolerance study used production
`psolve-performance-v01` (SHA-256
`ea6d17583c155ab1ddfb5f04a9342e2865364c660ed98c755ad6c71376d3eb65`)
and private audit v01 (`cd84845c7db047425e7c1d3c72078d72c0fed03f186c4b38b912c948e2cb31de`).
The final revision changes only uniform-row memory layout and scalar/vector
kernel specialization. Production `psolve-performance-v03` has SHA-256:

```text
399c2575c35bca19c8036bfffa3c7b1c23a30f61e5aceb73219bbdf85aeb8064
```

The final private audit executable has SHA-256
`2e3f5b4c235119b441752325660dc0debd7bdab025c7ca19ef53a522061b49c1`.
All 99 solver source/header/make files were matched to the isolated Horizon
build manifest for each measured revision; matching source archives are retained
with the ignored local investigation records. Builds used `ml reset` and the
required existing GSL/PROJ paths; linked libraries were checked. Independent
diagnostics are absent from production binaries. The remote build source and
binary were restored and verified against production v03 after the audit build.

### Algorithms and scope

The element-inverse preconditioner is the symmetric operator
`B=L^-1 H L^-1`, with `H=sum_e L_e M_e^-1 L_e + R + gamma*dt*D`.
Element row masses `L_e`, exact Hex8 or native P1 tetrahedral inverses, the
lumped CPML/nonporous remainder `R`, and active implicit damping are assembled
through the same constraints and MPI exchanges as the original mass. It adds
one operator application per preconditioner application. Its purpose is to
reduce total work, not merely the reported iteration count.

Independent three-point Gauss integration and dense element inversion agree
with the Hex8 factorization to `6.67e-16` maximum absolute error for a unit cell.
On an ideal homogeneous periodic Cartesian grid, the preconditioned spatial
condition number changes from 27 to `64/27 = 2.37037`. This is a mathematical
uniform-grid result, not a bound for arbitrary material, interfaces or VT cuts.
A small heterogeneous assembled algebra check also remained positive definite.

For a regular Hex8 cell, let `T+ = [[2,1],[1,2]]` and
`T- = [[2,-1],[-1,2]]`. With volume `V` and the local solid/fluid material
inertia block `R_e`, the exact factors are:

```math
M_e=\frac{V}{216}T_+^{\otimes3}\otimes R_e,\qquad
L_e=\frac{V}{8}I\otimes R_e,\qquad
L_eM_e^{-1}L_e=\frac{V}{8}T_-^{\otimes3}\otimes R_e.
```

Thus the preconditioner uses the same local coupling weights with a factor
of 27 and a sign change for odd corner separation. For scalar storage,
replace `R_e` by `1/M`. For each native P1 tetrahedron the corresponding
factor is `(V/4)*(5*I-1*1^T)`. These factors are used only in the preconditioner.
The periodic one-dimensional eigenvalue product is `(4-cos(theta)^2)/3`;
the three-dimensional tensor product gives the stated `64/27` condition number.

Previous solutions are cached separately for each inertia stage, the pressure
update, and the stabilization storage inverse. The initial RHS normalization
and mandatory recomputed true residual remain in force. The
[Chronopoulos-Gear single-reduction CG recurrence](https://petsc.org/release/manualpages/KSP/KSPCGUseSingleReduction/)
combines the two inner products in one blocking collective per ordinary
iteration. Residual verification/replacement adds reductions as needed. It is
not an asynchronous pipelined implementation. Existing CUDA-aware field
exchanges are preserved; the small reduction scalars still synchronize with
the host. No claim that CUDA-aware MPI alone removes that synchronization is made.

The `auto` operator builds partial assembled 27-neighbor nodal rows only when
physical porous Hex8 spacing, density and storage modulus are globally uniform.
Each row has one writer, avoiding element atomics. Shared-node partial rows
still undergo contribution exchange. Variable material/spacing and VT select
the general element fallback. This adds setup work and stencil memory. The final kernel stores rows in
coalesced neighbor columns and specializes scalar/vector applications at compile
time; compiler reports show no stack or register spills for either kernel.

### Independent checks and integration regressions

Jobs **94389/94390** completed 13 short optimized cases covering inviscid,
viscous, sealed and anisotropic flow; heterogeneous material; FEM and native
VT; acoustic no-flow/open-pore coupling, point and face-mass projections,
velocity-potential acoustics, porous–elastic contact, and the established
viscoelastic/nonlinear-solid execution cases. These are finite-output and
integration checks, not full source-pulse validation for every combination.

The independent two-node quadrature/inverse checks in **94392**, plus the
manufactured-time fixture in **94393**, produced 15 operator audit rows. The
maximum operator-product relative L2 error was **3.316e-16**; maximum
manufactured-inverse relative L2 error was **1.438e-10**. The latter is a
different norm from the requested `1e-10` preconditioned stopping norm.
The actual damped update retained second-order convergence: velocity errors
`1.51574134e-4`, `3.94957974e-5`, and `1.00901119e-5` at 100/200/400 steps,
with reduction factors **3.8377** and **3.9143**.

Host transport **94397** agrees with CUDA-aware transport to at most
**6.02e-12** relative L2 across the tested station groups. Default-path job
**94398** passed all seven physical cases, with maximum difference
**3.36e-14** relative L2 from the previous production defaults. Its two deliberate
invalid-formulation/one-iteration cases correctly aborted; the enclosing Slurm
job therefore reports failed by design.

The short records are mostly prearrival. Their tiny absolute amplitudes can
produce large relative discrepancies: for example, a VT flux difference of
order `4e-20 m/s` exceeds the reference tail itself. Such ratios do not establish
waveform error, and are not suppressed or reported as physical accuracy.
Full-pulse tolerance comparisons below address that separate question.

### Full-pulse tolerance comparison and performance measurement

Job **94395**, using v01, completed the original tight solve and all optimizations at
`1e-10`, `1e-8`, and `1e-6`, in an identical finite 256 m cube with uniform
4 m spacing, a 16 Hz central-frequency source, gamma 0.005, 50 microsecond
steps, 32 m default-profile CPML and 0.6 s duration. This is a comparison of
algebraic error against the same discrete model, not a new full-space Wolf
reference test. It measures full/late waveforms, FAS, amplitude-weighted phase,
and independently integrated physical-domain kinetic/frame/storage energy.
The requested 0.6 s duration produces 11,999 solver steps (station records
through 0.5998 s). Energy is sampled at step midpoints and excludes CPML auxiliary energy;
comparison does not imply exact energy conservation. The late window starts
at 0.35 s. Station text precision limits interpretation of extremely small
solid/pressure differences.

All four cases completed on one node/four GPUs. Physics times exclude setup
and the independent diagnostic energy sampling. Full process wall times were
353.65, 249.79, 200.99 and 141.93 s, respectively.

| Configuration | Tolerance | Physics time, s | Speedup | Mean inertia/storage iterations |
| --- | ---: | ---: | ---: | ---: |
| Previous algorithm | 1e-10 | 347.38 | 1.00x | 47.46 / 48.39 |
| All optimizations | 1e-10 | 243.99 | 1.42x | 18.17 / 18.29 |
| All optimizations | 1e-8 | 194.81 | 1.78x | 13.32 / 13.45 |
| All optimizations | 1e-6 | 135.98 | 2.55x | 8.43 / 8.71 |

Errors below are relative L2 **fractions**, not percentages, and take the worst
station/physical-quantity group. No gain, time shift or filter is fitted.
Phase is reference-amplitude-weighted RMS over 0 < f <= 80 Hz where the
reference amplitude is at least 1% of its spectral peak.

| Optimized tolerance | Full waveform | FAS amplitude | Phase RMS, degrees | Late waveform | Total physical energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1e-10 | 2.35e-8 | 2.70e-8 | 7.36e-7 | 8.08e-7 | 1.02e-11 |
| 1e-8 | 1.63e-7 | 1.17e-7 | 1.97e-6 | 8.74e-6 | 2.96e-9 |
| 1e-6 | 8.84e-6 | 6.29e-6 | 1.51e-4 | 6.01e-4 | 2.40e-7 |

The largest separate kinetic/frame/storage energy relative L2 errors are
`4.32e-11`, `7.78e-9`, and `6.49e-7`, respectively, so agreement of total
energy is not solely cancellation between components.

The corresponding maximum late absolute differences divided by the full
reference peak are `9.81e-8`, `1.26e-7`, and `8.65e-6`. Thus the 1e-6 worst
late-window L2 error is about **0.060%**, while the 1e-8 full-waveform error is
about **0.0000163%**. These results support testing relaxed tolerance for this
inviscid configuration; they do not establish a universal tolerance for longer,
strongly heterogeneous, nonlinear or coupled-interface runs. The solver default
remains 1e-10, and all performance-control defaults retain the previous path.

### Revised uniform kernel and four-node pilot

The first four-node sweep, **94394**, completed five incremental v01
configurations: reference, assembled rows, element-inverse preconditioning,
previous guesses, and single-reduction CG. Each uses a 2048 m cube with 4 m
cells, 256 m CPML and 50 timesteps, without private diagnostics. The mesh has
134,217,728 elements. Total elements per GPU match the 32-node 4096 m case,
but the smaller domain has a different CPML fraction; this is not a direct
32-node scaling measurement. Times below are the reported rank-zero solver
timers, not an independent maximum-rank timing. Process wall time includes
solver setup but excludes the separate mesh-generation phase.

| v01 configuration (1e-10) | Physics time, s | Process wall time, s | Mean inertia/storage iterations |
| --- | ---: | ---: | ---: |
| Previous algorithm | 39.10 | 122.76 | 55.26 / 56.28 |
| Rows only | 42.65 | 127.73 | 55.26 / 56.28 |
| Rows + element inverse | 22.49 | 107.44 | 15.62 / 16.10 |
| Above + previous guesses | 21.12 | 106.59 | 14.04 / 13.52 |
| Above + single reduction | 22.35 | 108.11 | 14.04 / 13.51 |

The first row kernel was **9.1% slower** than the original operator. The v03
revision addresses its noncoalesced table access and runtime component loops;
no performance improvement is assumed from compiler output alone. Also, fewer
collectives did not make single-reduction CG faster in this v01 four-node test.

Final-kernel job **94413** passed inviscid, FEM open-pore and point-projection
integration cases; their physical station groups differed from v01 by at most
**1.69e-13** relative L2. Two-node **94414** independently checked inviscid and
open-pore operators: maximum product error **3.30e-16**, manufactured-inverse
error **6.99e-11**. Final v03 pulse job **94415** used the element-inverse
preconditioner, previous guesses, assembled rows and **standard CG**, retaining
`1e-10`. It completed all 11,999 steps, with worst waveform/FAS errors
**2.35e-8 / 2.70e-8**, phase RMS **7.36e-7 degrees**, late waveform **8.11e-7**,
and total physical energy **1.03e-11**, against the same original tight solution.
Physics/process wall times were **242.49 / 248.56 s**, with mean inertia/storage
iterations **18.168 / 18.286**. Thus the revised kernel and standard-CG
combination also preserve the pulse within the measured algebraic errors.

Final production job **94416** completed all five v03 cases on the same
four-node allocation, with no private diagnostics:

| v03 configuration (1e-10) | Physics time, s | Process wall time, s | Mean inertia/storage iterations |
| --- | ---: | ---: | ---: |
| Previous algorithm | 39.19 | 123.32 | 55.26 / 56.28 |
| Rows only | 37.00 | 122.11 | 55.26 / 56.28 |
| Rows + element inverse + previous guesses; standard CG | **17.16** | **101.93** | **14.04 / 13.49** |
| Above + single reduction | 18.27 | 103.25 | 14.04 / 13.51 |
| Element inverse + previous guesses + single reduction; element operator | 20.61 | 103.72 | 14.08 / 13.50 |

The revised rows alone reduced physics time by **5.6%**. The best combination
was **2.28x faster** in physics time (**56.2% less time**), with **17.3% less
process wall time** for this short run. Single-reduction CG was **6.5% slower**
than standard CG with the other options fixed. Retaining single reduction but
using rows instead of the element operator reduced physics time by **11.3%**.
These are single-allocation pilot measurements, not confidence intervals or
measured 32-node/full-duration speedups. The complete 0.6 s pulse separately
showed a **1.43x** physics-time speedup with the recommended tight configuration.

For the next comparable uniform-mesh experiments, the measured choice is:

```text
poroelastic_mass_matrix = consistent
poroelastic_mass_preconditioner = element_inverse
poroelastic_mass_initial_guess = previous
poroelastic_mass_krylov = standard
poroelastic_mass_operator = auto
poroelastic_mass_tolerance = 1e-10
```

Keep the original defaults for general runs until the selected case is checked.
The tolerance sweep supports a separate `1e-8` comparison; its results do not
justify silently changing the tolerance default. Single reduction remains an
option for larger process counts, where the communication tradeoff may differ.

Single precision and an optimized full-Biot drained exposed surface have not
been run in this investigation (the unchanged reduced/drained default was
checked). The tests cover the listed fixtures, not every constitutive law,
terrain or interface arrangement;
existing full-Biot and porous-CPML validation gaps remain. These optimizations
change the linear algebra, not Q1 dispersion or the pressure nullspace.

The existing full Figure 5 baseline job **94379** and its `cuda-615-SC`
reservation were left unchanged. All new optimization jobs also use
`cuda-615-SC`. Pending jobs were reordered only through dependencies to allow
available nodes to run the independent accuracy tests first.

### Consistent-mass local fusion and profiling follow-up (2026-11-14)

The `fusion-v01` production executable has SHA-256
`625a2628bdd6b2ebf0232dcf2ebc4efc9dcf9f1146aed129f8c7ea7becca55ef`.
It was built on Horizon using `ml reset`, with resolved GSL/PROJ dependencies.
An isolated audit executable adds independent quadrature, manufactured-time and
physical-domain energy diagnostics; these hooks and the investigation helpers
are local-only and unavailable in a clean checkout. A 100-file source manifest
was checked after restoring the production build tree.

`HERCULES_BIOT_MASS_FUSION=1` (default) combines neighboring local update,
preconditioner, operator-finish and dot-product operations. `0` preserves the
separate-kernel comparison path. Both PCG variants and preconditioners retain
MPI contribution assembly, constraints and true-residual acceptance. Vectors
still cover harbored nodes, and solid/relative-fluid halo messages remain
separate. This is a performance change, not a change of Q1 approximation spaces
or physical equations. Existing production binaries and jobs were not replaced.

`HERCULES_GPU_MPI_PROFILE` now distinguishes consistent-mass mechanics,
consistent-mass storage and ordinary Biot fluid traffic. The old mechanical
mass tags were omitted from its counters, and storage was mixed with ordinary
pore traffic. `HERCULES_BIOT_MASS_PROFILE=1` additionally measures local GPU
intervals with events and CPU exchange/reduction waits. It allocates events at
setup, collects them at existing reduction waits and drains the final events at
teardown. CPU/GPU intervals overlap and must not be summed. No hardware memory
bandwidth measurement was made.

Job **94475**, four nodes/16 GPUs, ran identical 2048 m cube, uniform 4 m,
50-step cases in one IO_testing allocation. Physical/source inputs, tolerance
`1e-10`, standard CG, element-inverse preconditioning, previous guesses,
assembled rows and partition weight 4.38 were fixed. Production timing excluded
both private audits and profiling:

| Kernels | Physics, s | Inertia PCG, s | Storage PCG, s | Process wall, s |
| --- | ---: | ---: | ---: | ---: |
| Separate | 17.222129 | 8.998964 | 6.461244 | 102.972346 |
| Fused | 11.847048 | 6.524004 | 3.563379 | 96.873302 |

Fusion gave **1.454x** faster physics (**31.21% less time**) and **5.92%** less
process wall time in this setup-dominated pilot. Mean inertia/storage iterations
were 14.04/13.51 for both runs. The largest relative station-group L2 difference
was **1.30e-13**. This gain is relative to the already optimized consistent-mass
algorithm, not the original lumped integrator. It is one allocation, not a
confidence interval or a measured 32-node/full-duration speedup.

Separate profiled cases in that allocation measured these cumulative
rank-average GPU intervals (seconds):

| Work | Separate | Fused |
| --- | ---: | ---: |
| Inertia nodal/vector/dot passes | 4.416585 | 2.414532 |
| Storage nodal/vector/dot passes | 3.301061 | 0.891809 |
| Inertia mass/preconditioner operator | 1.060256 | 1.063087 |
| Storage mass/preconditioner operator | 0.918719 | 0.917926 |

The local vector/dot work fell by about **2.33x** in aggregate; operator work
was almost unchanged. These intervals support prioritizing vector-pass fusion,
but cannot establish HBM bandwidth saturation or isolated network cost.
Small floating-point changes can alter the last iteration of a solve; profiled
iteration counts differed slightly. Profiled physics times were 17.382/11.934 s;
all quoted production speedups use the unprofiled pair.

Complete small-domain pulse job **94477** ran all 11,999 steps of the established
0.6 s, 256 m cube case on one node/four GPUs, with uniform 4 m cells and 32 m
CPML. This is a separate accuracy/regression fixture; the Figure 5 production
domain is a 4096 m cube. Each fusion pair uses the same preconditioner, guess,
operator, CG recurrence and tolerance. The `1e-10` pair uses **standard CG**;
the `1e-6` pair inherited **single-reduction CG** from its earlier tolerance
fixture, as verified in both archived inputs and startup logs. An earlier
version of this record incorrectly described both pairs as standard CG. Thus
each off/on fusion speedup is controlled, but the difference between the two
speedups cannot be attributed to tolerance alone. Energy audit hooks
were active for both members; timing is reported separately from the production
pilot above. Errors below are fractions, comparing fused to separate kernels,
not comparisons against the continuum Wolf solution.

| Tolerance | Physics separate / fused, s | Speedup | Worst waveform L2 | Worst FAS L2 | Worst late waveform L2 | Physical energy L2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e-10` | 272.606 / 240.939 | 1.131x | 3.41e-11 | 2.40e-11 | 2.06e-9 | Equal to printed precision |
| `1e-6` | 142.604 / 131.671 | 1.083x | 8.34e-8 | 5.92e-8 | 3.12e-6 | 5.43e-10 |

Late means `t >= 0.35 s`. Amplitude-weighted phase RMS below 80 Hz, using
reference amplitudes above 1% of their peak, was at most 2.96e-10/1.16e-6 degrees
for the tight/loose pairs. Mean inertia/storage iterations were unchanged to
three decimals within each pair: 18.168/18.286 at `1e-10`, 8.430/8.709 at `1e-6`.
The small one-node pulse gives less speedup than the much larger four-node pilot;
neither establishes a production 32-node speedup. Output precision limits the
smallest reported differences. These checks preserve the selected tolerance;
they do not establish `1e-6` accuracy for arbitrary models.

Independent two-node audit job **94485** passed **33** mass-product/manufactured-
inverse checks spanning inviscid/viscous, heterogeneous/anisotropic, sealed,
FEM-drained, VT open/no-flow, point open-pore, and both CG/preconditioner paths.
The largest quadrature-product/inverse relative L2 errors were **3.30e-16 /
7.92e-10**. Job **94473** checked both CG recurrences and preconditioners; matched
separate-kernel job **94483** kept the lumped-preconditioner waveform differences
below **1.64e-13**. Production regressions **94478** covered elastic,
viscoelastic, acoustic, acoustic–elastic, reduced poroelastic, lumped full Biot,
porous–elastic/no-flow, porous–acoustic/open-pore with velocity potential, and full-Biot viscoelastic and
nonlinear-solid paths. Compared with their stored baselines, those production
station groups differed by at most **2.27e-13**. Pressure-formulation point
open-pore checks in 94481/94485 differed by at most **1.78e-13** from their
stored baseline; the viscous fixture in 94473
differed by **5.67e-12**. This is regression coverage of the listed fixtures,
not every constitutive law, topology or long-duration stability regime.

Matched two-node separate-kernel references **94492** gave a maximum relative
L2 difference **1.89e-13** for the drained FEM case. Near-zero fields in the
anisotropic and VT fixtures are more sensitive in relative terms: the worst
anisotropic relative-fluid velocity difference was **2.40e-4**, with a maximum
absolute difference **9.94e-23 m/s** against a **2.77e-19 m/s** reference peak.
VT open/no-flow pressure differences were **5.41e-5 / 2.93e-4** relative, with
maximum absolute differences **1.65e-17 / 9.18e-17 Pa** against reference peaks
about **3.1e-13 Pa**. These are tiny prearrival signals; do not describe every
component of these short fixtures as agreeing to machine precision. The
independent constrained-operator/inverse checks and the full pulse provide
separate evidence. Host-MPI job **94481** also checked inviscid, point open-pore
and VT open-pore paths; its inviscid comparison with the current one-node
CUDA-aware result differed by **1.98e-13**. Matched host-MPI VT off/on job
**94499** completed with the optional profiler enabled; its largest pressure
difference was **1.57e-4** relative and **4.59e-17 Pa** absolute against a
**3.08e-13 Pa** reference peak. This confirms the same near-zero-signal caveat
on the host backend. Single precision remains untested.

The manufactured damping test **94486** retained error reductions **3.838** and
**3.914** on halving the timestep. Its permeability is `6e-9 m²`, viscosity
`0.001 Pa s`. An additional strongly stiff test at `6e-13 m²` failed the expected
factor-four test at 100 and 50 microseconds: relative velocity errors were
**0.0624878848** and **0.0280979052**, ratio **2.22393393**. The same errors occurred
with fusion on (94476) and off (94487). This identifies an existing temporal
accuracy limitation at the tested steps, not a fusion-induced change. The
operator/inverse checks passed separately. The broad audit job 94476 was stopped
and replaced by operator-only job 94485; do not count the stiff temporal test as
validated or infer uniform second-order accuracy for arbitrary Darcy stiffness.

The 32-node tests use the full **4096 m cube**, uniform 4 m cells, **50 steps**
at 50 microseconds, the original Figure 5 source/receivers, gamma 0.005, 256 m
CPML and `1e-10` tolerance. Job **94474** compares partition weights **4.38, 2, 1**
with separate kernels, followed by profiling at 4.38 and 1. Job **94488**, dependent
on successful completion of 94474, compares separate/fused pairs for all three
weights in one allocation. Both request **IO_testing** and were queued at this
record's preparation. No 32-node speedup or preferred partition weight is yet
established. The parser default remains 4.38. Original jobs 94379, 94436 and
94445 and their reservations were left untouched.


### Compact mass vectors and combined solid/fluid packets (2026-11-15)

The `compact-v01` production executable has SHA-256
`1bd96de1acf1113e6a3442e799476ee8e41380f5fc1bab7936e7b3fd6577ca31`.
Horizon built it with `ml reset` and resolved GSL/PROJ dependencies. A 102-file
source manifest matches the implementation. The isolated `compact-audit-v01`
adds independent quadrature, manufactured solves, dense projected nodal-block
references and physical-domain energy checks. The private hooks, scripts and
fixtures remain ignored and are unavailable in a clean checkout.

Two independent runtime controls default to `1`:

- `HERCULES_BIOT_MASS_COMPACT` stores Krylov work/history only on nodes touching
  consistent porous cells plus the complete solid/pore constraint and MPI owner
  closure. Excluded nodes have independent nodal blocks, solved directly with
  the same mass/storage and damping stages. Their RHS norm remains in the
  full-system stopping normalization. Physical fields and station layouts stay
  full-sized. Nodes touching both physical and CPML cells remain iterative.
- `HERCULES_BIOT_MASS_HALO_BATCH` combines solid and relative-fluid payloads in
  one packet per neighbor per mechanical phase. Both CUDA-aware and host MPI
  use GPU packing/constraint kernels; host MPI stages only packets. Fixed peer
  unpack order, duplicate mappings and receive-reuse events are retained.

Setting both to `0` restores the previous full-vector/separate-exchange path.
The previous fusion option is independently selectable. These changes preserve
Q1 spaces and the mass/storage equations; they do not resolve pressure nullspaces
or the previously reported stiff temporal-accuracy limitation.

Four-node/16-GPU job **94524** compared all four combinations in one IO_testing
allocation, with fusion on, standard CG, element-inverse preconditioning,
previous guesses, automatic assembled rows and tolerance `1e-10`. Mesh/source,
partition weight 4.38 and physical inputs were identical: 2048 m cube, uniform
4 m cells, 256 m CPML, 50 steps at 50 microseconds. Private audits and profiling
were disabled.

| Compact / combined | Physics, s | Inertia PCG, s | Storage PCG, s |
| --- | ---: | ---: | ---: |
| Off / off | 11.446063 | 6.042852 | 3.655958 |
| On / off | 10.400583 | 5.209176 | 3.455776 |
| Off / on | 11.282936 | 5.858638 | 3.651506 |
| On / on | 10.182923 | 4.986082 | 3.458897 |

Both changes gave **1.124x** faster physics (**11.04% less time**) relative to
already fused kernels. Compaction alone gave **1.101x**, combining packets alone
**1.014x**. These are single-allocation observations, not confidence intervals.
Process wall time fell from **97.52 to 96.02 s** (**1.54%**) because setup
dominates this short run. The compact set retained 67,684,865 of 136,888,441 harbored node copies: **49.45%**.
Thus the work/history vector allocation fell by **50.55%**, before accounting
for new maps/packets; this is not a 50.55% reduction in total GPU memory.
The largest rank's retained set was 6,960,617 nodes, versus 10,662,984 full nodes
at the largest full-vector rank. Worst station-group relative L2 difference
among these three comparisons was **9.16e-14**. The pilot ends before arrivals
at the distant receivers, so complete pulses and independent operator tests
provide the accuracy evidence separately.

Complete inviscid pulse job **94541** ran all 11,999 steps of the 0.6 s,
256 m cube fixture with 4 m cells, 32 m CPML, 50 microsecond steps and the
16 Hz source on one node/four GPUs. Both tolerances use **standard CG**,
element-inverse preconditioning, previous guesses and assembled rows. Each
pair keeps all physical inputs, tolerance and fusion fixed, changing only the
two new switches. Energy hooks are active for both members; these timings are
separate from the production pilot. Errors are fractions against the existing
full-vector implementation, not the continuum Wolf solution.

| Tolerance | Physics off / on, s | Speedup | Worst waveform L2 | Worst FAS L2 | Worst late waveform L2 | Physical energy L2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e-10` | 247.406 / 218.170 | 1.134x | 1.07e-8 | 7.58e-9 | 1.79e-7 | 8.37e-12 |
| `1e-6` | 143.153 / 118.190 | 1.211x | 8.85e-6 | 9.77e-6 | 1.39e-3 | 1.71e-7 |

Late means `t >= 0.35 s`; the largest late absolute discrepancy normalized by
the full waveform peak was **1.09e-8 / 8.68e-6** for the tight/loose pairs.
Amplitude-weighted phase RMS below 80 Hz (bins above 1% reference peak) was
**7.51e-8 / 2.61e-4 degrees**. Energy uses 59 samples outside CPML; output
precision limits very small differences. Compaction eliminates direct nodal
blocks from Krylov iteration, so its recurrence and accepted solution need not
match full-vector CG bitwise even though both solve the same equations with
the same full-RHS tolerance. In particular, the loose-tolerance late-waveform
difference is **0.139%**, larger than the full-record relative difference.
These measurements do not establish `1e-6` accuracy for arbitrary models or
justify multiplying speedups from unmatched pilots.

Independent two-node audit **94540** passed all **33** quadrature-product and
manufactured-inverse checks for inviscid/viscous, heterogeneous/anisotropic,
FEM open/drained, VT open/no-flow, point open-pore and both CG/preconditioner
paths. Largest product/inverse relative L2 errors were **3.30e-16 / 7.92e-10**.
Direct excluded-block tests compare GPU mechanics against an independent dense
projected solve (including both damping stages) and storage against direct
scalar division: maximum solution L2 error **9.05e-17**, normalization error
**3.95e-15**. These setup tests use zero old increments; pulse comparisons
separately exercise history and accumulated updates. Temporal audit **94542**
retained halving ratios **3.838 / 3.914** for the moderate-drag manufactured
fixture. The earlier strongly stiff temporal-accuracy limitation is unchanged.

Production jobs **94523/94536** completed the standalone and coupled regression
fixtures: elastic, viscoelastic, acoustic, acoustic–elastic, reduced poroelastic,
lumped full Biot, consistent porous–elastic, acoustic pressure/velocity-potential
open pores, full-Biot viscoelastic and nonlinear skeletons, FEM drained terrain,
and VT open/no-flow. Host-MPI job **94537** completed six matching porous/interface
fixtures. These are regression checks of the named configurations, not blanket
validation of every constitutive law or long-duration interface stability.

Matched full-vector references **94546** help distinguish numerical differences
in tiny signals. The largest anisotropic relative-fluid velocity difference
was **1.85e-4** relative and **7.33e-23 m/s** absolute against a **2.77e-19 m/s**
peak. CUDA-aware VT open/no-flow pressure differences were **3.60e-4 / 8.17e-5**
relative, with absolute maxima **1.56e-16 / 5.55e-17 Pa** against approximately
**3.1e-13 Pa** peaks. FEM open/drained displacement differences were at most
**1e-24 m**, and the porous–elastic fixture's largest relative difference,
**3.51e-6**, corresponded to **8.60e-20 m**. These near-zero prearrival records
should not be described as bitwise or uniformly machine-precision agreement.

Two-node host-MPI job **94545**, with overlap disabled and profiling enabled,
completed all four compact/batch combinations and an unfused comparison on the
same VT open-pore fixture. With compact maps and identical mean inertia
iterations (**17.360**), combining packets halved rank-average mechanical send
messages from **423,423 to 211,711.5**, while preserving the **1.891e9-byte**
payload. The combined full-map variant sent **1.255e10 bytes**. Profiled wall
timings are not used for production speedup claims. The largest near-zero
pressure difference for the fused compact+combined path was **3.38e-4** relative,
**1.54e-16 Pa** absolute, against a **3.08e-13 Pa** peak.

The two-node shallow porous-layer test in **94544** retained 5,780 of 39,304
node copies, with rank retained counts **0–1,445**. All ranks, including empty
ones, participated in operator/inverse checks and propagation. Main elastic
receiver displacement/velocity histories differed by at most **2.51e-15**
relative L2. At a silent prearrival receiver, the full-vector solver retained
up to **1.58e-14 m/s**, while the compact solver returned zero within the shared
full-system tolerance. That component's relative difference is **1**, despite
its tiny absolute size; it is not evidence of a resolved wave discrepancy.

The complete viscous pulse pair in **94544** used two nodes/eight GPUs, the
same 256 m/4 m/32 m CPML geometry and 0.6 s record, standard CG and `1e-10`
tolerance, with viscosity **0.001 Pa s** and permeability **6e-13 m²**. Both
members ran private setup audits; energy sampling was disabled. Physics time
fell from **672.186 to 527.203 s**, **1.275x** faster (**21.57% less time**).
Worst waveform/FAS relative L2 differences were **5.25e-8 / 3.71e-8**, phase RMS
**3.60e-7 degrees**, late waveform L2 **1.43e-5**, and late peak difference
normalized by the full-record peak **1.58e-8**. The CPML direct-block dense
references passed with maximum solution/norm errors **9.31e-17 / 3.16e-15**.
This checks both damping-stage recurrence and late CPML behavior against the
existing discrete solver; it does not resolve the separate continuum or stiff
temporal-convergence questions.

Single precision, restart equivalence and a full-duration 4096 m production
run have not been independently rerun with these two switches. The startup
maps/work vectors are reconstructed from the mesh; physical state/output
layouts and existing full-Biot/CPML restrictions are preserved. This follow-up
establishes the listed algebraic and regression comparisons, not new continuum
accuracy or stability guarantees. Both documentation build targets, their
rendered parameter/runtime audits, lint and TypeScript checks passed.

The full-size 32-node job **94547** requests IO_testing after existing pilot
94488. It compares the same four combinations on the 4096 m cube, uniform 4 m
mesh, 256 m CPML and 50 steps, keeping weight 4.38 and tolerance `1e-10` fixed.
It was queued at preparation; no 32-node speedup is claimed here. Existing jobs
94379, 94436, 94445, 94474 and 94488 and their reservations were left unchanged.

## Completed Figure 5 consistent-Q1 comparison: job 94445 (2026-11-15)

Job **94445** completed with exit `0:0` on 32 nodes/128 GPU ranks under
`IO_testing`, with Slurm elapsed **2:31:26** and rank-zero physics time
**8,839.914 s**. The case `fig5_optimized_h4_t6_io` uses the `performance-v03`
executable from revision `85f630be5be25f47445f51f548b97204310d7d9e`, SHA-256
`399c2575c35bca19c8036bfffa3c7b1c23a30f61e5aceb73219bbdf85aeb8064`.
It predates the later vector-pass fusion, compact iterative vectors and combined
solid/fluid mass halos; this result does not validate those later changes on a
full-duration Figure 5 run.

The audited model is the inviscid full-Biot case: 4096 m cube, uniform **4 m**
elements, **1,073,741,824 elements**, zero dangling nodes, **50 microsecond**
timestep, **1.5 s** duration and gamma **0.005**. Consistent Q1 inertia/storage
apply to 770,703,360 physical cells; 303,038,464 porous CPML cells remain lumped.
CPML is **256 m** thick with reflection target `0.01`, order `2`, kappa `5`,
alpha `0.5 Hz`, alpha floor fraction `0.1` and cross fraction `0`.
The source rate is the authors' **16 Hz central-frequency** Ricker, delayed
0.07 s: all 700 samples at 1 ms spacing match exactly. The mass solver uses
`1e-6`, `element_inverse`, `previous`, standard CG and the automatic operator.
All eight solid/fluid histories contain 7,500 finite samples at 0.2 ms cadence.
The fluid output identities `q = phi * (vf - vs)` and `vr = q / phi` pass.

### Reference and comparison conventions

The reference is the authors' **Karpfinger semi-analytical full-space** archive
`homogeneous-fullspace/reference/16Hz_explosive_M1i.zip` from the already
downloaded Wolf dataset, not the Gar6More3D free-surface reference. The inner
ZIP SHA-256 is
`3ab14de9b21654bd07dcc7b2870c52ac0094d34a2de65b8339b20baa0a2a4c14`.
The authors' converter supplies the fixed scale `7025000000000`; neither
amplitude nor timing is fitted to the numerical result. All eight receiver
offsets match the author coordinates. Hercules velocity and Darcy flux use
their physical staggered time, printed time minus **25 microseconds**, and
are interpolated onto the original **3.125 ms** reference grid, from
**0.003125 to 1.496875 s**. No waveform filtering is applied.

Relative L2 below includes all three vector components and is normalized by
the reference vector norm over the same time window. Plots show radial
components. An independent outgoing modal Green-function calculation with a
smooth Ricker gives fast/slow speeds **3590.329 / 1055.976 m/s** and agrees with
the archived reference within **0.593% solid velocity / 0.613% Darcy flux** L2.
This checks normalization, orientation and arrival speeds; its smooth source
is not identical to the archived sampled-source processing. The archived
histories remain the reference for every error reported below.

### Full-record errors

| Receiver | Distance (m) | Solid velocity L2 (%) | Darcy flux L2 (%) |
|---|---:|---:|---:|
| x1 | 1000.000 | 9.034 | 9.383 |
| x2 | 600.000 | 6.328 | 5.987 |
| x3 | 600.000 | 7.090 | 6.362 |
| x4 | 1000.000 | 9.888 | 9.710 |
| d1 | 995.929 | 6.791 | 7.966 |
| d2 | 597.558 | 2.427 | 2.574 |
| d3 | 597.558 | 4.778 | 5.360 |
| d4 | 995.929 | 4.593 | 5.116 |

For comparison, earlier job **94139** used the same physical inputs, source
files, uniform mesh, timestep, gamma and CPML settings, with lumped Q1
inertia/storage. All changed active input keys are `poroelastic_mass_*`.
The executable revisions differ, so this is not a same-binary switch test.

| Metric, maximum across receivers | Lumped 94139 | Consistent 94445 |
|---|---:|---:|
| Full solid velocity L2 (%) | 122.881 | 9.888 |
| Full Darcy flux L2 (%) | 157.050 | 9.710 |
| Early fast-P solid velocity L2 (%) | 10.310 | 6.604 |
| Early fast-P Darcy flux L2 (%) | 10.320 | 6.606 |
| Slow-P solid velocity L2 (%) | 152.236 | 11.249 |
| Slow-P Darcy flux L2 (%) | 157.065 | 9.711 |
| Solid radial FAS L2, 1–48 Hz (%) | 58.446 | 7.909 |
| Darcy radial FAS L2, 1–48 Hz (%) | 15.053 | 7.921 |

The early window ends at `0.07 + r/c_fast + 0.16` seconds; the slow window
starts at `0.07 + r*(1/c_fast + 1/c_slow)/2`. The large improvement is in the
slow-P waveform and subsequent ringing. It is not an improvement in every
window: early fast-P solid errors at **x3/x4** increase from **2.33/3.78%** to
**4.92/6.60%**; their flux errors increase from **2.52/3.82%** to **4.98/6.61%**.
The remaining errors have not been assigned uniquely to any numerical term.

Full-record FAS uses a Tukey window with alpha `0.1`. Amplitude-only spectral
error does not capture phase error, explaining why the earlier lumped Darcy
FAS error is much smaller than its waveform error. Separate wavelet diagnostics
give maximum component envelope/phase misfits **6.173/1.573%** for solid
velocity and **5.910/1.420%** for Darcy flux. These use Morlet `w0=6`, 80
logarithmically spaced frequencies over 1–40 Hz and reference normalization;
they are not L2 percentages or a claim of identical analysis settings to the
paper. Agreement remains short of the paper's Figure 5 accuracy.

### Late oscillations and limitations

The late window starts at `0.19 + r/c_slow` seconds. The table reports radial
RMS divided by that field's full-record reference peak, not relative error
against the nearly silent late reference.

| Receiver / field | Late RMS, lumped (%) | Late RMS, consistent (%) | Consistent late peak / full reference peak (%) |
|---|---:|---:|---:|
| x3 solid velocity | 8.8760 | 0.02563 | 0.12719 |
| x3 Darcy flux | 10.0507 | 0.000492 | 0.002536 |
| d4 solid velocity | 19.5695 | 0.04023 | 0.14163 |
| d4 Darcy flux | 21.0176 | 0.010880 | 0.047746 |

The sampled windows are 0.759375–1.496875 s at x3 and 1.134375–1.496875 s
at d4. Across all eight receivers the largest remaining late peak is **0.142%**
of the full reference peak for solid velocity and **0.048%** for Darcy flux,
both at d4. Late oscillations are greatly reduced but are not zero. Late FAS
uses a Hann window and weighted-mean removal. The former dominant peaks near
34 Hz at x3 and 27 Hz at d4 have much smaller amplitudes; the remaining solid
tails peak near 20 and 19 Hz, respectively. The finite-window resolutions are
about 1.35 and 2.74 Hz; zero padding does not improve those resolutions.

This is substantial evidence that consistent inertia/storage improves this
inviscid 4 m calculation, but it does not establish convergence. Worst
full-record errors remain about **10%**, and worst slow-P solid error is
**11.25%**. A tighter-tolerance comparison with 94436 is not part of this
result; `1e-6` algebraic error has not been isolated from spatial, temporal,
source, stabilization and absorber errors. The finite-top return bound exceeds
1.55 s, but lateral/bottom CPML can still contribute; this comparison alone
does not identify the cause of the small remaining tails. The author archive
has no pressure history. Pressure was checked for complete finite output
(maximum absolute station value **104.064 Pa**), not analytical accuracy.

Reproduction inputs, the two local-only helpers, all receiver/component CSV
tables, JSON provenance and seven PNG/PDF figure pairs are preserved in the private
validation archive; its README specifies units,
dependencies, commands, windows and overwrite behavior. These artifacts and
raw author data are not included in a clean checkout. No solver, job or
reservation settings were changed for this analysis.

## Completed 32-node partition and fusion pilot: 94474/94488 (2026-11-15)

Both jobs completed with exit `0:0` under `IO_testing`: **94474** in **13:15**
and **94488** in **15:44**. The production `fusion-v01` binary has SHA-256
`625a2628bdd6b2ebf0232dcf2ebc4efc9dcf9f1146aed129f8c7ea7becca55ef`.
All eleven individual cases completed successfully, including the two profiled
cases in 94474. Runtime records confirm 32 nodes/128 GPU ranks, CUDA-aware MPI,
and the intended fusion/profiling switches. No private operator/energy audit
is enabled in these production timing runs.

The model remains the **4096 m cube, uniform 4 m mesh**, 1,073,741,824 elements,
zero dangling nodes, 256 m CPML, gamma 0.005, standard CG, element-inverse
preconditioning, previous initial guesses, assembled Q1 rows and **1e-10**
tolerance. Each case advances **50 steps at 50 microseconds**, only **2.5 ms**
of simulation time. Active input files differ only in `cpml_partition_weight`;
the source files are byte-identical. Each matched fusion pair has identical
inputs and changes only `HERCULES_BIOT_MASS_FUSION` in its runtime controls.

The following times are **maximum-rank Compute Physics**, excluding setup and
output. Both members of each pair are unprofiled and run in the same 94488
allocation. All speedups below are measured ratios from these matched runs,
not products of speedups from different meshes or tolerances.

| CPML partition weight | Separate kernels (s) | Fused kernels (s) | Fusion speedup at this weight | Fused speedup versus weight 4.38 separate |
|---|---:|---:|---:|---:|
| 4.38 | 27.690 | 19.558 | 1.416x | 1.416x |
| 2 | 21.058 | 15.561 | 1.353x | 1.779x |
| 1 | 17.658 | 13.342 | 1.323x | 2.075x |

Weight **1 plus fusion** is fastest among these tested combinations, using
**51.82% less physics time** than weight 4.38 with separate kernels. Weight 1
alone gives **1.568x** over the separate-kernel weight-4.38 baseline. The
independent unprofiled baselines in 94474 are **27.690, 21.189 and 17.791 s**
at weights 4.38, 2 and 1; their repeats in 94488 differ by at most **0.75%**.
There is only one fused sample at each weight, so this is not an uncertainty
estimate for repeated production runs.

The mesh statistics support the load-balance explanation. Maximum elements per
rank fall from **16,390,717** at weight 4.38 to **10,756,096** at weight 2 and
**8,388,608** at weight 1. Maximum harbored nodes fall from **16,597,534** to
**10,949,632** and **8,520,321**. Weight 1 distributes equal element counts,
reducing the largest node vectors processed in synchronous mass iterations.
At weight 4.38, partitions with more consistent physical cells carry much
larger iterative vectors than CPML-heavy partitions.

The separate-kernel profiles in 94474 show maximum-rank cumulative inertia
vector-kernel GPU time falling from **8.151 to 4.342 s** between weights 4.38
and 1, and storage vector-kernel GPU time from **5.918 to 3.157 s**. Maximum
inertia operator GPU time falls from **4.092 to 2.101 s**. Rank-average CPU
Allreduce intervals fall from **2.249 to 0.638 s** for inertia and **2.538 to
0.787 s** for storage. These intervals include rank waiting; CPU exchange and
GPU intervals overlap and must not be added as independent costs. They support
reduced imbalance, not a claim that all measured waiting is network overhead.

Each case completed **50 inertia solves and 100 storage solves**, with maximum
15 iterations and worst verified relative residual below `1e-10`. Mean inertia
iterations are **14.04–14.08**, and storage **13.44–13.51**. The gains therefore
do not come from relaxing the solver tolerance or a large iteration reduction.
All 16 station files per case contain 13 finite samples through 2.4 ms. This
window precedes physical receiver arrivals: in 94488 the maximum solid velocity
is only **8.47e-51 m/s** and maximum pressure **1.99e-44 Pa**. Largest same-weight
fusion differences are **1.30e-50 m/s** in solid velocity, **8.70e-51 m/s** in
Darcy flux and **1.45e-44 Pa** in pressure. These near-zero outputs are not a
resolved waveform-accuracy test; use the earlier full-pulse and operator audits
for that separate question.

These results identify weight 1 as the preferred **tested pilot setting for
this consistent-Q1 model**. They do not establish a universal default for
elastic, acoustic, mixed, viscous, nonlinear or topographic models, nor a
full-duration speedup for the 1.5 s Figure 5 calculation. In particular, the
2.075x gain must not be multiplied by unrelated tolerance/pilot speedups.
Short-case process wall times include mesh/setup costs (193.565 to 90.255 s
for the two extremes) and are not used as the propagation-speed metric.
The binary predates compact mass vectors and combined mechanics packets;
those are the separate 94547 comparison. No defaults, running inputs or job
reservations were changed.

The local-only `analyze_partition32.py` and its documented JSON/CSV outputs are
preserved in the private validation archive. They audit
completion, provenance, matching inputs, profiles, mesh balance, residuals and
station data. They are not included in a clean checkout.

### Force-stage cost ratio is distinct from the best tested partition weight

The elastic-derived default **4.38** was based on a force-stage host event-wait
estimate, not complete full-Biot timestep costs. Repeating that style of
calculation with 94488's equal-element ranks gives an effective porous CPML
**force-stage** cost of approximately **13–14.4 times** the ordinary porous
element force cost. It does not establish that ratio for the complete
consistent-mass timestep.

Use the `Force stiffness core wait only` column of the CPML rank/load table
for `partition32_w1_baseline`. All ranks contain **8,388,608 elements** and
advance **50 steps**. CPML-free ranks **34 and 35** both report **0.13 s**.
For a rank with CPML fraction `f`, the two-type linear work estimate is

```math
\frac{T(f)}{T(0)}=(1-f)+r_{\mathrm{force}}f,
\qquad
r_{\mathrm{force}}=1+\frac{T(f)/T(0)-1}{f}.
```

| CPML elements/rank | Exact CPML fraction | Printed force-core wait (s) | Inferred CPML / ordinary force cost |
|---|---:|---:|---:|
| 0 | 0 | 0.13 | reference |
| 2,097,152 | 0.25 | 0.52 | 13.0 |
| 3,670,016 | 0.4375 | 0.84–0.85 | 13.5–13.7 |
| 5,242,880 | 0.625 | 1.17 | 13.8 |
| 6,029,312 | 0.71875 | 1.38 | 14.4 |

For example, the last row gives
`1 + (1.38 / 0.13 - 1) / 0.71875 = 14.3779`. The printed rank values are
rounded to **0.01 s**, so those digits are arithmetic, not measurement
precision. CPML face/edge/corner composition also changes the number of active
stretch directions. These are selected-rank host waits for queued force work,
not isolated CUDA-event kernel durations; fixed overhead is not separately
subtracted. Fusion changes the mass kernels, and the matching fused case
reports essentially the same force-stage times for the printed CPML ranks.

For **the full consistent-mass update**, the existing logs do not provide a
defensible CPML/ordinary per-element cost ratio. In the weight-1 fused case,
rank-zero `Biot inertia PCG` and `Biot storage PCG` take **7.484783 + 4.325908 =
11.810691 s**, **88.68%** of its **13.318284 s** `Compute Physics` timer.
The storage total includes both pressure-update and stabilization solves;
the enclosing stabilization timer must not be added again. These PCG timers
include communication and waiting and are not costs attributable solely to
the interior cells.

Review of the archived `fusion-v01` sources confirms that consistent element
integrals exclude CPML, while the iterative vectors still span **all harbored
nodes**, including CPML nodes. Thus CPML pays common nodal-vector and diagonal
work but does not pay the same consistent-row work as the physical interior.
Pressure continuity/storage and other CPML history updates also occur outside
the force stage. The logged PCG aggregates do not separate local GPU time by
these categories or provide per-rank correlations with cell composition.
Using the **13–14.4 force ratio** as a complete-timestep partition weight would
therefore omit the dominant cost and is not justified.

Consequently, **weight 1 is the fastest tested configuration**, while
**CPML/ordinary whole-step cost = 1 is not a measured result**. Calibrating the
latter requires per-rank local GPU durations for consistent rows, shared
nodal work, and porous CPML updates, with communication/waiting separated.
No such new instrumentation or run is included in this analysis. The local
`analyze_cpml_cost.py` records the force estimates, contributing rank rows and
full-step limitation in `cpml-cost-analysis.json`; it changes no solver inputs,
defaults or jobs.

## Completed 32-node compact-vector/combined-packet pilot: 94547 (2026-11-15)

Job **94547** completed all four cases with exit `0:0` under `IO_testing`,
32 nodes/128 GPU ranks, in **13:15** elapsed. It used the archived production
`compact-v01` binary, SHA-256
`1bd96de1acf1113e6a3442e799476ee8e41380f5fc1bab7936e7b3fd6577ca31`.
The fetched runtime records confirm CUDA-aware MPI, **fusion enabled**,
and mass/MPI profiling and private audits disabled in every case. Compact
vectors and combined solid/fluid packets are independently toggled.

All four cases have identical active physical input settings, byte-identical
source files and byte-identical per-rank mesh statistics. They use the
4096 m cube, uniform 4 m elements, 256 m CPML, gamma 0.005, standard CG,
element-inverse preconditioning, previous initial guesses, assembled rows and
`1e-10` tolerance. **Partition weight stays at 4.38.** The only differing runtime
controls are `HERCULES_BIOT_MASS_COMPACT` and `HERCULES_BIOT_MASS_HALO_BATCH`.
The record is **50 steps at 50 microseconds**, or **2.5 ms**.

| Compact / combined | Physics maximum rank (s) | Change in physics time | Inertia PCG rank zero (s) | Storage PCG rank zero (s) |
|---|---:|---:|---:|---:|
| Off / off | 19.340 | reference | 10.675 | 7.227 |
| On / off | 20.497 | **5.98% more** | 11.383 | 7.710 |
| Off / on | 18.938 | **2.08% less** | 10.330 | 7.168 |
| On / on | 20.080 | **3.83% more** | 10.981 | 7.699 |

This full-size pilot **does not reproduce the compact-vector speedup** measured
on four nodes. Combining packets alone is fastest here (**1.021x**), but there
is only one measurement per configuration, so the small benefit should not be
treated as a precise production speedup. Both changes together give **0.963x**
relative to the already fused baseline. Process wall time is dominated by
setup: maximum-rank wall times are **183.980, 186.967, 186.056 and 186.384 s**
in the table's order. Those are not propagation-only timings.

Compaction retains **786,595,988 of 1,094,454,628 harbored node copies**, reducing
work/history vector length by **28.13% globally**. However, the **maximum rank
remains at 16,597,534 nodes**, identical to the maximum full-vector size.
The minimum retained rank falls to **815,425**, from **4,438,463** full nodes.
Thus the savings do not reduce the largest iterative vector in this weighted
partition. This supports the interpretation that the synchronous solve remains
limited by interior-heavy ranks, while compact indexing/direct-block processing
adds work. The job has no per-kernel profile to isolate those overheads, so that
causal explanation remains an inference. The 28.13% figure is a work/history
vector reduction, not a measurement of total GPU-memory reduction.

All four cases complete **50 inertia and 100 storage solves**. Mean inertia
iterations are **14.080** in every case; storage means are **13.43–13.49**.
Maximum iterations are **15**, and every reported worst true relative residual
is below `1e-10`. Each case has all 16 station files, 13 finite samples through
2.4 ms. The receivers have no resolved wave arrivals: the largest baseline
solid velocity is **1.06e-50 m/s**, and the largest baseline pressure is
**1.08e-44 Pa**. Across the three switch comparisons, maximum absolute
differences are **2.71e-50 m/s** for solid velocity, **1.81e-50 m/s** for Darcy
flux and **4.36e-44 Pa** for pressure. These tiny prearrival fields do not
establish full-wave agreement; the earlier complete pulse and independent
operator tests remain the separate numerical evidence.

The weight-1 fused result from 94488 (**13.342 s**) remains faster than any of
these weight-4.38 results. That cross-job observation does not measure compact
or combined-packet effects at weight 1; **no weight-1 combination is tested in
94547**, and the gains from the two jobs must not be multiplied. A full-duration
run with the new switches is also not established. No solver defaults, case
inputs or running-job reservations were changed after inspecting these results.

The local-only `analyze_pilot32.py`, `pilot32-analysis.json` and two CSV tables
in the private validation archive retain completion,
source/input/runtime checks, timers, map counts and station comparisons. The
helper's README documents dependencies, units, commands and overwrite behavior;
these investigation artifacts are unavailable in a clean checkout.

## Per-rank GPU cost calibration preparation record (2026-11-15)

The maintainer requested direct GPU measurements to distinguish consistent
matrix work, common nodal operations and porous CPML updates before inferring
a complete-timestep partition weight. An **isolated diagnostic build** now
records these categories per rank. Production solver sources, accepted public
inputs and defaults are unchanged; this is instrumentation for the present
investigation, not a newly validated physics path.

The local-only `prepare.py` preserves
a 103-file source snapshot, SHA-256 manifest and source archive. `cost-v03`
was built on Horizon with `ml reset`, the existing site configuration and the
required GSL/PROJ paths; its library check found no missing dependencies.
The executable SHA-256 is
`62bcbbac84549ccd3e8f66ec02b2d1e06f26de92fd9b2596f267b7fc68d864a0`.

The diagnostic exports per-rank CSV records containing:

- Physical/CPML cell counts, full/iterative node counts and consistent rows.
- Separate inertia/storage GPU durations for consistent operators, vector
  operations, dot products, direct blocks and scalar transfers.
- Solid-frame stiffness and CPML kernels; porous traction and continuity;
  stabilization; nodal updates/resets; packet/constraint and station kernels.
- CPU exchange, readiness and global-reduction intervals, **reported separately
  because they overlap GPU work and include waiting for other ranks**.

The additional GPU events use a fixed pool allocated during initialization.
Their completion is queried after existing mass-reduction waits; no new
per-step synchronization is introduced. Final teardown drains remaining
events. CUDA-aware MPI exchange paths are preserved. GPU event intervals are
service-cost measurements, not a measurement of the synchronized timestep's
critical path.

Three private diagnostic modes compare disabled instrumentation, original
mixed launches with profiling, and profiled launches split into physical/CPML
index sets. The split mode preserves the element equations but changes atomic
accumulation order and launch layout. Complete-pulse comparisons and original
versus split timing controls are therefore required before interpreting its
attribution. Setup reads the live device CPML map once and checks counts against
the consistent-mass setup and GPU CPML state.

An initial diagnostic revision incorrectly used a host map already released
by initialization. Jobs 94561/94562/94563 were canceled to replace that revision;
its split-category results are excluded. This was a profiler classification
error, not a change to the solver's CPML equations. Existing production jobs
and their reservations were not changed.

Following the maintainer's reservation update, these replacement jobs were
submitted under **`cuda-615-SC`**:

| Job | Nodes / GPU ranks | Purpose |
|---|---:|---|
| 94565 | 1 / 4 | Complete 0.6 s pulse, diagnostic modes 0/1/2 |
| 94566 | 4 / 16 | Short 2048 m cube pilots, profiling/split controls and weights 1/4.38 |
| 94567 | 32 / 128 | Production-size pilots, tolerance `1e-10`, weights 1/2/4.38 and unprofiled controls |
| 94568 | 32 / 128 | Matching production-size pilots at tolerance `1e-6` |

The production-size pilots retain the 4096 m cube, uniform 4 m mesh, 256 m
CPML, gamma 0.005 and 50 steps at 50 microseconds. They use fused standard CG,
element-inverse preconditioning, previous guesses and assembled rows, with
compact vectors and combined packets disabled to match the 94488 comparison.
The complete pulse uses a 256 m cube and 32 m CPML. These are homogeneous,
flat, inviscid full-Biot cases; neither four-node timing nor this profiling
exercise validates other material/interface/topography formulations.

The local analysis audits completion, residuals, rank counts and runtime
controls, separates GPU components from CPU waiting, and retains raw rank
tables plus component/composition fits. Pulse comparisons include time
histories and Fourier-amplitude spectra. At preparation the corrected runs
were pending; the completed measurements and selection are recorded below. A fitted
additive GPU-cost ratio must also be checked against unprofiled full-step
timings; it is not automatically the optimum scalar partition weight.

## Completed Figure 5 tolerance comparison: 94436 versus 94445 (2026-11-15)

Job **94436 completed successfully**, Slurm state `COMPLETED`, exit `0:0`,
on 32 nodes/128 GPU ranks under `IO_testing`. Its elapsed time was **4:45:02**,
ending 2026-11-15 at 01:33:30 Horizon local time. All eight solid/fluid receiver
pairs contain the expected 7,500 finite samples through 1.4998 s. The output
identities relating Darcy flux, relative fluid velocity and total fluid
velocity pass.

This provides the previously missing **same-binary tolerance comparison**.
Both jobs use `performance-v03`, SHA-256
`399c2575c35bca19c8036bfffa3c7b1c23a30f61e5aceb73219bbdf85aeb8064`.
The only changed active input is `poroelastic_mass_tolerance`: `1e-10` in
94436 and `1e-6` in 94445. All source files, mesh/schedule statistics and
relevant runtime controls match exactly. Both retain the 4096 m cube, uniform
4 m mesh, 256 m CPML, gamma 0.005, 50 microsecond timestep, 16 Hz source and
1.5 s simulation, with consistent physical-cell inertia/storage and lumped
CPML. The solver uses standard CG, element-inverse preconditioning, previous
guesses and assembled rows. The executable predates vector fusion, compact
vectors and combined mass halos; this comparison does not validate those
later changes over a full-duration Figure 5 run.

### Completion, iteration counts and timing

| Quantity | 94436: `1e-10` | 94445: `1e-6` |
|---|---:|---:|
| Slurm elapsed | 4:45:02 | 2:31:26 |
| Compute Physics, maximum rank (s) | 16,862.419 | 8,843.495 |
| Compute Physics, rank zero (s) | 16,859.176 | 8,839.914 |
| Inertia solves / mean iterations / maximum | 30,000 / 14.836 / 15 | 30,000 / 6.112 / 9 |
| Storage solves / mean iterations / maximum | 60,000 / 14.008 / 15 | 60,000 / 5.945 / 9 |
| Worst inertia true relative residual | 9.99766e-11 | 9.99702e-7 |
| Worst storage true relative residual | 9.99975e-11 | 9.99989e-7 |

The looser tolerance is **1.907x faster in maximum-rank physics time** and
**1.882x faster in Slurm elapsed time**, reducing physics time by 47.55%.
These are measured ratios between the two completed jobs, not compounded
speedups from different pilots or a repeated-run uncertainty estimate.

### Accuracy against the authors' reference

Both cases were reanalyzed against the same authors' inviscid Karpfinger
semi-analytical archive, with the fixed author amplitude scale, original
3.125 ms reference grid and physical velocity/flux times `t - dt/2` used in
the earlier 94445 analysis. No fitted amplitude/time corrections or waveform
filters were introduced. Vector waveform L2 includes all three components;
FAS and late-peak metrics below use radial projections.

| Metric, maximum across eight receivers | 94436: `1e-10` | 94445: `1e-6` |
|---|---:|---:|
| Full solid velocity L2 (%) | 9.88748 | 9.88776 |
| Full Darcy flux L2 (%) | 9.71009 | 9.71028 |
| Slow-P solid velocity L2 (%) | 11.24828 | 11.24852 |
| Slow-P Darcy flux L2 (%) | 9.71057 | 9.71076 |
| Solid radial FAS L2, 1–48 Hz (%) | 7.90944 | 7.90881 |
| Darcy radial FAS L2, 1–48 Hz (%) | 7.92072 | 7.92104 |
| Late solid peak / full reference peak (%) | 0.141216 | 0.141632 |
| Late Darcy peak / full reference peak (%) | 0.047749 | 0.047746 |

The remaining late peaks occur at d4 for both fields. Its late solid RMS is
0.040274% versus 0.040230% of the full reference peak, and its late flux RMS is
0.01088045% versus 0.01088005%. The dominant late spectral bins remain
18.75 Hz for solid velocity and 26.5625 Hz for Darcy flux. These finite-window
spectral bins do not imply that level of frequency resolution; the late
window at d4 resolves only about 2.74 Hz. The late window begins at
`0.19 + r/c_slow`, as before.

### Direct difference between tolerances

The following uses **all 7,500 native samples** at 0.2 ms spacing, not the
coarser analytical-reference grid. Relative L2 is the norm of the numerical
difference divided by the corresponding `1e-10` history norm; maxima are taken
over receivers. These are tolerance differences, not analytical errors.

| Field | Maximum waveform difference L2 (%) | Maximum 1–48 Hz FAS difference L2 (%) |
|---|---:|---:|
| Solid displacement | 0.003250 | 0.002546 |
| Solid velocity | 0.011268 | 0.003291 |
| Relative fluid displacement | 0.003193 | 0.003171 |
| Darcy flux | 0.004727 | 0.004623 |
| Pore pressure | 0.004723 | 0.004633 |

The largest absolute pressure difference is **0.002409 Pa**. There is no
analytical pressure history in the author archive, so this establishes
tolerance agreement for pressure, not its analytical accuracy.

**For this benchmark, `1e-6` is supported by a full-duration comparison with
`1e-10`.** Tightening the solve does not materially improve the roughly 10%
reference discrepancy or remove the late oscillations. Therefore iterative
tolerance is not their dominant cause here. This comparison does not separate
spatial, temporal, source, stabilization and absorbing-boundary errors, and
does not establish convergence or a universal tolerance for other cases.

The private tolerance study retains
`analyze.py`, `compare_tolerance.py`, usage/dependency documentation, accounting
records, JSON audits, receiver/component CSV tables, and six PNG/PDF figure
pairs covering waveforms, slow arrivals, full/late spectra and all receivers.
The artifacts and author data are unavailable in a clean checkout. No solver
settings, running jobs or reservations were changed for this analysis.

## Completed original consistent-mass Figure 5 run: 94379 (2026-11-15)

Job **94379** completed with exit `0:0` on 32 nodes/128 GPU ranks under
`cuda-615-SC`, with Slurm elapsed **7:21:55**. It used `consistent-v04`, SHA-256
`5f3f6b296735c83fe309bff8502cd904063897161f700ac63292bcb748cad476`.
The model, source files, 4 m mesh, 50 microsecond timestep, gamma, CPML and
1.5 s record match the subsequent 94436/94445 comparisons. This original
implementation used the original consistent-mass solver controls, before the
later preconditioner, initial-guess and assembled-row optimization controls.
The tolerance is `1e-10`. The complete station data are finite and satisfy the
solid/fluid output identities.

| Full-duration run | Maximum-rank physics (s) | Speedup relative to 94379 |
|---|---:|---:|
| 94379: original consistent, `1e-10` | 26,283.386 | 1.000x |
| 94436: optimized solves, `1e-10` | 16,862.419 | 1.559x |
| 94445: optimized solves, `1e-6` | 8,843.495 | 2.972x |

These are matched full-duration measurements, not extrapolations from the
four-node pilot. The jobs used different node allocations and executable
revisions, so these observed run-time ratios do not isolate allocation effects
from implementation effects. The 94379 reference comparison gives worst full-record vector
L2 errors **9.887476% solid velocity** and **9.710088% Darcy flux**, effectively
the same as 94436. Direct comparison against 94436 at all native samples gives
maximum relative L2 differences of **1.354e-5% solid velocity**, **4.505e-7%
Darcy flux**, and **9.427e-7% pressure** across the eight receivers. Largest
absolute differences are `1e-12 m/s`, `2.890e-14 m/s`, and `1e-5 Pa`,
respectively; finite station-print precision contributes to these tiny values.

This closes the full-duration check of the earlier solve optimizations for
this inviscid benchmark. It does not establish full-duration validation of the
subsequent fusion, compaction or combined-packet changes. The latter require
the separate final optimized run. The fixed author reference, normalization,
time alignment and metric definitions are those documented for 94436/94445.
Accounting, audits and comparison artifacts are retained locally.

## Double-precision mass tolerance default: `1e-6` (2026-11-15)

The maintainer approved changing the optional consistent-mass solver's double-
precision `poroelastic_mass_tolerance` default from `1e-10` to **`1e-6`** after
the completed 94436/94445 comparison. Explicit input values still override the
default. Single precision retains `1e-5`, its existing admissible lower bound.
The accepted ranges, true-residual acceptance check, iteration limit and
default `lumped` mass formulation are unchanged. This does not introduce a
consistent-mass path for reduced `brick_up`; that existing restriction remains.

The production-source `default-v01` executable was built on Horizon with
`ml reset`, the existing site configuration and required GSL/PROJ paths, with
no missing runtime libraries. SHA-256:
`bcc451df38c4f679c86989cd9e50f7e91f1926b7f5b00a83590ea4ca1774cca3`.
Its source manifest differs from the preceding `compact-v01` manifest only
in `psolve.cc`, for this initializer change. The private subunit CPML partition
experiment is a separate executable and is not part of this default change.

### Default selection and applicable-family regression checks

Job **94572** completed all **38 runs** with exit `0:0` under `cuda-615-SC`,
one node/four GPU ranks, in **10:19**. Nineteen paired inputs omit the tolerance
or specify `1e-6`, with otherwise identical inputs and runtime controls.
Fourteen consistent-mass pairs cover inviscid/viscous, heterogeneous,
anisotropic and sealed materials; FEM/VT geometry; open-pore and no-flow
acoustic interfaces, point projection and acoustic-potential coupling;
poroelastic–elastic coupling; poroelastic–viscoelastic and nonlinear
elastoplastic coupling; and drained FEM. Five inactive-control pairs cover
lumped full Biot, acoustic, acoustic–elastic, elastic and viscoelastic cases.
Every active case reports `1e-6`, converges within the selected true-residual
tolerance and produces the expected station outputs. Undefined nonporous
pressure retains its existing output convention.

Eighteen pairs meet the provisional maximum-difference threshold of 0.1% of
each field's peak, with absolute floors for near-zero fields. The short viscous
pair is **retained as a threshold flag**, rather than silently relaxed: its
20 ms record precedes the earliest physical arrival at the nearest receiver
(128 m from the source), and has extremely small numerical tails. For example,
its velocity difference is `7.9e-13 m/s`, while some flux differences are of
order `1e-12 m/s`. Relative percentages against those near-zero records are not
full-wave accuracy measures. Both inputs report the same tolerance, but GPU
accumulation and tolerance-based stopping do not guarantee bitwise repeats.
These short cases check initialization, convergence and applicable paths;
they do not establish reference accuracy for every coupling or constitutive law.

### Complete viscous pulse follow-up

To assess resolved signals, job **94576** ran the existing **0.6 s viscous
pulse**, with viscosity `0.001 Pa s`, under identical fused/compact/combined-
packet settings: omitted `1e-6`, explicit `1e-6`, and explicit `1e-10`.
All three completed with exit `0:0`, valid residuals and finite output, on one
node/four GPU ranks under `cuda-615-SC`, in **10:16** total.

| Field | Omitted versus explicit `1e-6`, maximum L2 (%) | Either `1e-6` versus `1e-10`, maximum L2 (%) |
|---|---:|---:|
| Solid displacement | 0.0000151 | 0.010873 |
| Solid velocity | 0.0000449 | 0.017768 |
| Relative fluid displacement | 0.004630 | 0.083228 |
| Darcy flux | 0.0002132 | 0.077745 |
| Pore pressure | 0.0007893 | 0.018657 |

These are numerical-run differences normalized by the comparison history,
maximized across receivers, not analytical errors. The corresponding worst
full-record FAS differences against the tight run are **0.01261% velocity**,
**0.05632% Darcy flux**, and **0.02050% pressure**. Maximum weighted phase RMS
for these three fields is below **0.0071 degrees**. The full pulse distinguishes
the short prearrival sensitivity from agreement of resolved waveforms. Together
with Figure 5, it supports the requested default for these tested cases;
tighter explicit tolerances remain available for more demanding error targets.

The parser/documentation audit now checks the precision-specific initializer
and rendered defaults. Frozen-lockfile installation, lint, TypeScript, Sites
production/rendered audits and GitHub Pages export/rendered audits pass.
Those documentation checks are separate from the successful Horizon build
and runtime cases. Local `default-analysis.json` preserves all short-case flags;
`viscous-analysis.json` preserves the complete-pulse waveform/spectral metrics.

## Completed per-rank GPU cost measurements (2026-11-15)

The corrected diagnostic jobs **94565–94568** completed with exit `0:0`
under `cuda-615-SC`. Their elapsed times were 18:09, 5:50, 14:26 and 12:08,
respectively. This supersedes the pending status in the preparation record
above. The binary and measurement modes are those recorded there. No private
profiling hooks are added to the production solver.

The 0.6 s pulse checks show maximum field-history relative L2 differences of
`1.3472e-10` for the original profiled launches and `1.2163e-10` for the split
profiled launches, relative to the unprofiled run (fractions, not percentages).
The largest corresponding FAS differences are below `9.6e-11`. Profiling does
change execution time: the small pulse's maximum-rank physics time rises from
226.177 s to 409.633/403.972 s. Profiled elapsed times therefore must not be
used to claim a production speedup.

The production-size experiments distinguish consistent-matrix applications,
common nodal work, ordinary-element skeleton/traction/continuity work, and
CPML histories/traction/continuity work. The following fitted GPU service
costs use the weight-1, split-launch measurement, normalized per cell per
step. Common work includes the fitted nodal cost times the measured number
of node copies per cell, plus stabilization. All figures are **nanoseconds per
cell per step**, not the elapsed time of an individual cell execution.

| Component | Tolerance `1e-10` | Tolerance `1e-6` |
|---|---:|---:|
| Consistent operator, ordinary cells | 9.789 | 5.657 |
| Ordinary skeleton, traction and continuity | 0.601 | 0.601 |
| Common nodal work and stabilization | 8.218 | 5.102 |
| CPML histories, traction and continuity | 5.605 | 5.602 |
| Ordinary cell including common work | 18.609 | 11.360 |
| CPML cell including common work | 13.823 | 10.704 |
| CPML / ordinary ratio | **0.743** | **0.942** |

Repeating this component fit at partition weights 2 and 4.38 gives ratios
0.744/0.758 at `1e-10` and 0.944/0.963 at `1e-6`. Component-fit relative RMS
residuals in these cases are below 3%. An independent composition fit to the
original unsplit `1e-10` launches gives **0.776**, with relative RMS 0.845%.
The matching original-launch `1e-6` control, job **94584**, completed with
exit `0:0` in **1:50** on the same 32-node reservation and gives **0.980**, with
relative RMS **0.752%**. Thus the original-kernel check independently supports
a cost ratio close to one at the new tolerance.
The four-node controls give 0.737–0.751 from split components and 0.792 from
original launches. The few-percent method dependence is retained explicitly;
these are conditional cost estimates, not a universal coefficient.

The ordinary cells become cheaper at `1e-6` because the iterative work falls,
while the direct CPML update cost barely changes. This explains why the cost
ratio moves toward one. The earlier **13–14.4 force-stage ratio** excludes
consistent inertia/storage and common vector work and must not be substituted
for this broader cost estimate.

GPU packet/constraint work and station output are recorded separately; CPU
exchange/readiness/reduction intervals overlap GPU work and are not added to
these costs. Consequently, this additive model does not predict the exact
synchronized timestep optimum. Whole-step, unprofiled measurements provide the
necessary check:

| 32-node pilot | Weight 4.38 physics (s) | Weight 1 physics (s) | Observed speedup |
|---|---:|---:|---:|
| `1e-10`, fused, full vectors, separate halos | 19.023254 | 13.289764 | 1.431x |
| `1e-6`, fused, full vectors, separate halos | 11.964935 | 8.557816 | 1.398x |

The complete `Solver` timer includes communication and output outside the
`Compute Physics` regions. It falls from 19.585864 to 14.071194 s at `1e-10`
(**1.392x**) and from 12.528078 to 8.981334 s at `1e-6` (**1.395x**).
Partition selection uses this complete-loop measurement; physics timings are
also retained for comparison with the earlier records. Total wall time adds
mesh construction and setup and is a separate quantity.

Each pair runs within the same allocation, uses original unprofiled kernels,
and differs only in partition weight. These are 50-step, prearrival timing
pilots; their speedups must not be multiplied by ratios from other runs or
reported as measured full-record speedups. The next comparison includes
compact vectors and combined halos and brackets weight 1 from both sides.
Only the complete final Figure 5 run can establish their full-record timing
and waveform agreement. Raw per-rank data, fits and plots are retained in the
private validation archive, unavailable in a clean
checkout. The scope remains flat, homogeneous, inviscid full Biot with
consistent physical Q1 mass and lumped CPML, on this mesh and allocation size.

## Completed weight and compact-vector selection: 94575/94582 (2026-11-15)

Jobs **94575** and **94582** completed all eleven unprofiled cases with exit
`0:0`, on 32 nodes/128 GPU ranks under `cuda-615-SC`, in **17:45** and
**7:02**, respectively. They retain the same 50-step Figure 5 timing model,
standard CG, element-inverse preconditioner, previous guesses, assembled rows
and fused kernels. Combined solid/fluid halo packets are enabled throughout.
The tolerance is omitted, selecting the new double-precision `1e-6` default.
All cases have identical non-weight physical inputs/source files and finite
station output; every true residual is below `1e-6`. Mean inertia/storage
iterations are 7.280/6.690 in all cases.

A separate private `partition-v01` executable permits weights below one for
this experiment, retaining ordinary-cell integer cost 100 and using rounded
`100*weight` for CPML cells. Its SHA-256 is
`fc6e91f5cb117beebfd0999da7b6f73cd1c4b19c9476e0bb7cb578a3b8496ca9`.
This changes partition ownership, not the physical equations. The public
parser still requires `cpml_partition_weight >= 1`; its default remains 4.38.
No subunit-weight feature or default change is released from this experiment.

| CPML weight | Full vectors: complete Solver (s) | Compact vectors: complete Solver (s) |
|---:|---:|---:|
| 0.25, private | — | 9.723635 |
| 0.5, private | 9.721745 | 9.573076 |
| 0.8, private | 9.378427 | 9.546501 |
| **1** | **8.918996–8.936533** | 9.194043 |
| 1.25 | 9.439891 | 9.653414 |

The three weight-1 full-vector measurements are 8.936533 and 8.920729 s in
94575 and 8.918996 s in 94582, a spread below 0.2%. Their corresponding
maximum-rank `Compute Physics` times are 8.501697, 8.475086 and 8.474170 s.
The complete `Solver` timer, including communication/output outside the
physics regions, is the selection metric. These observations identify
**weight 1 with full vectors as the fastest tested configuration**; they do
not establish a continuous optimum or a universal setting for all models.

Compact vectors retain about 72% of all node copies. At weight 1, however,
the maximum retained rank still has all 8,520,321 nodes. Lowering the weight
reduces that maximum to 8,081,383, 7,369,217 and 6,769,409 at 0.8/0.5/0.25,
but the complete-loop measurements remain slower. Reduced vector storage or
an additive GPU cost estimate alone therefore does not establish an elapsed-
time improvement. The compact path remains available and its public default
is unchanged; the selected configuration is an explicit benchmark override.

The full-duration confirmation is job **94592**, submitted under
`cuda-615-SC` on 32 nodes/128 GPU ranks, using the production `default-v01`
binary recorded above. It preserves the complete 1.5 s physical model, source,
receivers, 4 m mesh, 50 microsecond step, gamma 0.005 and 256 m CPML from 94445.
The only parameter-file differences from 94445 are omitting the tolerance
(previously explicit `1e-6`) and adding `cpml_partition_weight=1`. Runtime
settings explicitly enable fusion and combined halos, disable compact vectors,
and disable both profilers. The source manifest matches all 102 current
production files. The private partition/profiling executables are not used
for this confirmation.

At this preparation point 94592 is pending; no full-record accuracy or speedup
is inferred from the short pilots. The final comparison must audit all native
samples, fixed-reference waveforms/spectra, late windows, pressure agreement,
true residuals and completed full-run timers. Local inputs, provenance,
`weight-analysis.json` , CSV/plots and readers are retained in the private
validation archive, unavailable in a clean checkout.

## Completed full-duration optimized confirmation: 94592 (2026-11-15)

Job **94592 completed successfully** with exit `0:0`, on 32 nodes/128 GPU
ranks under `cuda-615-SC`, in **1:15:57** Slurm elapsed time, ending at
04:39:51 Horizon local time. This supersedes the pending status in its
preparation record. It uses the production `default-v01` executable, SHA-256
`bcc451df38c4f679c86989cd9e50f7e91f1926b7f5b00a83590ea4ca1774cca3`,
with the physical inputs and selected controls documented above. The omitted
tolerance is logged as `1e-6`. No reservation or input of an existing run was
changed; all new overnight jobs used `cuda-615-SC`.

The complete audit verifies the expected 1,073,741,824 elements, zero hanging
nodes, 303,038,464 porous CPML cells, identical source files, and all eight
solid/fluid receiver pairs with **7,500 finite samples** through 1.4998 s.
The fluid output identities pass. The 30,000 inertia solves average 6.111
iterations (maximum 9), with worst true residual `9.99854293e-7`. The 60,000
storage solves average 5.939 iterations (maximum 9), with worst true residual
`9.99967704e-7`. Both satisfy the selected tolerance.

### Measured complete-record performance

| Run | Slurm elapsed | Maximum-rank physics (s) | Complete Solver loop (s) |
|---|---:|---:|---:|
| 94379: original consistent, `1e-10` | 7:21:55 | 26,283.386 | 26,330.142 |
| 94436: optimized solves, `1e-10` | 4:45:02 | 16,862.419 | 16,914.951 |
| 94445: optimized solves, `1e-6` | 2:31:26 | 8,843.495 | 8,896.895 |
| **94592: fused, combined halos, weight 1, `1e-6`** | **1:15:57** | **4,444.796** | **4,456.752** |

Relative to 94445, the observed physics speedup is **1.990x**, and the complete
solver-loop speedup is **1.996x**. Relative to the original 94379, they are
**5.913x** and **5.908x**. These are measured full-record ratios, not products
of unrelated pilot gains. Slurm elapsed includes launch/setup/teardown; the
Solver timer includes communication and output outside the Compute Physics
regions. The jobs use different allocations and executable revisions, so the
observed overall ratios do not isolate every optimization or allocation effect.

### Numerical agreement and reference accuracy

All native output samples are compared directly with 94445, without fitted
shifts, gains or filtering. Maximum differences across receivers are:

| Field | Waveform L2 (%) | All-frequency FAS L2 (%) | Weighted phase RMS (degrees) |
|---|---:|---:|---:|
| Solid displacement | 0.00004086 | 0.00003275 | 0.00001044 |
| Solid velocity | 0.00009266 | 0.00006584 | 0.00001285 |
| Relative fluid displacement | 0.00015309 | 0.00016625 | 0.00002285 |
| Darcy flux | 0.00005257 | 0.00005209 | 0.00001640 |
| Pore pressure | 0.00005387 | 0.00005286 | 0.00001171 |

Native FAS differences use all Fourier bins; phase uses `0 < f <= 48 Hz`,
weighted by baseline amplitude squared, retaining components at least 1% of
that field's maximum baseline FFT amplitude. These are numerical consistency
metrics, not analytical errors. The source archive has no analytical pressure
or displacement histories for this comparison.

The fixed author Karpfinger reference comparison retains the established
normalization and `t-dt/2` alignment for velocity/flux. Worst full-record
vector L2 errors are **9.887768% solid velocity** and **9.710282% Darcy flux**.
Worst radial 1–48 Hz FAS errors are **7.908809%** and **7.921039%**, respectively,
using the same Tukey 0.1 window and zero-padding convention as the earlier
reference comparisons. These errors remain essentially those of 94445;
the optimization preserves the computed solution rather than improving
spatial/source/stabilization accuracy.

Late windows begin at `0.19 s + r/c_slow`. The largest late velocity/flux
peaks are **0.141648% / 0.047747%** of the corresponding full author-reference
peaks. Late maximum differences from 94445, normalized by its full numerical
peak, are only **0.00004468% velocity**, **0.000004028% Darcy flux**, and
**0.000003015% pressure**. The existing late oscillations persist; these
measurements do not attribute their origin or establish convergence. The
complete tail shows numerical agreement with the earlier stable run, not
removal of the underlying Figure 5 discrepancy.

This closes the full-duration confirmation of the **selected** fused,
combined-halo, weight-1, full-vector configuration with the new double default.
Compact vectors are explicitly **off** in this run; its result does not supply
full-duration validation of compact vectors. Their complete-pulse and coupled
integration evidence remains separate. The measured weight and timings apply
to this flat, homogeneous, inviscid consistent-Q1 case and allocation size;
no universal physics-dependent partition default is introduced. Existing
poroelastic CPML/interface restrictions remain unchanged.

The local-only `Hercules_Full_Biot_Overnight_Validation_2026-11-15.pdf` includes the
default checks,
per-rank cost calibration, repeated weight comparisons, complete runtime/
accuracy charts, waveforms and full/late spectra. Its JSON/CSV data, readers,
source manifests and scheduler records are retained locally; these investigation
artifacts and author archives are unavailable in a clean checkout. The
paper/dataset provenance is Wolf et al., *JCP* 455 (2022), 110886,
DOI `10.1016/j.jcp.2021.110886`, and dataset `10.5281/zenodo.5236133`, which
credits Martin Galis for the homogeneous histories produced with Karpfinger's
code. No unrequested community publication or merge accompanies this work.

## Optional extrapolated guesses and aggregate halo packing (2026-11-15)

The consistent-mass solver now accepts `linear` and `quadratic` in
`poroelastic_mass_initial_guess`, in addition to `zero` (unchanged default)
and `previous`. Each inertia stage, pressure increment and stabilization solve
has its own fixed-step history. Startup uses available lower orders; a
nonfinite or worse-than-zero extrapolation is discarded. Every production solve
retains its original projected-RHS normalization and recomputed true-residual
acceptance. This changes the algebraic initial guess, not the Q1 spaces,
physical equations, damping integration or pressure stabilization.

`HERCULES_BIOT_MASS_HALO_PACK=1` optionally gathers all peers and unpacks their
contributions in one GPU launch per exchange phase. Its default is `0`.
Setup groups receives by destination node in fixed peer/mapping order. Existing
CUDA-aware and host transports, receive-reuse fences, full/compact vectors,
separate/combined mechanics packets, scalar storage and field constraints remain.
Physical fields, mesh/source/station indexing and the mass operators are unchanged.

The immutable `predictor-v01` and `combined-v01` executables were built on
Horizon with `ml reset`, checked dependencies and archived source manifests.
Their SHA-256 values are respectively
`03a09fc2f103a6c9263436a1482dc5c805dc7c87da409607c479169989defa40`
and `29306c344e162b67a15d8283a54fc99f08e94edd82f3b0ef69c4e0736c287c9c`.
Job **94733** completed six 0.6 s inviscid/viscous pulse cases on one node,
under `cuda-615-SC`, comparing previous, linear and quadratic guesses at `1e-6`.
The inviscid previous/linear/quadratic complete Solver times were
**113.806 / 76.346 / 53.652 s**. Average inertia/storage iterations fell from
**8.430/8.709** to **2.551/2.047** with quadratic prediction; all checked true
relative residuals remained below `1e-6`. Worst inviscid waveform differences
from previous guesses were **0.00328% v**, **0.000305% q** and **0.000313% p**
(relative vector L2 across each receiver, then maximum over receivers).

Job **94741** separately compared aggregate packing with quadratic guesses.
The inviscid complete Solver time fell from **54.168 to 48.849 s (1.109x)**.
These are small-domain pulse measurements, not full Figure 5 or 32-node
speedups. No independent pilot speedups are multiplied. The dedicated
arbitrary-field packet audit in **94760/94761** found exactly zero difference
from per-peer kernels in all 88 tested scalar/mechanics exchange phases,
covering full/compact maps, FEM/VT geometry, a zero-active-rank fixture and
CUDA-aware/host MPI. These are packet-operator checks, not physical interface
accuracy measurements.

**Weak viscous fluid signals need a tolerance-specific comparison.** The first
`1e-6` packing pair differed by 0.634% q L2. In follow-up **94758**, repeating
an unchanged unpacked `1e-6` run differed by 0.871% q L2, whereas the packed
versus unpacked `1e-10` pair differed by only **1.36e-9% q L2**. This supports
sensitivity to iterative tolerance/reduction roundoff rather than packet-map
corruption. In **94766**, both `1e-8` variants differed from their matching
`1e-10` reference by at most **0.00231% q L2**, **0.000424% v L2** and
**0.000274% p L2**. The `1e-6` default remains unchanged; it is not a universal
relative-error bound for weak fluid fields.

The broad short regressions completed standalone and coupled acoustic,
elastic, viscoelastic, nonlinear elastoplastic and poroelastic cases, including
heterogeneous/anisotropic/sealed material, FEM/VT, point/open/no-flow acoustic
interfaces and both MPI backends. They verify successful execution and
convergence, not blanket waveform accuracy: short prearrival records can have
large relative errors in very small signals. Completed longer coupled comparisons
and 32-node resolved-wave partition tests are recorded below.
Single-precision and restart propagation are not established by these tests.

Private experiments with other stopping norms, periodic verification and
uniform-box scratch ordering are isolated from production. In particular,
production does not skip per-solve true-residual verification. Local helpers,
raw outputs, failed/replaced private prototypes and source/binary manifests
remain in the private validation archive, unavailable in
a clean checkout. No default partition weight, physical source/material or
existing job reservation was changed.

Resolved 0.3 s comparisons in **94759/94775/94777** completed FEM open-pore,
VT open/no-flow, drained FEM, heterogeneous/anisotropic/sealed material and
coupling to elastic, viscoelastic and nonlinear elastoplastic solids. For the
FEM/VT acoustic interfaces and drained FEM, worst receiver differences in all
five native fields were below **0.0012% L2**, with complete-loop speedups of
**1.88–2.10x**. These are before/after numerical comparisons, not independent
physical interface or constitutive-law calibrations.

Other material/solid-coupling fixtures showed larger `1e-6` differences, so
they are not covered by that accuracy statement. The anisotropic on-axis
fluid signal is nearly zero (peak q about `1e-13 m/s`, maximum difference
about `3.1e-15 m/s`); retain absolute values alongside its large percentages.
In **94790**, tightening the same previous-versus-quadratic comparisons to
`1e-8` reduced worst v L2 differences to **0.00391% heterogeneous**,
**0.00777% anisotropic** and **0.00822% sealed**; p differences were below
**0.00544%**. Sealed relative-fluid fields remained exactly zero. These
fixtures preserve their existing Krylov and operator settings.
The freshly rebuilt production `combined-v03` binary used for these checks is
`d6b3f8da7c6682a173f93e4560483c26dfa8675cbce232f1e27b32c2881e671f`.
Its only source difference from `combined-v01` is size-safe grid arithmetic
for aggregate packet launches; both CPU and CUDA solver units were rebuilt.

Private stopping-rule experiments **94751/94752** retain diagnostic work and
therefore provide numerical/energy evidence rather than speedup measurements.
The separate **94787/94788** experiments actually omit selected sweeps and
include paired all-solve audits. Their combined lumped-RHS normalization and
every-tenth-solve verification mode measured **48.598 → 40.420 s (1.202x)**
inviscid and **91.062 → 81.638 s (1.115x)** viscous on the complete small
0.6 s pulses. Both independent audits found zero violations of the original
`1e-6` criterion. The timing variant's printed residuals are explicitly
accepted recurrence/checked residuals, not guaranteed recomputed residuals
for every solve. These private results do not change production acceptance.
The changed normalization alone gave no clear viscous gain (**1.002x**):
average storage iterations increased from **3.963 to 4.617**, offsetting
the saved normalization stencil. These measurements quantify the complete
algorithm rather than assuming that removed sweeps translate directly to speed.

The private uniform-box scratch-ordering/tensor-stencil prototype passed all
eight direct mass/preconditioner product audits in **94776**, with exactly
zero difference from stored rows on the tested arbitrary fields. Small pulse
timings improved by 2.7% inviscid and 5.5% viscous. However, the four-node,
16-GPU **94778** comparison (2048 m cube, uniform 4 m mesh, 0.3 s) measured
**504.542 versus 503.910 s (1.0013x)**: no meaningful complete-loop gain at
the production case's average element count per GPU. Native v/q/p differences
were below **0.000013% L2**. This restrictive prototype remains local-only;
it does not add a production mesh-ordering or implicit-stencil option.

Independent full-size audit **94798** checked all **35,997** accepted solves
through 0.6 s on 32 nodes, with zero violations of the original `1e-6`
criterion (worst **9.292803e-7**). Its u/v/w/q/p differences from the strict
production control were below **0.0012% L2**. Its extra verification work is
excluded from performance claims. Completed full-duration experimental and
production confirmations are recorded separately below.

The remaining `1e-8` solid-coupling follow-ups in **94795** completed with
worst v/q/p receiver L2 differences of **0.00272/0.00292/0.00326%** for
elastic no-flow, **0.00679/0.00812/0.00528%** for nonlinear elastoplastic
coupling, and **0.00323/0.00378/0.00362%** for viscoelastic coupling.
The corresponding complete-loop gains were **3.49x / 2.47x / 2.55x** against
previous guesses with per-peer packing at the same `1e-8` tolerance. These
are resolved 0.3 s numerical comparisons; they do not replace independent
physical validation of the constitutive models.

The five distinct completed 32-node, 0.6 s cases in **94747** use quadratic
guesses and retain every production true-residual check. Complete Solver
times are **1163.445 s** at weight 1 with per-peer packing, **1170.265 s**
at weight 1 with aggregate packing, **1211.408 s** at weight 1.5,
**1227.534 s** at weight 2, and **1467.198 s** at weight 4.38 (the latter
three use aggregate packing). Weight 1 is the fastest tested partition;
aggregate packing has no clear complete-loop benefit at this scale. These
are whole-step measurements, not a universal per-element CPML/interior cost
ratio. The per-peer weight-1 configuration was selected for the separate
full-duration production confirmation **94833**, recorded below. The weight-1
aggregate repeat completed in **1166.857 s**,
only 0.29% different from its first run. Both packed timings remain close
to the per-peer control; they do not establish a packing benefit at this
scale. All six cases satisfy the original true-residual criterion.

The separate full-duration private case **94811** completed 1.5 s on 32
`IO_testing` nodes in **2511.361 s Solver / 2587.422 s total wall time**.
This is **1.775x** faster than the completed **4456.752 s** Solver time of
94592, and **2.510x** the **1000.715 s** lumped Solver time of 94139. These
are measured full-run ratios, not products of independent pilot gains.
Mean inertia/storage iterations were **1.776/1.890**. It combines quadratic
guesses, weight 1, per-peer packing, the private lumped-RHS normalization
and every-tenth-solve verification. Its reported residuals are accepted
residuals; this timing run does not independently verify every solve.

All eight station pairs contain 7,500 finite samples and satisfy the fluid
output identities. Worst native differences from 94592 were
**0.01170% v**, **0.00514% q**, and **0.00513% p** L2. Against the fixed
author reference, worst v/q L2 errors were **9.88858/9.71145%**, compared
with **9.88777/9.71028%** for 94592. Thus the optimization preserves the
existing reference agreement; it does not remove spatial discretization
error or late oscillations. The private full-duration run is separate from
the production confirmation and the shorter all-solve audits.
Against the available full `1e-10` numerical reference **94436**, the same
private run differs by **0.00446% v**, **0.00216% q**, and **0.00216% p** L2.
This comparison includes solver-control/partition differences and measures
agreement with the tightly solved discrete model, not spatial convergence.

Fresh private audit binary `criteria-fast-v05` adds explicit finite/nonnegative
checks and correct zero-denominator handling to the independent diagnostic.
It was rebuilt in both CPU and CUDA translation units (SHA-256
`0058f4895deb434077f10f1feeeec504a9cfe93ef03823a92c47b14aa7254e21`).
Jobs **94823–94827** completed all **13** audit cases: the large 0.6 s model,
inviscid/viscous pulses, and ten resolved 0.3 s material/geometry/coupling
fixtures. Across **827,955** accepted solves, the independent audit found
**zero violations**, with worst original-norm residual **9.987292e-7**.
Reported per-role audit counts equal the inertia/storage solve totals.
Audit timings are excluded from speedup claims. General coupled private
tests use standard CG and full vectors on both sides of each comparison.

This algebraic result does not establish equal relative waveform accuracy
for vanishing fields. For example, the anisotropic on-axis strict flux
peaks near `1e-13 m/s`, while the changed criterion produces a difference
of order `1e-10 m/s`; raw results retain those absolute values and large
relative percentages. Production retains its original norm and recomputed
acceptance residual for every solve.

The paired large-case timing retry **94801** completed both 0.6 s cases on
the same 32-node CUDA allocation. With the same private executable, source,
physical inputs, quadratic guesses and weight 1, original checks took
**1163.650 s**, while lumped-RHS normalization plus every-tenth-solve
verification took **1017.155 s**: **1.144x** faster (12.59% less Solver time).
Mean inertia/storage iterations increased slightly, **1.899/1.967 →
1.967/2.029**, while omitted sweeps supplied the net saving. Worst native
differences were **0.000608% v**, **0.000745% q**, and **0.000660% p** L2;
all five fields were below **0.00166%**. This is the incremental gain after
quadratic prediction, not an independently multiplicative earlier pilot.
The private log labels both variants' residuals as accepted; the original
variant explicitly selects norm 0 and verification interval 1, while the
experimental variant selects norm 1 and interval 10. Independent all-solve
audits are the separate jobs described above.

**Full production confirmation 94833 is complete.** It used the production
`combined-v03` binary, quadratic per-role time prediction, weight 1, full
vectors, standard CG, element-inverse preconditioning, assembled rows, fused
passes, combined solid/fluid packets and per-peer packing. All original RHS
normalization and per-solve true-residual checks were retained at `1e-6`.
The 1.5 s run on 32 nodes took **2892.591 s Solver / 2880.261 s maximum-rank
physics / 2968.251 s total wall time**, with Slurm elapsed **00:49:52**.
The job was submitted to `cuda-615-SC`; the maintainer reassigned it to
`IO_testing` and confirmed that change. The agent made no reservation changes.

Compared with 94592's **4456.752 s** complete Solver timer, this is a measured
**1.541x speedup / 35.10% runtime reduction**. It is **2.891x** the historical
lumped Solver timer of 94139. These full-run ratios span allocations/revisions
and do not isolate every prior optimization. Mean inertia/storage iterations
fell from **6.111/5.939 to 1.747/1.861**. All **90,000** solves passed their
original true-residual check; worst inertia/storage residuals were
**9.999989e-7 / 9.999904e-7**.

All eight receiver pairs contain 7,500 finite samples and satisfy the fluid
output identities. Worst native differences from 94592 were **0.01152% v**,
**0.00517% q**, and **0.00516% p** L2; corresponding FAS differences were
**0.00968% / 0.00515% / 0.00515%**. Against the full `1e-10` numerical
reference 94436, differences were **0.00460% v / 0.00224% q / 0.00224% p**.
Against the fixed author histories, worst v/q errors were
**9.88861% / 9.71151%** L2; late peaks were **0.14116% / 0.04776%** of
reference full peaks. The performance changes preserve the existing
reference agreement and do not establish spatial convergence or removal
of late oscillations.

The measured recommendation for this homogeneous Figure 5 case is quadratic
time prediction with weight 1, full vectors and per-peer packing. Q1 fields,
physical equations and outer time integration are unchanged. Aggregate packing
remains optional; the private stopping rules and restrictive uniform-grid
prototype are not production options. No new default is inferred from these
case-specific measurements. The local PDF report, all numerical comparisons,
source/binary manifests and failed/replaced attempts are retained in the private
validation archive. The final PDF is
`Hercules_Biot_Predictor_Performance_Validation_2026-11-15.pdf`; these local
artifacts are unavailable in a clean checkout.

## Production solver defaults after maintainer approval (2026-11-15)

The completed production options above are now defaults for the optional
consistent-mass solver: `poroelastic_mass_preconditioner=element_inverse`,
`poroelastic_mass_initial_guess=quadratic`, `poroelastic_mass_operator=auto`,
and `HERCULES_BIOT_MASS_COMPACT=0`. Existing explicit settings still override
them. The mass formulation remains `lumped` by default. Double/single tolerances
remain `1e-6`/`1e-5`, with standard CG, fused local kernels, combined solid/fluid
halo packets, per-peer packing, and the original RHS norm and true-residual
check on every solve. The private stopping-rule experiments are not promoted.
Earlier sections describe the defaults at the time of those experiments.

The global `cpml_partition_weight` stays **4.38**. The homogeneous Figure 5
recommendation uses an explicit **1**; it is not a calibrated choice for mixed
elastic/porous models. The current partitioner applies one coefficient to all
CPML materials and a common base cost to all ordinary cells. Neither number is
established as optimal for mixed physics; compare complete-step timing on the
intended model before selecting a case-specific override.

The defaults executable was rebuilt on Horizon with `ml reset`, rebuilding
both solver translation units and checking library resolution. SHA-256:
`206d308cb5242c4649cc91e8d785bef3b520ee7cdd1fd0b145e39ce31627364c`.
Jobs **95015** and **95016**, each one node/four ranks on `cuda-615-SC`,
completed in **13:39** and **10:51** with exit status zero. They compare
explicit production settings against omitted inputs and unset runtime flags:
42 CUDA-aware executions cover 19 short standalone/coupled fixtures and two
complete 0.6 s pulses; eight host-MPI executions repeat heterogeneous, VT
open-pore and both pulse pairs. Short fixtures include elastic, acoustic,
acoustic-elastic, porous-acoustic, porous-elastic, viscoelastic and nonlinear
solid paths, sealed/drained boundaries and FEM/VT geometry. These establish
execution and default selection, not new long-time or analytical validation.
All **425,572** consistent-mass solves passed their true-residual criterion.
The source, physical inputs, station sets and selected paths were checked.

Across all receivers, complete inviscid pulse differences were below
**0.000077% v**, **0.000013% p** and **0.0000046% q** relative L2 on both
backends. Viscous differences were at most **0.00143% v**, **0.0148% p** and
**0.678% q**. The weak viscous flux peaks near **4.13e-9 m/s**, with maximum
absolute difference **7.15e-11 m/s**. This retains the previously observed
weak-field sensitivity at `1e-6`, not bitwise reproducibility or a universal
observable-error bound; the earlier `1e-8` comparisons remain relevant when
that weak flux requires tighter agreement. Near-zero prearrival receiver
percentages are retained alongside absolute errors and all-receiver norms.

Both documentation builds and parser/rendered audits passed. Local-only
fixtures, source/build manifests, logs and individual station comparisons are
in the private validation archive, unavailable in
a clean checkout. No historical benchmark input or reference was changed.

## CPML history bypass and cached predictor products (2026-11-15)

Two algebraic optimizations extend the preceding production defaults:

- Shared CPML filters with zero feed bypass history loads/stores while retaining
  their direct stretch factor. Porous storage handles the old/new memory tails
  the same way, and full-Biot CPML skips Darcy history updates for exactly zero
  component drag. History allocation and active cascades are unchanged. This
  covers the shared elastic, acoustic and porous operators and interface factors;
  it does not change absorber equations or remove interface/geometry restrictions.
- Consistent-mass prediction retains verified `Ax` and `BAx` per solution history
  and solve role, allowing the starting residual and preconditioned residual to
  use the same linear combination as the predictor. The original RHS norm,
  consistent stabilization inverse, CG recurrence and final true-residual check
  remain. An apparently converged cached guess receives a fresh residual before
  acceptance. Damping-coefficient changes invalidate products. Each stored
  solution requires two extra vectors; setup reports the maximum per-rank cost.

`HERCULES_BIOT_MASS_PREDICTOR_CACHE` accepts `0|1`, consistently across ranks.
It defaults to **1 in double precision and 0 in single precision**; `0` recomputes
starting products. It is inactive for lumped mass or zero guesses. Single-precision
caching is experimental. Batched inertia/stabilization solves, one-round halo
assembly, periodic acceptance checks and lumped stabilization are not introduced.
The global CPML partition weight remains 4.38.

The final double-precision executable SHA-256 is
`ceb773353e940d18e13d450a989b854ca2a26851f68030ba6477dcbd4b55f154`.
Horizon builds use `ml reset` and resolved default-module libraries. The porous
CPML header also received a C-compatible typedef/C++ guard so C query utilities
can rebuild; its data layout is unchanged. The earlier performance executable
SHA-256 is `9484c990d642a18d631ba75b45748fc220a85eafac400939c2d67afd211ae9b3`.
Its numerical kernels are identical: the final source change only makes the
cache default precision-specific. Performance inputs select the cache explicitly.

### Correctness and waveform evidence

The CPU filter audit compares 24,064,000 values over double/float histories,
eight active-axis masks, non-unit inactive kappa, zero/nonzero alpha, strided
storage, component Darcy drag, and elastic/acoustic/porous/interface terms.
Maximum scaled discrepancy is **1.04e-13**; 5,131,013 values differ in bits.
Independent rational-filter algebra and C-header compilation checks pass.
Bitwise identity is not claimed.

Jobs **95098** and **95102**, each two nodes/eight ranks, exercise CUDA-aware and
host MPI. Their private audit manufactures projected histories on the actual
distributed operators; checks cached residual/product identities for all predictor
orders; deliberately supplies a falsely converged cache; tests coefficient-key
invalidation and zero RHS; and independently reevaluates every accepted residual.
All **206,488** audited accepted solves pass, including 102 manufactured solves.
Private audit code is excluded from production and timing binaries.

Jobs **95099**, **95100** and **95101** compare standalone/coupled physics,
FEM/VT geometry, sealed/drained/open/no-flow boundaries, heterogeneous and
anisotropic materials, both acoustic formulations, applicable viscoelastic and
nonlinear-solid paths, alternative operators/preconditioners/CG variants,
full/compact vectors and fused/separate or batched/unbatched communication paths.
Short fixtures establish execution and equivalence at their selected tolerance;
nonlinear examples remain synthetic and do not establish calibrated yielding.
Resolved 0.3 s FEM and VT interface pairs at `1e-8` differ by at most
**4.42e-6% / 2.98e-6%** all-receiver solid-velocity L2, respectively.

The complete 0.6 s inviscid pulse with both optimizations differs from the preceding
production baseline by **0.0000883% v / 0.0000113% p / 0.00000550% q** all-receiver
relative L2. The viscous pulse differs by **0.00134% v / 0.00728% p / 0.650% q**;
its tiny flux remains tolerance-sensitive at `1e-6`, as in the preceding repeated-run
evidence. Repeating cache off/on at `1e-10` reduces differences to
**8.04e-9% v / 2.75e-13% p / 4.76e-11% q** for the inviscid pulse and
**1.13e-8% v / 1.55e-7% p / 1.77e-9% q** for the viscous pulse. These are numerical
before/after comparisons, not new analytical or spatial-convergence results.

The completed suite contains **121 successful double-precision executions**,
including the earlier ten-case smoke test and four final default-selection runs.
The homogeneous **10 s / 128 m CPML** baseline/optimized pair also completes:
all-receiver differences are below **3.74e-10%** v/p/q L2, and successive one-second
station peaks decay throughout the record. The last-second v/q peaks are about
0.040% of the full-record peaks. This case uses lumped mass and no acoustic
interface, isolating the CPML shortcut; it cannot relax coupled-interface guards.

### Measured small-case cost

Job **95101**, one node/four GPUs, runs all variants in the same allocation.
These pulses use a **256 m cube, uniform 4 m cells and 32 m CPML**; they are much
smaller than the 4096 m production model. Complete Solver times in seconds:

| Variant | Inviscid | Viscous |
|---|---:|---:|
| Prior production baseline | 53.743 | 101.489 |
| CPML bypass only | 53.354 | 101.307 |
| Product cache only | 46.032 | 90.358 |
| Both | 46.050 | 90.069 |
| Baseline / both | **1.167x** | **1.127x** |

Caching reduces starting-product work while leaving mean iteration counts nearly
unchanged. The CPML shortcut has no resolved substantial whole-loop benefit in
this small case. Do not multiply these gains by earlier speedups or treat them as
32-node measurements. The completed production-scale comparison is recorded
below; the originally proposed 77 ms/step was a forecast, not a measurement.

### Completed 32-node production comparison: 95105 / 95113

Job **95113** completed with exit `0:0` in **41:37** on 32 nodes / 128 GPUs.
It advances the same inviscid full-Biot Figure 5 model through **1.5 s / 30,000
steps**: a 4096 m cube, uniform 4 m cells, 256 m CPML, 16 Hz source-rate peak,
50 microsecond timestep, stabilization `0.005`, consistent mass, tolerance
`1e-6`, quadratic prediction, element-inverse preconditioning and case-specific
CPML partition weight `1`. Its executable is the pinned earlier performance
binary `9484c990d642a18d631ba75b45748fc220a85eafac400939c2d67afd211ae9b3`.
The later acoustic solid-kernel skips are absent from these executables.

The first three 0.6 s cases in **95105** completed and were downloaded and
verified. At the maintainer's request to avoid duplicate completed combined
runs, its fourth case was intentionally cancelled after starting. It is excluded
from results, not classified as a solver failure. The first 0.6 s of 95113 supplies
the combined comparison. Both jobs were externally moved to `IO_testing`;
that reservation was preserved. No full CPML checkpoint/restart was attempted.

The following uses complete-step monitor wall-time differences from **0.0025
to 0.5975 s** (11,900 steps), excluding initialization and the first 50 steps:

| Variant | ms/step in common window | Baseline / variant |
|---|---:|---:|
| Prior production baseline | 98.210 | 1.000x |
| CPML bypass only | 90.857 | 1.081x |
| Product cache only | 81.319 | 1.208x |
| Both, prefix of full run | 81.303 | 1.208x |

The three 95105 cases share one allocation; 95113 uses another. Monitor times
have 0.1 s resolution. These are whole-step measurements, not exact per-stage
0.6 s counters. The combined result has no resolved additional whole-loop gain
over caching alone in this comparison. Do not multiply the individual gains.
All-receiver common-prefix differences from the baseline are at most
**0.000258% v / 0.0000414% p / 0.0000527% q** relative L2 across the three variants.

Full-duration timings use the same 30,000-step workload:

| Timer | Prior quadratic run 94833 (s) | Combined 95113 (s) |
|---|---:|---:|
| Complete Solver | 2892.591 | 2392.201 |
| Maximum-rank Compute Physics | 2880.261 | 2379.135 |
| Rank-zero inertia PCG | 1227.434 | 1021.523 |
| Rank-zero storage PCG, both solves | 745.472 | 591.479 |

Complete Solver cost falls from **96.42 to 79.74 ms/step**, a **1.209x speedup**
or **17.30% reduction in solver time**. Rank-zero PCG timers include waiting for
other ranks; they are not independent critical-path costs and must not be summed
with maximum-rank timers. The Slurm elapsed time includes allocation/launch work
outside the Solver timer. Against original consistent run **94379** (26330.141 s
Solver, tolerance `1e-10`), the cumulative Solver speedup is **11.007x**. This
historical ratio includes all intervening controls, tolerance and implementation
changes, across different allocations. It is not the isolated gain from caching
or CPML bypass. The new Solver time remains **2.390x** the historical lumped
run 94139 (1000.715 s); that older binary has neither latest optimization.

All **90,000** mass/storage solves pass their original fresh true-residual checks
at `1e-6`; worst inertia/storage residuals are **9.99979e-7 / 9.99981e-7**. Mean
iterations are **1.743 / 1.859**, nearly unchanged from 94833's 1.747 / 1.861.
Inertia reuses 29,999 starting residuals with zero initial refreshes; storage
reuses 59,998 with two initial refreshes. Final verification is retained on every
solve. Product caching adds at most **3.047 GiB per rank** in this case.

Eight complete receiver pairs contain 7,500 samples each, all finite; physical
inputs and source files match the earlier full runs. Maximum per-receiver
relative L2 differences over the complete native histories are:

| Field | From 94833 (%) | From tight `1e-10` run 94436 (%) |
|---|---:|---:|
| Solid velocity | 0.000523 | 0.004591 |
| Darcy flux | 0.0000207 | 0.002225 |
| Pore pressure | 0.0000202 | 0.002224 |

Against 94833, full-record FAS differences are at most **0.000338% v /
0.0000159% q / 0.0000158% p**. Amplitude-weighted phase RMS differences in
the 0–48 Hz band (components at least 1% of the reference spectral peak) remain
below **0.00000657 degrees** for v/q/p. In the arrival-based late window
`t >= 0.19 s + receiver_distance / slow_P_speed`, maximum changes normalized by
each full-record field peak are below **0.000271% v / 0.00000905% q /
0.00000566% p**. These are numerical differences, not physical error bounds.

The archived author/Karpfinger semi-analytical histories use the established
fixed source normalization and physical velocity time `t - dt/2`, with no fitted
gain, time shift or filtering. Worst reference L2 errors are **9.88861% v /
9.71150% q**, and 1–48 Hz FAS errors are **7.91048% / 7.92261%**. These are
essentially unchanged from 94833. Late v/q peaks remain **0.14120% / 0.04776%**
of the corresponding full reference peaks; existing late oscillations persist.
The author archive contains no independent pressure history. This validates
optimization equivalence for the tested model, not spatial convergence, a new
boundary-reflection bound or all coupled configurations.

### Precision limits and provenance

Single-precision short VT cache off/on comparisons showed approximately 5% v and
12% p differences, while repeated cache-disabled runs also changed by about
3.27% v and 6.61% p despite passing their algebraic tolerance. Thus these data do
not isolate caching as the cause, and do not establish single-precision waveform
accuracy. The conservative default remains recomputation. Job **95135** checks
explicit versus omitted precision-specific controls and their selected runtime
paths; it is not a single-precision accuracy certificate.

The single-precision inviscid pulse grows and aborts with caching disabled, both
with CPML and with dashpots. A corrected earlier-CPML control, **95136**, reproduces
both failures with both optimizations disabled. This broader precision limitation
is separate from the validated double-precision comparisons; its root cause and
repair are outside this change. The first attempted old-CPML control, **95118**,
was excluded because historical header timestamps reused a stale kernel object.
The corrected build forcibly recompiles all three affected translation units.
Failed/cancelled attempts and the correction are retained, not treated as passes.

Both documentation build targets and all 75 rendered/parser audits pass.
Exact inputs, source/binary hashes, logs, every-receiver waveform/FAS differences,
private audits and preparation/build scripts remain local-only, unavailable in a
clean checkout.
No existing physical benchmark input or analytical reference was changed.

## Batched mass solves and owner-complete rows — 2026-11-15 investigation

Two new execution experiments are **off by default**:
`HERCULES_BIOT_MASS_BATCH_SOLVES=1` pairs the first inertia solve with the
old-pressure stabilization inverse, retaining independent PCG recurrences and
fresh residual checks; `HERCULES_BIOT_MASS_ONE_ROUND=1` assembles complete uniform
Q1 rows at their owners and exchanges only their input columns during each
application. The latter retains the established two-exchange path for VT,
variable physical material/spacing, compact vectors and geometric constraints.
Single-reduction CG retains sequential solves. Neither switch changes physical
parameters, the stabilization operator, the RHS norm or the `1e-6` tolerance.
The second damped inertia stage and subsequent pressure update retain their order.

The comparison baseline includes the acoustic solid-kernel skips from `32398ef`
as well as the preceding CPML bypass and predictor cache. Every executable is
pinned separately; submitted inputs are immutable. Builds use Horizon's default
modules after `ml reset`, including both communication backends. The final
`joint-v02` executable has SHA-256
`bc42033288070db6c9e2a5998e470c8f9f78953881ed2d7da832bce1fb94b4b4`.

### Small-domain evidence

Job **95303** runs four variants of the same complete 0.6 s inviscid pulse on
**two nodes / eight GPUs**, in one allocation: 256 m cube, uniform 4 m cells,
32 m CPML, 50 microsecond timestep, 11,999 steps and tolerance `1e-6`.
These are not timings for the production-size Figure 5 model.

| Execution | Complete Solver (s) | Speedup from sequential |
|---|---:|---:|
| Sequential, two exchanges | 70.987488 | 1.000x |
| Batched solves only | 57.677068 | 1.231x |
| Owner-complete rows only | 51.240872 | 1.385x |
| Both | 45.036924 | 1.576x |

All solves pass the unchanged fresh residual check. With both enabled,
all-receiver aggregate L2 differences from sequential are **0.00008873% solid
velocity / 0.00001226% pore pressure / 0.000004563% Darcy flux**. At tolerance
`1e-10`, job **95310**'s inviscid comparison reduces these differences to
**8.46e-10% / 2.75e-13% / 3.68e-11%**, respectively. Its combined Solver time is
109.845 s against 197.512 s sequential. These measured ratios must not be
multiplied together or projected directly to 32 nodes.

Private setup-only audits in **95306** (CUDA-aware MPI) and **95309** (host MPI)
compare six manufactured distributed applications against the original operator:
mechanics, scalar and joint fields, each with mass and element-inverse patch.
Worst relative L2 discrepancy is below **5.24e-17**. Audit instrumentation is
absent from performance binaries. Job **95320** additionally completes a
single-rank run with empty peer schedules and the new joint profiling counters.

The complete viscous pulse is more sensitive at `1e-6`. In **95307**, the combined
variant differs from sequential by **0.001160% solid velocity / 0.009439% pressure /
0.5710% flux** aggregate L2. Flux peaks at only `4.13e-9 m/s`; maximum absolute
flux difference is `5.29e-11 m/s`. Against the tight sequential pulse, the two
`1e-6` variants have **0.2938% / 0.4338%** flux errors. A residual tolerance is
not a separate relative-error guarantee for each physical field. The completed
`1e-10` viscous comparisons in **95310** give the following aggregate flux
differences from the tight sequential solution:

| Execution | Complete Solver (s) | Flux difference (%) |
|---|---:|---:|
| Sequential | 337.582604 | 0 |
| Batched solves only | 283.086823 | 1.63e-9 |
| Owner-complete rows only | 228.157188 | 2.04e-9 |
| Both | 201.344498 | 1.57e-9 |

With both enabled, tight solid-velocity/pressure differences are **2.08e-9% /
8.71e-9%**, and all fresh residuals remain below `1e-10`. This supports algebraic
equivalence and identifies tolerance sensitivity in the weak-flux comparison;
it does not make the loose-tolerance results bit-identical or establish their
physical accuracy. The final-binary loose viscous pair in **95317** has a
**0.4141%** flux difference, consistent with that sensitivity.

Short comparisons in **95300**, **95307**, **95308** and **95309** cover the
refactor-disabled baseline, supported FEM/VT acoustic–porous interfaces,
pressure/velocity-potential acoustics, porous–solid coupling, anisotropic/sealed
permeability, physical damping, nonlinear-solid paths, and reduced/lumped or
nonporous controls. Solver-control comparisons cover both preconditioners,
predictors, standard/single-reduction CG, element/row operators, compact/full
vectors, fused/separate kernels and packet choices. Their scope is numerical
before/after regression, not a new physical validation or release of an interface
guard. In particular, **95312**'s longer interface fixtures select the intended
single-reduction fallback; **95318** repeats them with standard CG to exercise
actual batching. Final-binary pulse and zero-stabilization checks are **95317**.
With zero stabilization, no batched solve executes and there is exactly one
storage solve per step; the requested owner-row path remains active. Its completed
inviscid/viscous comparisons give **1.559x / 1.489x** Solver speedups; final
inviscid aggregate differences are **0.00007579% v / 0.00001013% p /
0.000003846% q**. The completed FEM standard-CG interface pair in 95318 gives
**1.208x** with differences **0.000003449% v / 0.0000008235% p /
0.0000001270% q**. The complete 0.3 s VT pair gives **1.181x** (501.866 s to
425.030 s), with differences **0.000002927% v / 0.0000004357% p /
0.00000006551% q**. Both interface pairs use tolerance `1e-8`, execute actual
batching, and explicitly fall back to two-exchange operators. Their distinct
11,999 / 59,999 step counts reflect the inherited FEM / VT timesteps. The
single-reduction controls in 95312 complete with both experiments inactive,
as intended. These **116** smaller cases pass executable identity, complete-record,
finite applicable fields and applicable fresh-residual checks. The four
production-size comparisons are reported below.
Single-precision accuracy, a new absorption bound and unsupported physical
combinations are not established by these double-precision comparisons.

### Full-size validation status

**95311** completes the 32-node, 1,073,741,824-element setup and 50-step checks.
The owner-row path is active and the residual checks pass. The combined Solver
is slightly slower during this short, quiet interval (6.847 s against 6.586 s);
this is not a resolved-pulse performance measurement.

**95322** completed on **32 nodes under `cuda-615-SC`** at 20:38:01 Central
on November 15, after 1:13:06 elapsed. It compares all four variants over 0.6 s
using the same 4096 m cube, uniform 4 m mesh, 256 m CPML, 16 Hz source,
50 microsecond timestep, `gamma=0.005`, tolerance `1e-6` and partition weight 1.
All variants share one allocation and the same executable.

| Execution | Complete Solver (s) | Solver speedup | Change in Solver time |
|---|---:|---:|---:|
| Sequential | 998.923693 | 1.000x | baseline |
| Batched solves only | 967.543097 | 1.032x | 3.14% reduction |
| Owner-complete rows only | 1048.376412 | 0.953x | 4.95% increase |
| Both | 1014.958098 | 0.984x | 1.61% increase |

The small-domain gains do **not** transfer directly to this scale. Batching is
only modestly faster; owner-complete rows regress here, both alone and combined.
The 32-node result does not justify making either experiment a default.

All four complete records pass source/input equivalence, binary identity,
finite applicable station fields and fresh-residual checks. Across the three
variants, worst per-receiver v/p/q waveform L2 differences from sequential are
below **0.000726%**, FAS differences below **0.000606%**, and amplitude-weighted
phase RMS below **0.000073 degrees** (0–80 Hz, components above 1% of the spectral
peak). Differences after 0.35 s remain below **0.000806%** of each corresponding
full-record field peak. This fixed 0.6 s comparison window is not the full
arrival-based late-time test or an independent physical error bound.

**95347** completed the full **1.5 s** confirmation on 32 CUDA-reservation nodes
at **22:17:41 Central on November 15**, after **1:23:11** elapsed. Batching runs
first and the current sequential baseline second in the same allocation,
reversing the 95322 order. The immutable `joint-v02` inputs preserve the same
physics; only `HERCULES_BIOT_MASS_BATCH_SOLVES` differs.

| Full 1.5 s execution | Complete Solver (s) | Average Solver (ms/step) |
|---|---:|---:|
| Sequential | 2440.030108 | 81.334 |
| Batched solves only | 2365.970210 | 78.866 |

This is **1.0313x faster**, or **3.04% less Solver time**, consistent with the
modest gain in 95322. It is not the larger small-domain gain and does not justify
enabling owner-complete rows. Keep both experimental switches off by default;
batching remains an optional, workload-dependent choice. The existing cache,
CPML and other established execution defaults are unchanged.

Both complete 30,000-step cases contain eight solid/fluid station pairs with
7,500 finite samples each. Every one of the 30,000 inertia and 60,000 storage
solves passes its original fresh true-residual check at `1e-6`. Mean iterations
are 1.743 / 1.859 for both executions. Worst inertia/storage residuals are
`9.99994965e-7 / 9.99971975e-7` batched and
`9.9998758e-7 / 9.99985165e-7` sequential.

Maximum per-receiver differences from the matched sequential baseline are:

| Field | Complete waveform L2 (%) | FAS, 0–48 Hz (%) | Late difference / full baseline peak (%) |
|---|---:|---:|---:|
| Solid velocity | 0.0004011 | 0.0000102 | 0.0002064 |
| Darcy flux | 0.0000180 | 0.0000124 | 0.00000706 |
| Pore pressure | 0.0000195 | 0.0000118 | 0.00000607 |

The native, untapered full-record FAS comparison uses 0–48 Hz; amplitude-weighted
phase RMS is below **0.00000644 degrees** for these three fields, retaining
components at least 1% of the baseline spectral peak. The late window starts at
`0.19 s + receiver_distance / slow_P_speed`, with velocity time `t - dt/2`.
These numerical differences are not physical reflection-error bounds.

Batched results also retain the earlier production and tight numerical agreement:
worst v/q/p waveform differences from 95113 are **0.0004454% / 0.00001592% /
0.00001942%**, and from tight run 94436 are **0.004598% / 0.002222% /
0.002221%**. The author-reference errors remain about **9.89% velocity /
9.71% flux**; optimization equivalence does not establish mesh convergence.

### Removing pressure stabilization

At the maintainer's request, **95354** repeats the full 1.5 s, 32-node case under
`IO_testing` with `poroelastic_pressure_stabilization=0` instead of `0.005`.
All other parsed physical/source inputs and the executable are unchanged.
The run completed successfully in **35:28** Slurm elapsed; complete Solver time
is **2024.548856 s**. Batching remains selected but does not execute, because
there is no stabilization inverse to pair with inertia. The run performs
30,000 inertia and 30,000 storage solves, rather than 30,000 and 60,000.
All residuals pass `1e-6` and complete station records are finite. Removing
stabilization changes the discrete numerical dissipation, so this timing is not
an optimization-only speedup; the allocation also differs from 95347.

Finite completion does not imply acceptable waveforms. The same fixed
author/Karpfinger reference, source normalization `7.025e12` and physical
velocity time `t - dt/2` give the following worst-receiver results. No gain,
time shift or filter is fitted. The late window starts at
`t >= 0.19 s + receiver_distance / slow_P_speed`; late peaks are radial
amplitudes normalized by the corresponding full reference pulse peak.

| Quantity | Batched `gamma=0.005`, 95347 | `gamma=0`, 95354 |
|---|---:|---:|
| Solid-velocity waveform L2 error (%) | 9.88861 | 9.08551 |
| Darcy-flux waveform L2 error (%) | 9.71150 | 173.86423 |
| Solid-velocity FAS error, 1–48 Hz (%) | 7.91048 | 6.13551 |
| Darcy-flux FAS error, 1–48 Hz (%) | 7.92261 | 97.73142 |
| Late solid-velocity peak / reference peak (%) | 0.14125 | 4.25351 |
| Late Darcy-flux peak / reference peak (%) | 0.047757 | 31.82761 |

The slightly smaller solid-velocity L2 error hides a much larger fluid error.
There is no independent pressure history in the author archive. Native
arrival-window spectra put the dominant gamma-zero late flux around **18–33 Hz**;
less than **0.045%** of its windowed spectral power lies above 80 Hz. The
80 Hz plots therefore include the principal oscillations. These 1.5 s records
do not establish unbounded growth or uniquely distinguish discretization and
absorbing-boundary effects.

The short follow-up **95370** checks tolerance sensitivity on two CUDA nodes
using complete 0.6 s pulses in a 256 m cube with 4 m cells. The gamma-zero pair
changes only the PCG tolerance from `1e-6` to `1e-10`. A second comparison uses
the tight stabilized pulse from 95310; its numerical implementation is unchanged
between the pinned `joint-v01` and `joint-v02` binaries.

| All-receiver aggregate L2 difference (%) | Solid velocity | Pore pressure | Darcy flux |
|---|---:|---:|---:|
| Gamma zero: tolerance `1e-6` vs `1e-10` | 0.002321 | 0.000500 | 0.000258 |
| At `1e-10`: gamma zero vs `0.005` | 6.36880 | 112.23283 | 124.23717 |

Thus loose algebraic tolerance is not the main cause of the gamma-zero change
in this small control. This is not a full-domain tight-tolerance certificate or
a unique diagnosis of the oscillations. Retain `gamma=0.005` for the current
Figure 5 setup; no universal stabilization calibration or production default is
changed by this experiment.

Direct comparison with the completed **95347 sequential baseline** permits
only the declared stabilization change in parsed inputs and the inert batching
selection in runtime controls; source files are identical. Worst per-receiver
gamma-zero v/q/p waveform differences are **10.17637% / 178.64651% /
296.97068%**. Native 0–48 Hz FAS differences are **8.31126% / 109.09149% /
213.29733%**, and weighted phase RMS differences reach **3.11414 / 82.45256 /
87.98614 degrees**, respectively. Arrival-window differences normalized by
full baseline component peaks are **4.53380% / 80.29889% / 49.09466%**.
Those baseline-normalized differences are distinct from the radial
author-reference late amplitudes in the table above.

All **125** cases in this investigation are complete and pass their applicable
executable/input identity, record-length, finite-field and residual audits.
This includes gamma-zero cases with poor physical-reference agreement: passing
these consistency audits is not an accuracy certificate. The final seven-page
local PDF includes performance, all-receiver summary metrics, methods and
representative waveform, FAS and late-window figures. Both documentation build
targets and the required lint, TypeScript and rendered/parser audits pass.

Detailed manifests, inputs, logs, all-receiver waveform/FAS comparisons and local
helpers are retained in the private validation archive and
are unavailable in a clean checkout. No established physical benchmark input or
analytical reference is modified by this investigation.

## Homogeneous Wolf reruns with current defaults (2026-11-16)

The requested **homogeneous Figures 5 and 7 comparisons are complete**. The
inviscid Figure 5 record reuses the sequential case in **95347**; the three new
jobs completed with exit `0:0`. Completion and algebraic convergence are checked
separately from the reference errors below. The free-surface pressure remains
less accurate than the solid motion near the drained boundary.

| Job | Homogeneous setup | Uniform spacing | Record / timestep | Nodes | Slurm elapsed | Complete Solver time |
|---|---|---:|---:|---:|---:|---:|
| 95347, sequential | Figure 5, inviscid full space | 4 m | 1.5 s / 50 microseconds | 32 | Reused case within allocation | 2440.030 s |
| 95427 | Figure 7, drained free surface | 8 m | 2 s / 100 microseconds | 4 | 22:04 | 1248.023 s |
| 95430 | Figure 7, drained free surface | 4 m | 2 s / 50 microseconds | 32 | 56:09 | 3269.919 s |
| 95433 | Figure 5, viscous full space | 4 m | 1.5 s / 50 microseconds | 32 | 1:04:55 | 3792.249 s |

All allocations use four GPU ranks per node. 95427 and 95430 ran on
`cuda-615-SC`; 95433 was submitted there and externally reassigned to
`IO_testing`. The reassignment was preserved. These timings compare different
physics, record lengths and allocations; they are not optimization-only speedups.

### Models and execution controls

Each model is a **4096 m cube**, with **256 m CPML** at the four lateral faces
and bottom, and no top CPML. There is no graded mesh or source refinement.
The 8 m and 4 m meshes contain **134,217,728** and **1,073,741,824** cells,
respectively, with no hanging nodes. Physical porous Q1 mass/storage are
consistent; CPML uses the existing lumped update. Default CPML profiles use
reflection target `0.01`, polynomial order `2`, maximum kappa `5`, maximum alpha
`0.5 Hz`, alpha-floor fraction `0.1` and zero MPML cross damping. Thickness,
`poroelastic_pressure_stabilization=0.005` and `cpml_partition_weight=1` are
explicit benchmark choices; the latter two global defaults remain `0.05` and
`4.38`.

The new runs omit mass-performance option lines and environment overrides to
exercise tolerance `1e-6`, maximum 200 iterations, element-inverse
preconditioning, quadratic prediction, standard CG, automatic rows, fused
passes, full vectors, combined solid/fluid packets and cached predictor
products. Aggregate packet kernels, batched inertia/stabilization and
owner-complete one-round rows remain off. Every solve retains fresh residual
verification. The reused inviscid case explicitly selects the same performance
options. All four cases pass finite-field, record-length, mesh, receiver,
source/executable identity and residual audits.

| Case | Inertia / storage solves | Mean inertia / storage iterations | Maximum inertia / storage iterations |
|---|---:|---:|---:|
| Figure 5 inviscid | 30,000 / 60,000 | 1.743 / 1.859 | 9 / 9 |
| Figure 5 viscous | 60,000 / 60,000 | 2.202 / 2.485 | 13 / 9 |
| Figure 7, 8 m | 20,000 / 40,000 | 1.763 / 1.942 | 9 / 9 |
| Figure 7, 4 m | 40,000 / 80,000 | 1.710 / 1.977 | 9 / 9 |

Worst fresh relative residuals are below `1e-6` in every case. The viscous run
uses two Darcy SDIRK inertia stages; its additional solve and different
iteration counts explain why its workload differs from the inviscid run.
These checks do not bound continuum error or independently establish tolerance
insensitivity for the weak viscous fluid signal.

The pinned `joint-v02` executable has SHA-256
`bc42033288070db6c9e2a5998e470c8f9f78953881ed2d7da832bce1fb94b4b4`.
Its recorded source hashes matched the workspace at preparation. The later
parser-only commit `9dcbd95` makes consistent mass the conditional default for
active full Biot; these runs explicitly selected consistent mass and predate
that parser change. No new CUDA/MPI executable was built for that change, and
these jobs do not test omission of the mass-matrix key. Default modules after
`ml reset` and retained GSL/PROJ paths resolve the pinned executable's libraries.

### Figure 5: full-space solid velocity and Darcy flux

The source rate has **16 Hz central frequency**, distinct from the approximately
**48.5436 Hz** mesh-control input. Its fixed position is `(2049,2049,3073)` m;
the four axial and four diagonal reference receivers are unchanged. The remote
top surface has a conservative earliest-return lower bound of **1.5517 s**,
beyond the 1.5 s record. This is a finite-domain full-space approximation, not
an assertion of zero CPML reflection.

The viscous case changes only viscosity among the physical inputs, from `0` to
`0.001 Pa s`; source files, receivers, mesh, timestep and stabilization are
identical. Each case uses its own supplied Karpfinger reference archive. The
fixed source normalization is `7.025e12`, and physical velocity/flux time is
`t - dt/2`. No fitted gains, shifts or filtering are applied. Darcy flux is the
volume-weighted relative-fluid velocity, not that quantity divided by porosity.

| Receiver | Inviscid velocity L2 (%) | Inviscid flux L2 (%) | Viscous velocity L2 (%) | Viscous flux L2 (%) |
|---|---:|---:|---:|---:|
| x1 | 9.03489 | 9.38466 | 3.50521 | 4.85883 |
| x2 | 6.32884 | 5.98741 | 1.86903 | 2.54731 |
| x3 | 7.09047 | 6.36302 | 4.85849 | 5.93528 |
| x4 | 9.88861 | 9.71150 | 6.77225 | 8.59158 |
| d1 | 6.79223 | 7.96727 | 0.70498 | 1.15493 |
| d2 | 2.42747 | 2.57450 | 0.92319 | 1.28111 |
| d3 | 4.77914 | 5.36117 | 1.75997 | 2.14196 |
| d4 | 4.59348 | 5.11735 | 2.07271 | 2.40212 |

L2 errors use the full three-component history at author sample times over
`0.001–1.499 s`. Radial FAS uses a Tukey `0.1` taper and `1–48 Hz`; phase RMS
is weighted by reference spectral power where amplitude exceeds 1% of its peak.
Plots extend to 80 Hz but reported FAS scores use the stated band.

| Worst receiver metric | Inviscid velocity | Inviscid flux | Viscous velocity | Viscous flux |
|---|---:|---:|---:|---:|
| Waveform L2 (%) | 9.88861 | 9.71150 | 6.77225 | 8.59158 |
| Radial FAS L2, 1–48 Hz (%) | 7.91048 | 7.92262 | 0.90412 | 1.14193 |
| Weighted phase RMS (degrees) | 4.00849 | 4.29688 | 3.84399 | 4.90741 |

The viscous amplitudes agree more closely than the waveform scores alone
suggest; phase differences remain. This comparison does not uniquely separate
spatial, temporal, source-discretization, boundary or algebraic errors. Small
absolute viscous flux does not establish relative accuracy. No matching
full-domain tighter-tolerance viscous control was run in this refresh, and the
Karpfinger archive contains no independent pressure history.

### Figure 7: free-surface refinement and reference checks

Both runs retain Table 7's `z>0` material, the original **5 Hz** source rate at
**500 m depth**, and zero applied total/fluid traction at the top. The drained
condition is imposed naturally in the weak formulation: **surface pressure
DOFs are not strongly clamped to zero**. The top receiver measures the resulting
pressure-trace error; no stronger full-Biot pressure-clamp option was used.
The reduced-formulation pressure penalty does not change this full-Biot path.

Reference receivers have horizontal offsets `(500,400)` m and depths **0.5 m**
and **500 m**. Additional 1 m and 0 m receivers are diagnostics, without matching
archived subsurface references at 1 m. The source remains at `(2049,2049,500)` m
in both meshes; its depth lies on an element face at 4 m. Timestep halves with
spacing, and source position within a cell changes. This is joint spatial/time
refinement, not an isolated spatial convergence-rate measurement. Mesh-control
frequencies are approximately **28.8523 Hz / 57.7045 Hz**, not the 5 Hz source
central frequency.

The supplied Gar6More3D files contain displacement. Comparisons preserve the
fixed author scale `1e10/(2*pi^2*5^2)` and time convention `+0.001 s`. Velocity
is derived with a centered NumPy gradient because the original conversion
helper is absent; identical derivative processing is not claimed. Archived
SeisSol velocity and pressure supply a **secondary numerical comparison**, not
an independent analytical pressure solution.

| Reference / field | Depth | 8 m L2 error (%) | 4 m L2 error (%) |
|---|---:|---:|---:|
| Gar6More3D displacement | 0.5 m | 1.88860 | 0.74436 |
| Gar6More3D derived velocity | 0.5 m | 2.42617 | 0.91080 |
| Gar6More3D displacement | 500 m | 3.79821 | 3.52855 |
| Gar6More3D derived velocity | 500 m | 4.47020 | 4.12712 |
| SeisSol velocity | 0.5 m | 2.46063 | 0.98667 |
| SeisSol velocity | 500 m | 1.51790 | 0.82084 |
| SeisSol pressure | 0.5 m | 52.40021 | 16.35827 |
| SeisSol pressure | 500 m | 0.93088 | 0.28701 |

The full-window metrics retain the documented reference quadrature artifact
near 1.5 s. Restricting the 4 m derived-velocity comparison to before 1.3 s gives
**0.90474% / 4.03313%** at the shallow/deep receivers, so removing that portion
does not eliminate the deep discrepancy. At 4 m, displacement/velocity FAS
errors over `1–12.5 Hz` are **0.51857% / 0.63367%** shallow and
**0.63819% / 0.69463%** deep.

With the same fixed conventions, the archived SeisSol velocity itself differs
from the differentiated Gar6More3D reference by **0.23890% shallow** and
**3.84574% deep**. This contextualizes the persistent deep discrepancy; it does
not identify which reference or processing is responsible, remove the primary
reference error, or prove Hercules correct. No fitted alignment is introduced.

At the exact surface, pressure peak decreases from **0.473057 Pa** to
**0.170401 Pa**, a **2.776x** reduction; RMS decreases from **0.082493 Pa** to
**0.029680 Pa**. The 4 m peak is **0.52790%** of the deep-station pressure peak.
The 0.5 m pressure difference also decreases substantially but remains
**16.36%** against SeisSol. Near-surface pressure therefore remains a validation
gap despite the much closer solid motion. Pressure histories are retained in
the local PDF, including a common-axis refinement comparison and arrival zoom.
No independent Figure 7 pressure/flux reference is available in these archives.

### Timestep, stabilization and scope limits

The generic roughly 1 ms startup “critical timestep” is not the final
consistent-mass check. The latter includes a hard-coded wave ceiling
`0.08*h/c_fast` (`0.02*h/c_fast` for VT), introduced with the original
consistent-mass implementation. No derivation or controlled timestep sweep
justifying those coefficients was found. They are **software-enforced heuristic
ceilings, not measured critical timesteps**. For these 4 m cases they are
89.1283 microseconds (Figure 5) and 75.3499 microseconds (Figure 7), above the
chosen 50 microseconds. The Figure 7 8 m ceiling is 150.6999 microseconds.
Passing those checks does not establish that 50 microseconds is necessary or
optimal. The guard and timesteps were not changed in that refresh. The
[subsequent timestep study](#homogeneous-timestep-investigation-2026-11-16)
below measures stability and temporal accuracy separately using private builds;
the production guard remains unchanged.

Default stabilization `0.05` was used in early graded/lumped Wolf jobs
93629–93631. Those runs do not provide a matched `0.05` versus `0.005` comparison
for the present uniform consistent-mass setup; the
[subsequent matched gamma study](#matched-stabilization-investigation-2026-11-16)
below supplies that comparison. This refresh does not recalibrate
the global stabilization default, establish a boundary-reflection bound, or
certify elastic, acoustic or heterogeneous coupling from homogeneous data.

Figure 6 and LOHp (Figure 8 onward) remain outside the homogeneous-only scope.
They need material-dependent hydraulic properties; LOHp also needs its
prescribed combined stress/pressure line source. All required homogeneous
references are already downloaded; no further download was needed.

The completed local report is
`Hercules_Wolf_Homogeneous_Validation_2026-11-16.pdf`. Immutable inputs, logs,
source/executable hashes, reservation history, all-receiver metrics and plots
are retained in the private validation archive. These
investigation artifacts/helpers are local-only and unavailable in a clean
checkout. Submitted inputs and existing benchmark/reference outputs are
unchanged.

## Homogeneous timestep investigation (2026-11-16)

This follow-up isolates the timestep in the **inviscid Figure 5** case before
examining the Figure 7 surface condition and then stabilization. It retains the
uniform 4 m mesh, 4096 m cube, original 16 Hz source rate and eight receivers,
256 m default-profile CPML, consistent physical inertia/storage, stabilization
`0.005`, mass tolerance `1e-6` and the previously selected performance defaults.
The top remains mechanically free and hydraulically undrained, without top
CPML. It lies beyond the original 1.5 s reference window, but can reflect later.

### Private guard relaxation and independent stability estimate

The experiment uses an isolated, hashed `dtguard-v02` executable. It relaxes the
two native consistent-mass wave ceilings to permit the sweep; it does not change
the physical operators, VT/lumped ceilings, dissipation guard or CPML updates.
The production source and input interface are unchanged. Consequently the larger
steps below are **not accepted simply by editing a production input file**; no
public bypass or new default is introduced by this investigation. A 50-microsecond
control checks the private build against the pinned original executable.

An independent periodic uniform-Q1 recurrence gives a limiting step of
**501.466 microseconds** for this material and stabilization. The limiting
three-axis corner mode agrees with

```math
\Delta t_{\mathrm{corner}}
=h\sqrt{\frac{\rho-\rho_f^2/m}{3(\lambda+4\mu)}},
\qquad m=\frac{\mathcal{T}\rho_f}{\phi}.
```

Here `lambda` and `mu` are drained skeleton moduli, `rho` is mixture density,
and `m` is Biot relative-fluid inertia density. The formula describes this
discrete corner mode; it is not a general critical-step formula for arbitrary
materials, stabilization, meshes or boundaries. The symbol excludes CPML,
free surfaces, solver tolerance and MPI constraints.

Job **95568** supplies a finite-domain stress test: a uniform 4 m, 256 m cube
with 32 m CPML. Steps of 450, 480 and 500 microseconds complete 10 s with finite,
decaying station records. At **510 microseconds**, records grow exponentially
before an algebraic norm failure. The measured amplitude growth, approximately
685–704 per second in the fitted interval, is consistent with the periodic
prediction of about 722 per second; finite-box modes need not match exactly.
Thus the observed onset is bracketed between **500 and 510 microseconds** in
this stress model. An earlier startup rejection at 300 microseconds came from
a remaining software guard, not a measured physical instability. Failed MPI
cleanup was cancelled after preserving the evidence.

### Full Figure 5 temporal differences and cost

Job **95569** runs the full 1.5 s comparisons on 32 Horizon nodes. The reference
for temporal differences is the completed 50-microsecond case in **95347**,
whose Solver timer is 2440.030 s. All completed cases pass input/source/binary,
uniform-mesh, finite-station and per-solve true-residual audits. Table entries
are the **worst receiver relative L2 changes from that numerical baseline**,
not errors against an analytical solution.

| Step (microseconds) | Displacement (%) | Solid velocity (%) | Darcy flux (%) | Pressure (%) | Solver (s) | Solver speedup |
|---|---:|---:|---:|---:|---:|---:|
| 100 | 0.06506 | 0.10412 | 0.12121 | 0.12114 | 1449.277 | 1.684x |
| 200 | 0.29696 | 0.44235 | 0.54437 | 0.54440 | 814.117 | 2.997x |
| 300 | 0.66507 | 0.98373 | 1.20924 | 1.20953 | 635.158 | 3.842x |
| 400 | 1.16929 | 1.72338 | 2.11588 | 2.11638 | 537.738 | 4.538x |
| 450 | 1.47239 | 2.16664 | 2.65961 | 2.66024 | 477.000 | 5.115x |

**200 microseconds is the working accuracy choice for this study**, the largest
tested step keeping all four fields within 1% of the 50-microsecond baseline.
The 1% criterion is a study target, not a physical law. This step has
`c_fast*dt/h = 0.17952`, about 40% of the periodic limiting step. The stable
threshold and the accuracy margin must not be conflated with one universal
replacement for the production `0.08` coefficient.

At 200 microseconds, worst solid-velocity, Darcy-flux and pressure FAS changes
from the baseline over 1–48 Hz are **0.19722%, 0.01640% and 0.00768%**; weighted
phase changes are **0.22785, 0.31157 and 0.31175 degrees**. Original source
samples, source normalization and receiver positions are unchanged. Displacement
and pressure use time `t`; velocity and flux use `t-dt/2`. No gain or time shift
is fitted.

Against the fixed author references, the 200-microsecond worst velocity/flux
L2 errors are **9.77081% / 9.45229%** (Karpfinger), and pressure is **8.97097%**
against the archived SeisSol numerical record. The 50-microsecond values are
**9.88861% / 9.71150% / 9.24811%**. Slightly smaller reference errors at larger
steps can result from cancellation of temporal and spatial errors; they do not
establish improved spatial convergence. SeisSol pressure is a secondary
numerical comparison, not an independent analytical pressure solution.

### Full-domain ten-second control

Job **95602**, 32 nodes on `IO_testing`, completes the full 4096 m model at
400 microseconds through **10 s**. Slurm elapsed is **1:13:48**, Solver time
4323.720 s. All 25,000 inertia and 50,000 storage solves meet the `1e-6` true
residual threshold; all station records are finite. Its first 1.5 s agrees with
the separate CUDA short run to at most **0.00134%** across displacement,
velocity, flux and pressure. No late exponential growth is observed at the
stations. This is a configuration-specific stability observation, not a global
energy proof or certification of other interfaces.

Later arrivals must be interpreted using the actual boundary geometry. The top
is reflective, not absorbed. Image-source slow-P top-return guides at 5.402 and
5.918 s lie close to observed pressure peaks at 5.414 and 5.930 s. This supports
a reflected-pulse interpretation, without uniquely isolating mode conversion
or every returning path. Pressure RMS across the station components declines
from **1.05051 Pa** at 5–6 s to **0.56370**, **0.03265**, **0.001383** and
**0.00008447 Pa** in successive one-second bins. Reference scoring remains
restricted to the original 1.5 s window.

### Figure 7 timestep qualification

Figure 7 retains **its own Table 7 material, original 5 Hz source at 500 m depth
and author receivers**, rather than substituting a shallow Figure 5 model.
Its periodic-Q1 limiting step at 4 m is **432.461 microseconds**, different from
Figure 5. The natural drained-boundary controls qualify 200 microseconds against
the original smaller steps before changing the boundary treatment.

| Job | Spacing | Baseline step (microseconds) | New step (microseconds) | Solver time (s) | Solver speedup |
|---|---:|---:|---:|---:|---:|
| 95631 | 4 m | 50 | 200 | 1018.558 | 3.210x |
| 95632 | 8 m | 100 | 200 | 692.608 | 1.802x |

Both complete the original 2 s record with finite fields and passing true
residual checks. Across the 4 m subsurface receivers, changes in displacement,
velocity, Darcy flux and pressure are below **0.014%**; the small nonzero surface
pressure trace changes by **0.151%**. Corresponding 8 m changes are below
**0.026%** subsurface and **0.094%** at the surface. The 4 m surface pressure
peak remains **0.17036 Pa**, and its 0.5 m pressure error against SeisSol remains
**16.36%**. The surface discrepancy is therefore not explained by this timestep
change.

These results apply to homogeneous inviscid, native Hex8, consistent-mass cases.
They do not validate a larger production ceiling for viscosity, changing
stabilization, heterogeneous/elastic/acoustic interfaces, hanging nodes,
FEM/VT terrain, restart or nonlinear/attenuating skeletons. Production guards
and defaults remain unchanged. Local-only inputs, pinned builds, hashes,
completed records and analysis are in the private validation archive; these
artifacts are not
available in a clean checkout. The subsequent boundary and gamma assessments
must retain their own scope and reference qualifications.

## Private drained-surface experiment (2026-11-16)

After the separate Figure 7 timestep qualification, this experiment retains
200 microseconds, stabilization `0.005`, Table 7 material, the original 5 Hz
source, uniform meshes and the author receiver positions. It compares the
production natural zero-fluid-traction condition with **strong zero excess
surface pressure** in a private executable. No production option or default
is added.

The private implementation removes owned top pressure degrees of freedom from
both consistent scalar solves, including the inverse used by stabilization.
The physical pressure and stabilization auxiliary therefore share zero surface
trace. Zero initial pressure plus zero surface increments preserves `p=0`;
normal owner/ghost communication preserves that value. This is a restriction of
the trial/test spaces, not a post-step pressure clamp or artificial loss of
storage capacity. The auxiliary trace changes the stabilization boundary closure
and is part of the experiment, rather than an assumed independent consequence
of the physical condition.

Solid and relative-fluid mechanical unknowns remain active. Normal Darcy flow
is allowed: the strong/natural surface-normal peak-flux ratios over the reference
window are **1.00356** at 8 m and **1.00086** at 4 m. The boundary has not become
sealed. Setup reports **263,169** and **1,050,625** constrained pressure unknowns,
respectively; the surface diagnostic pressure remains exactly zero throughout
each completed record.

### Reference-window accuracy

All rows below use the original **2 s** reference window. Pressure errors are
against archived **SeisSol numerical data**, not an independent analytical
pressure solution. Displacement/velocity use the fixed Gar6More3D conventions
already recorded above. The 1 m and surface stations are diagnostics only.

| Mesh / treatment | Job | Shallow displacement L2 (%) | Shallow velocity L2 (%) | Shallow pressure L2 (%) | Deep pressure L2 (%) |
|---|---:|---:|---:|---:|---:|
| 8 m natural | 95632 | 1.88982 | 2.42840 | 52.39279 | 0.92976 |
| 8 m strong | 95669 | 1.78556 | 2.29131 | **1.92036** | 0.80975 |
| 4 m natural | 95631 | 0.74782 | 0.91617 | 16.36272 | 0.29151 |
| 4 m strong | 95670 | 0.74133 | 0.90704 | **0.51874** | 0.30308 |

Shallow means 0.5 m depth and deep means 500 m. Strong-condition shallow
pressure FAS errors over 1–12.5 Hz are **1.86329%** at 8 m and **0.43453%** at
4 m. The improvement persists with refinement while preserving the physical
source position; changing that position within the source cell prevents claiming
an isolated spatial convergence order from these two meshes.

At 4 m, strong versus natural displacement/velocity differences are below
**0.03447% / 0.04454%** L2 across the stations. The largest Darcy-flux change is
**0.35525%**, at the surface. These are differences between schemes, not errors
against an independent flux reference. The deep pressure difference between
schemes is only **0.03773%**. The original deep-reference discrepancy remains;
this boundary change does not explain or resolve it.

### Long-duration check and supported scope

The 8 m strong case **95669** completes **10 s**, with Slurm elapsed **1:12:39**
and Solver time **4281.563 s**. All 50,000 inertia and 100,000 storage solves
meet the `1e-6` true-residual threshold; all station fields remain finite. No
late exponential growth is observed at the stations. Later pulses prevent a
monotonic-decay claim: pressure RMS is **0.00126651 Pa** at 3–4 s,
**0.00025393 Pa** at 7–8 s, **0.00012939 Pa** at 8–9 s and
**0.00003030 Pa** at 9–10 s. The matched natural 10 s control **95685** also
completes, in **1:12:41** with Solver time **4283.433 s**, finite fields and all
50,000/100,000 inertia/storage residual checks passing. Its pressure RMS over
the same windows is **0.00126793**, **0.00025286**, **0.00012841** and
**0.00002992 Pa**. The later pulses occur in both treatments and do not show an
exponential growth mechanism introduced by the strong constraint. Across the
one-second windows from 2 to 10 s, maximum strong-minus-natural velocity, flux
and pressure differences are **0.000397%, 0.000141% and 0.000135%** of each
field's natural first-2s peak. These small differences do not identify every
returning path or prove global energy stability.

The 4 m strong case **95670** covers **2 s**, not 10 s. It completes in **18:52**,
with Solver time **1033.062 s** and all 10,000/20,000 inertia/storage solves
passing. Its different allocation from the natural run prevents treating the
small timing difference as an isolated boundary-overhead measurement. The
coarser long test does not establish fine-mesh long-time stability or global
energy conservation.

The private build is limited to flat, all-porous native Hex8, consistent mass,
full vectors, sequential solves and the usual two-round mass exchange. It
rejects topography, acoustic/mixed nonporous cells and hanging nodes. Arbitrary
initial pressure, restart, FEM/VT terrain, elastic/acoustic interfaces,
alternative mass communication/solver paths, viscosity, attenuation and nonlinear
constitutive response are not certified by these tests. Extending this method
to those paths requires explicit constraint, initialization, time-integration,
GPU/MPI and output assessments; homogeneous pressure agreement is insufficient.

The pinned `surface-strong-v01` executable differs from `dtguard-v02` only in
the private pressure-space restriction. Its SHA-256 is
`79614b76bd66b9234dc4f80cbd976dfa37bc4ffb3f88d21f425d28a233cb54ce`.
Both builds, source hashes, input hashes and complete station records remain
local investigation artifacts. Production `poroelastic_surface_pressure_bc`
behavior remains unchanged; the reduced-formulation pressure penalty has not
been made effective in full Biot by this experiment.

## Matched stabilization investigation (2026-11-16)

After the timestep and surface studies, the Figure 5 inviscid gamma comparison
retains the uniform 4 m mesh, 4096 m cube, 256 m CPML, original 16 Hz source,
all eight author receivers, consistent physical inertia/storage and optimized
defaults. Every case uses **50 microseconds** and a **1.5 s** record with the
same pinned `joint-v02` executable; only `poroelastic_pressure_stabilization`
changes. At gamma `0.05`, the unchanged consistent stabilization guard is
**58.45623 microseconds**, so a common 200-microsecond step would be rejected.
This sweep therefore isolates gamma without relaxing that dissipation guard.

Archived jobs **95354** (`0`) and **95347** (`0.005`) pass the same input,
source, executable, mesh, station and residual audits. The zero-gamma run's
selected batching control is inactive because there is no stabilization solve;
its timer and solve counts confirm that distinction. Velocity and Darcy-flux
errors below use the fixed Karpfinger reference; pressure uses the secondary
SeisSol numerical records. No source scale, time alignment or filter is fitted.

| Gamma | Job | Worst velocity L2 (%) | Worst Darcy-flux L2 (%) | Worst pressure L2 (%) | Late flux peak (%) | Solver (s) |
|---|---:|---:|---:|---:|---:|---:|
| `0` | 95354 | 9.08551 | 173.86423 | 289.49178 | 31.99759 | 2024.549 |
| `0.001` | 95741 | 7.49253 | 7.23309 | 7.95859 | 1.07972 | 2423.184 |
| `0.0025` | 95742 | 8.22253 | 7.78633 | 7.96979 | 0.13526 | 2483.446 |
| `0.005` | 95347 | 9.88861 | 9.71150 | 9.24811 | 0.04818 | 2440.030 |
| `0.01` | 95743 | 13.65899 | 14.58413 | 13.63940 | 0.03613 | 2493.565 |
| `0.02` | 95745 | 20.52898 | 23.54191 | 22.49261 | 0.02529 | 2482.948 |
| `0.05` | 95736 | 34.49334 | 41.51166 | 40.69306 | 0.01545 | 2539.771 |

Late flux is the native radial peak after `0.19 s + radius/c_slow`, divided by
the full reference radial peak; it is not the reference-grid late-error metric
used in some earlier tables. All receivers contribute to each worst value.
The native sampling explains the small difference between the zero-gamma late
peak here and the earlier 31.83% reference-grid result.

The default-gamma job **95736** completes in **44:01**, with finite station
records and all 30,000 inertia / 60,000 storage true-residual checks passing.
The 1–48 Hz radial FAS errors at `0.05` are **23.70% velocity, 41.46% flux and
40.65% pressure**, versus **7.91%, 7.92% and 7.06%** at `0.005`. Its smaller
late tails accompany substantially attenuated slow-P pulses and worse full
waveforms. This demonstrates why selecting gamma solely by tail suppression
or algebraic convergence is insufficient.

Gamma `0.0025` improves velocity, flux and pressure L2/FAS agreement relative
to `0.005`, with a somewhat larger late flux peak. Its worst radial FAS errors
are **6.55669%, 4.67218% and 3.95159%**, respectively. Gamma `0.01` attenuates
the signal more strongly and gives worse full-waveform errors. Both intermediate
jobs complete with finite fields and all 30,000/60,000 inertia/storage residual
checks passing. Its separate 200-microsecond confirmation is recorded below.

Gamma `0.001` also completes with finite fields and all residual checks passing.
It lowers the worst full-waveform velocity/flux errors slightly further, with
radial FAS errors **5.83258%, 2.57925% and 2.36669%**, but leaves a **1.07972%**
late flux peak at the far diagonal receiver, versus **0.13526%** at `0.0025`.
The pressure L2 scores of these two values differ little. Separate candidate
checks at 200 microseconds are recorded below; agreement and residual tail
suppression are distinct selection criteria.

The final `0.02` case **95745** completes in **43:06**, with finite station
records and all 30,000/60,000 inertia/storage true-residual checks passing.
Its radial velocity/flux/pressure FAS errors are **15.03086% / 23.20660% /
22.09392%**. All seven gamma comparisons are complete.

For this 4 m, 16 Hz inviscid case, **`0.0025` is the recommended compromise
between pulse agreement and late-flux suppression**. Gamma `0.001` gives the
smallest tested worst-receiver full-waveform and FAS errors, but its late-flux
peak is about eight times larger. Thus `0.001` is preferable if full-pulse error
alone is the criterion; `0.0025` addresses the user's additional concern about
late oscillations while improving the full-pulse scores over `0.005`. This
choice is not a fitted optimum and does not improve every receiver or metric.
The global default remains `0.05`.

Both candidate values pass their own 200-microsecond confirmations below.
These results do not calibrate gamma or a larger timestep for heterogeneous,
viscous, topographic or coupled cases.

### Combining stabilization and timestep

The matched `gamma=0.0025` / 200-microsecond confirmation **95827** completes
the full 1.5 s record in **15:34**, with Solver time **832.386 s**. Against its
own 50-microsecond baseline **95742**, Solver time improves by **2.984x**.
Worst receiver temporal L2 changes are **0.28871% displacement, 0.43122%
velocity, 0.52510% Darcy flux and 0.52512% pressure**. All four pass the study's
chosen less-than-1% temporal-change criterion. All 7,500 inertia and 15,000
storage true-residual checks pass; every station record is finite.

The 200-microsecond worst author-reference velocity/flux errors are
**8.07684% / 7.45635%**; secondary SeisSol pressure error is **7.64182%**.
These remain separate from the temporal differences. A slight decrease in
reference error with a larger step can be cancellation of temporal and spatial
errors, rather than better spatial convergence. The source, receivers, mesh and
physical operators remain fixed; only timestep/output cadence and the previously
qualified private wave guard differ. The dissipation guard is unchanged.

The matched `gamma=0.001` / 200-microsecond confirmation **95829** also completes,
in **16:00**, with Solver time **845.612 s** and **2.866x** Solver speedup against
its own 50-microsecond baseline. Temporal displacement/velocity/flux/pressure
changes are **0.28300% / 0.42401% / 0.51286% / 0.51318%**. All finite-field and
7,500/15,000 inertia/storage true-residual checks pass. Its worst reference
velocity/flux/pressure errors are **7.33124% / 6.87639% / 7.62967%**. Both
candidate gammas therefore meet the same temporal-accuracy criterion.

These confirmations cover 1.5 s; they do not extend the separate `gamma=0.005`
ten-second observations to different stabilization values or change production
guards. The Figure 7 boundary experiment also retains its own `gamma=0.005`
material/source scope.

## Advisory timestep checks and stabilization default (2026-11-16)

Following the completed study, the requested production-policy update changes
`poroelastic_pressure_stabilization` from `0.05` to **`0.0025` when omitted**.
Explicit values, including zero, remain effective. This default applies to
active full-Biot stabilization with either mass matrix and flat/FEM/VT geometry;
it remains inactive for `brick_up` and fully sealed full-Biot material. The
homogeneous Wolf calibration is not a universal error bound for other materials,
interfaces, attenuation or nonlinear constitutive response.

Projected setup timestep limits are now **advisory**. A requested step above
the mesh/CPML estimate, VT wave estimate, full-Biot wave/stabilization estimate,
or consistent-mass wave/stabilization/stiffness-damping estimate produces a
solver-rank-zero warning in **stderr and the monitor file**. It does not abort
or reduce the requested timestep. Invalid input, unsupported configurations,
nonpositive mass matrices, PCG failures and other numerical errors retain their
existing error handling. No timestep-dependent physical operator is changed.

The hard-coded native/VT consistent wave coefficients remain `0.08`/`0.02`,
now warning thresholds. They multiply `h/c_fast`. For the Figure 5 material and
4 m spacing, 200 microseconds corresponds to approximately `0.1795*h/c_fast`.
The factor `200/501.466 ≈ 0.40` instead measures the fraction of the periodic
limiting step. Substituting `0.4` for `0.08` in the existing formula would give
about 446 microseconds, not 200. No universal new coefficient is inferred.

The preceding private surface-pressure experiment is unchanged by this update.
Its zero surface trace for both pressure and the stabilization auxiliary is
not promoted to a general production option. Original job inputs and reports
retain the defaults/guard policy of their pinned executables; the earlier
statements that production defaults were unchanged describe that study stage.

The isolated Horizon `warnings-v01` build uses the current source manifest and
default modules. Its executable SHA-256 is
`4221378a07b637ff8c24a229eba8c67f2b6a8c4239f1cac508c82afba709ded0`;
the comparison `joint-v02` executable is
`bc42033288070db6c9e2a5998e470c8f9f78953881ed2d7da832bce1fb94b4b4`.
Old-build comparisons explicitly set `gamma=0.0025`; new-build comparisons
omit it. Explicit `0.005`/`0.05`, zero-stabilization, sealed and reduced-`brick_up`
controls verify overrides and inactive paths. The original fixtures and source
records are preserved.

Jobs **96147, 96151 and 96157** completed 59 runs on two nodes/eight CUDA-aware
MPI ranks per allocation. Coverage includes flat/FEM/VT geometry, both mass
matrices, viscous/inviscid/sealed fluid, scalar acoustic pressure and velocity
potential, acoustic–elastic, porous–acoustic and porous–elastic contacts,
and representative attenuation/nonlinear-solid paths. These are execution and
before/after regression controls, not new accuracy validation for every
constitutive model. The distributed operators and both communication backends
are unchanged; this follow-up exercised CUDA-aware MPI only.

Five deliberately over-estimate controls retained their supplied step and
emitted each expected warning once in stderr and once in the monitor file,
with none in stdout:

| Control | Requested step (microseconds) | Relevant estimate (microseconds) | Record |
| --- | ---: | ---: | --- |
| Consistent, default gamma | 200 | 89.1283 | Complete small-domain pulse, approximately 0.6 s |
| Lumped, default gamma | 300 | 278.526 | Four steps |
| Consistent, explicit gamma `0.05` | 70 | 58.4562 | Four steps |
| VT consistent, gamma zero | 300 | 128.100 VT; 25.6200 consistent | Four steps |
| Elastic mesh estimate | 10000 | 9707.95 | Four steps |

All completed runs passed input/source/binary provenance, record-length,
finite-physical-field and applicable PCG true-residual checks. Ordinary station
pressure remains unavailable (`NaN`) in dry solid cells, with identical
before/after availability masks. The over-estimate tests establish warning
routing and continuation, not stability of arbitrary large timesteps.

The comparison screen retains a pre-existing weak-viscous-flow limitation
instead of declaring every `1e-6` result equivalent. For a short strongly
viscous control, the old/new aggregate Darcy-flux difference is **9.77%** at
`poroelastic_mass_tolerance=1e-6`, with reference peak `2.957e-10 m/s`;
solid velocity and pressure differ by **0.00619% / 0.000162%**. Repeating the
same old/new executables gives **8.78% / 8.30%** flux differences. A same-build
explicit-versus-omitted gamma control similarly differs by **10.48%**. At
`1e-10`, old/new flux difference falls to `4.84e-10` relative L2, solid
velocity to `5.12e-14`, and printed pressure agrees exactly. The tolerances
measure the global algebraic residual, not each weak physical component's
relative waveform error.

The complete 0.6 s viscous pulse initially differs by **0.5448%** in flux,
**0.0160%** in pressure and **0.001224%** in solid velocity at `1e-6`.
The matched `1e-10` controls in jobs **96162/96163** reduce old/new velocity,
pressure and flux differences to **`2.13e-11` / `4.68e-11` / `3.23e-11`**
relative L2. Repeating the same old/new executables at `1e-6` gives
**0.5770% / 0.5341%** flux differences, comparable to the old/new difference.
Together with the short-control repeats, these results support
existing tolerance/roundoff sensitivity rather than a changed physical
operator or gamma-parsing error. The default tolerance is unchanged by this
update; weak viscous fluid quantities still need their own tolerance study.

The final two jobs completed all four pulse controls, bringing the follow-up
to **63 completed runs**. All execution/provenance checks and the five warning
controls pass; the weak-viscous `1e-6` comparisons retain their recorded
relative-error flags. The tight-tolerance controls resolve those flags as
pre-existing sensitivity rather than evidence of a changed operator. Horizon
compilation/linking, documentation lint/type checks, both site builds and
rendered/parser audits pass. The default/reference table was also checked on
a narrow viewport. This update changes setup policy and the omitted gamma
value, not surface-pressure enforcement or physical evolution kernels.

## Suggested consistent-mass wave coefficient (2026-11-16)

Following the maintainer's request, the full-Biot consistent-mass Hex8 wave
estimate changes from `0.08*h/c_fast` to **`0.2*h/c_fast`** for flat and FEM
geometry. The VT coefficient remains `0.02`; lumped-mass estimates are unchanged.
For the Figure 5 material at 4 m spacing, the new wave suggestion is
**222.821 microseconds**. The prior 200-microsecond accuracy comparisons remain
the evidence for that tested step, not a measurement at exactly 222.821
microseconds or a universal accuracy guarantee.

The estimate remains advisory: larger supplied timesteps warn to stderr and
the monitor file and continue without alteration. Stabilization and physical
stiffness-damping estimates can still lower the combined suggestion. The
coefficient is internal, not a new input key or an automatic timestep choice.
Earlier tables retain the estimates from their pinned binaries.

The isolated `warnings-v02` executable has SHA-256
`bc3549f3f7ee18dd91bb660f1121d50c2dfd57f883c9177c2ef250cd6e4c333d`.
Job **96215** completed five controls on two Horizon nodes/eight CUDA-aware
MPI ranks in **1:50**, exit zero. The 200-microsecond small-domain pulse
reports **222.821 microseconds** and no timestep warning; a four-step
300-microsecond control reports warnings and preserves the requested step.
The FEM estimate increases by exactly the coefficient ratio **2.5**.
The explicit-`0.05` stabilization estimate remains **58.4562 microseconds**;
the zero-stabilization VT consistent estimate remains **25.6200 microseconds**.
Each expected warning appears once in stderr and monitor, absent from stdout.

All five controls pass input/source/binary provenance, finite-output and
applicable true-residual checks. Fixed-step comparisons against `warnings-v01`
have maximum aggregate relative L2 differences **1.46e-7** for the complete
pulse and **1.58e-6** for FEM, with stabilization/VT controls at roundoff.
No physical kernel or time-integrator change is included. The large-step
four-step controls check continuation only. Horizon compilation/linking,
documentation lint/type checks and both site builds/rendered audits pass.

## Opt-in weak full-Biot drained surface (2026-11-16)

The new `poroelastic_surface_pressure_enforcement=weak` option adds a numerical
surface drainage flux for `full_biot` with `poroelastic_surface_pressure_bc=drained`.
The default remains `natural`; reduced `brick_up` retains its existing per-step
compliance projection. The existing penalty factor supplies dimensionless beta,
with default `1`. It is not a measured hydraulic leakage coefficient.

For exposed permeable normals, `Y=beta*c_fast/M` times the squared projection
onto permeable coordinate directions. Assemble the positive, row-lumped surface
matrix `R` after pressure constraints. With `B=dt*R/2`, the pressure increment
satisfies `(C+B)*delta_p=F_step-2*B*p_old`. The coefficient `Y` is independent of
timestep. Both lumped and consistent storage include `B`; the stabilization
auxiliary retains `C`. Natural zero fluid-pressure traction and normal relative
fluid motion remain in the momentum equations. The numerical surface flux is
`q_n*=q_n+Y*p`; station fluid files still report the volume field `q`.
The added semi-discrete loss is `p^T R p >= 0`, not a proof of stability of the
whole time-discrete CPML/coupled system.

Flat/FEM surfaces use exposed mesh subfaces; native VT uses positive triangle
quadrature on reconstructed linear terrain contours inside retained P1
tetrahedra. Explicit seabed coverage excludes submerged quadrature points.
Internal acoustic/elastic contacts and sealed truncations receive no dry-surface
penalty. Curved terrain and shoreline quadrature require refinement. A partially
retained tetrahedron with no reconstructable vertex contour aborts explicitly.
Existing full-Biot, hydraulic-interface, nonlinear/attenuation and CPML
restrictions remain. The private strong experiment above is a different closure:
it removes both pressure and auxiliary surface unknowns and is not selected
by the new option.

`HERCULES_BIOT_SURFACE_DIAGNOSTICS=1` optionally writes integrated storage/RHS/
penalty-outflow balance and penalty energy loss each step. It adds a device
readback and MPI reduction and defaults to `0`. The integrated diagnostic uses a
different norm from the PCG acceptance test; their numerical thresholds are not
interchangeable. Its full format and limitations are in the
[Extended README](../../doc/Extended%20README.md#weak-full-biot-drained-surface).

The clean isolated Horizon build used current default modules. The current
physics implementation matches source snapshot `weak-v02`, executable SHA-256
`6d9d3aebfb70aa53f9df41ba220df833badea6bc87a589312a8197ded0255a4e`.
A first private `weak-v01` smoke exposed a missing term in the lumped scalar-MPI
overlap path. The corrected version evaluates that term inside each local nodal
update and audits the complete interior/boundary update. The failed private
lumped run is not validation evidence for the released path.

Completed implementation checks include extracted production face/VT quadrature
against exact planar areas, wet-surface exclusion, projected total flux,
positive row lumping, midpoint balance/energy, and second-order drainage-only
timestep convergence. These host checks do not validate GPU/MPI propagation.
Horizon checks exercise full/lumped mass, compact direct blocks, fused/unfused
kernels, both preconditioners/CG variants, element/row operators, cached guesses,
batched inertia/stabilization, owner-complete rows, standalone and coupled
flat/FEM/VT models, selected nonlinear dispatch, Rayleigh damping,
BKT/BKT2/BKT3/BKT3F attenuation, and both MPI backends. These are execution, algebraic-balance and before/after checks, not
independent reference validation for every physical combination.

The attenuation fixtures use synthetic Qp=Qs=50 and compare before/after
default behavior as well as weak execution. Their first input attempt omitted
required Q controls and aborted at parsing; corrected job 96346 completes all
12 cases with weak balance residuals below `3.0e-9`. These are regression
controls, not calibrated attenuation benchmarks.

The corrected lumped overlap/sequential pulse balance residuals are about
`1.7e-17` relative to the diagnostic scale. Combined one-node CUDA-aware versus
two-node host comparisons give maximum aggregate differences of `2.99e-6`
solid velocity, `5.61e-7` pressure and `1.23e-7` Darcy flux in relative L2.
The porous/elastic integrated balance decreases from `2.16e-6` at PCG tolerance
`1e-6` to `1.27e-10` at `1e-10`. A tiny viscous fluid signal remains
sensitive to the loose solve tolerance: a before/after default-path difference
of 9.2% of its `2.88e-10 m/s` peak signal falls to `6.84e-10` relative L2 at
`1e-10`. This resolves that implementation comparison, not the accuracy of
weak viscous signals at the default tolerance. Zero beta recovers the natural
condition within solve/reduction roundoff. Invalid enforcement names, reduced
formulation, undrained selection and nonfinite penalty factors are rejected.

### Matched Figure 7 results

All 13 large comparisons finish with complete finite receiver records and
true PCG residuals within `1e-6`. The uniform `4096 m` cube, original 5 Hz
source at 500 m depth, original receiver coordinates, gamma `0.0025`,
consistent mass/storage, and 256 m CPML are fixed. The CPML defaults are
reflection target `0.01`, order `2`, kappa maximum `5` and alpha floor fraction
`0.1`; both meshes log alpha maximum `3.14 rad/s` and automatic reference
speed `4246.852 m/s`. Partition weight is explicitly `1`. Thickness and
partition weight are case overrides. Default enforcement remains `natural`.

The table gives relative L2 differences (%) against archived **SeisSol numerical**
pressure and Darcy-flux data. Flux norms include all three components. Surface
peak is at the diagnostic receiver with the same horizontal source offset
`(500, 400) m`; it is not a maximum over the entire surface.

| h (m) | Treatment | Shallow p, 0.5 m (%) | Shallow q (%) | Deep p, 500 m (%) | Surface pressure peak (Pa) |
|---:|---|---:|---:|---:|---:|
| 8 | natural | 39.12644 | 0.83127 | 0.72336 | 0.3535117 |
| 8 | weak | 6.64966 | 0.81339 | 0.72040 | 0.07307303 |
| 8 | strong | 1.57317 | 0.36422 | 0.66497 | 0 |
| 4 | natural | 11.99759 | 0.34246 | 0.29916 | 0.1244364 |
| 4 | weak | 2.91413 | 0.19521 | 0.29820 | 0.03397712 |
| 4 | strong | 0.48795 | 0.29765 | 0.30709 | 0 |

The weak rows use beta `1` and `200 microseconds`. Pressure and flux improve
under refinement in this configuration. Weak versus natural solid-velocity
histories differ by at most `0.00260%` across the four 4 m receivers. The
physical pressure update reuses the existing solve; it adds no extra global
solve. Optional surface balance diagnostics are off in these large runs.

The `strong` rows use the prior private prototype, which constrains both
pressure and stabilization-auxiliary traces. It was not a released input option at the time of that study
and is not the infinite-penalty version of this unchanged auxiliary space.
The subsequent public strong option is documented below.
Its executable SHA-256 is
`79614b76bd66b9234dc4f80cbd976dfa37bc4ffb3f88d21f425d28a233cb54ce`.

### Penalty and timestep sensitivity

At 8 m and 200 microseconds, holding every other input fixed:

| beta | Shallow p L2 (%) | Shallow q L2 (%) | Surface pressure peak (Pa) |
|---:|---:|---:|---:|
| 0.25 | 19.51361 | 0.65586 | 0.1855908 |
| 1 | 6.64966 | 0.81339 | 0.07307303 |
| 4 | 1.01869 | 0.92411 | 0.02104434 |
| 16 | 1.31400 | 0.95965 | 0.005458444 |

A penalty that minimizes pressure error need not minimize fluid-flow error.
These are two-second tests; the larger penalties have not received the
ten-second or 4 m checks. No universal optimal penalty is selected.

| Weak timestep compared with 100 microseconds | Shallow p change (%) | Shallow q change (%) | Surface p change (%) |
|---:|---:|---:|---:|
| 200 microseconds | 0.007009 | 0.011296 | 0.187999 |
| 400 microseconds | 0.036316 | 0.053254 | 0.541040 |

The 100/200/400-microsecond comparison uses beta `1` on the 8 m mesh.
All time-domain scores use raw fields with their recorded staggering; spectra
use the same 10% Tukey taper and `1–12.5 Hz` band. These checks do not establish
a universal critical timestep. The existing timestep recommendations remain
advisory; no new timestep rejection is introduced.

### Ten-second controls

The 8 m natural and beta-1 weak cases both complete 50,000 steps at
200 microseconds. Below are peak magnitudes during 8–10 s divided by the
corresponding first-two-second peak. Flux magnitudes combine three components.

| Treatment | Receiver depth (m) | Late p / pulse p | Late q / pulse q |
|---|---:|---:|---:|
| weak | 0.5 | 5.001486e-05 | 4.113223e-05 |
| weak | 500 | 2.509561e-05 | 2.404826e-05 |
| natural | 0.5 | 5.051708e-05 | 4.076859e-05 |
| natural | 500 | 2.507721e-05 | 2.402806e-05 |

All receiver samples are finite. Windowed late amplitudes and histories are
retained in the local report. Shallow pressure and flux show a small late
rebound, so decay is not monotonic. The first two seconds agree with the
short controls to within `0.00049%` relative L2 across the four receivers
and all recorded fields. These station records do not prove whole-domain
energy stability, arbitrarily long-time stability, or coastal-interface
stability. The CPML/interface restrictions remain unchanged.

### Reference and validation limits

The primary solid-motion reference remains the supplied Gar6More3D displacement
with fixed author scale `1e10/(2*pi^2*5^2)` and time convention `+0.001 s`.
Velocity is differentiated with the same centered NumPy gradient used in the
earlier comparisons; the authors’ derivative helper is absent. No gain, delay
or corrective filtering is fitted. The original reference feature near 1.5 s
is retained, and before-1.3 s scores are reported separately.

| Weak beta-1 mesh | Shallow displacement L2 (%) | Shallow derived velocity L2 (%) | Deep derived velocity L2 (%) |
|---:|---:|---:|---:|
| 8 m | 1.73870 | 2.21702 | 4.43663 |
| 4 m | 0.73734 | 0.90050 | 4.13353 |

The archived SeisSol velocity itself differs from differentiated Gar6More3D
by `0.23890%` shallow and `3.84574%` deep. The deep gap persists
before 1.3 s; it is not explained away by excluding the later reference feature.
It remains a reference/processing limitation, not proof of Hercules accuracy.
Pressure and Darcy flux use a secondary numerical reference because no
independent analytical versions are supplied for Figure 7.

The quantitative surface study is flat, homogeneous and inviscid. FEM/VT,
variable material, mixed acoustic/elastic/porous, selected nonlinear, viscous
and attenuation checks exercise implementation and transport; they do not
establish independent surface/scattering accuracy for all those combinations.
Curved-terrain/shoreline convergence and a complete CPML energy budget remain
outside the demonstrated scope.

Large jobs: `96315` (4 m closures), `96316` (8 m closures), `96317`
(penalty/timestep), `96343`/`96344` (ten-second weak/natural). All used
`cuda-615-SC` and current default Horizon modules. No existing reservation
was changed. Small checks and their resolved historical failures are recorded
above; uninvestigated audit failures: zero.

The large numerical comparisons use `weak-v02`, SHA-256
`6d9d3aebfb70aa53f9df41ba220df833badea6bc87a589312a8197ded0255a4e`.
The final `weak-v03` change only initializes the optional diagnostic file
during setup, including when `number_output_stations = 0`; no numerical
kernel, equation or quadrature changes. Its SHA-256 is
`547254b54b3502670ad058f516312c06088021d54b0a663f19e6c89882f66a8d`.
Job `96422` checks diagnostics on/off against v02, custom/default directories
without station output, and lumped audit output. All five runs pass;
waveforms and integrated outflow/energy agree within `1e-4` relative L2.
The optional audit adds no file or readback when disabled.

Investigation inputs, source snapshots, regressions, outputs and helper scripts
remain local-only and are
unavailable in a clean checkout. Jobs and source/binary/input hashes are retained
there; no production binary or pre-existing Horizon job was replaced.

## Strong drained surface with FEM and VT (2026-11-16)

`poroelastic_surface_pressure_enforcement = strong` was introduced as an opt-in
for active `full_biot` with `poroelastic_surface_pressure_bc = drained`;
`natural` was the default during this study. The subsequent geometry-dependent
default change is recorded below. The implementation restricts pressure and stabilization-auxiliary
trial/test spaces to zero on the positive support of exposed permeable dry
traces, including positive hanging-constraint parents. Normal relative fluid
motion remains free. Both consistent and lumped scalar updates use this
restriction; it is not a post-solve pressure clamp. Mechanical inertia, Darcy
SDIRK stages and constitutive/interface equations retain their existing form.
There is no additional per-step MPI exchange. Weak penalty/balance controls are
inactive for strong; setup reports the fixed-node count and dry permeable area.

Conforming flat/FEM faces receive ordinary nodal Dirichlet conditions. Native
VT uses an **experimental conservative support restriction**: every pressure node supporting a
reconstructed cut trace is fixed. This also suppresses pressure variation
throughout affected cut tetrahedra. It is **not a minimal cut-trace or fitted
boundary method**. Nonmatching faces and shorelines can extend fixed support
outside the exposed dry patch. The implementation is usable as an opt-in
boundary approximation; zero trace alone does not establish near-surface
accuracy. Existing groundwater, restart, VT/CPML and interface restrictions
remain. Strong-boundary reactions and a whole-system fluid-content/energy audit
are not currently output.

### Build, regression and transport evidence

The isolated Horizon clean double-precision CUDA/MPI build used default modules (`ml reset`),
retained GSL/PROJ dependencies and binary SHA-256
`550682fb135338750c9c3b159e80361ed1d8fbcfbecbff5107ffbdfd54801b0b`. The GPU structure changed, so every object was rebuilt.
There are **149 accepted run records** and **3
isolated expected-rejection checks**. Inputs, sources, executable hashes,
station counts, finite fields and reported true PCG residuals are audited.
Baseline natural controls use pre-change `weak-v03`; weak before/after and
explicit-gamma private-strong controls are retained separately.

The matrix covers elastic/acoustic/reduced controls, full-Biot standalone and
acoustic/elastic contacts, both acoustic formulations, field-specific/shared
geometric paths where supported, FEM/VT, heterogeneous/anisotropic/sealed
material, both mass matrices, selected nonlinear dispatch, viscous Darcy
stages, Rayleigh and BKT-family attenuation. Compact/full vectors, fused and
unfused updates, batch solves, one-round communication where eligible,
single-reduction/standard CG and both mass operators/preconditioners are
exercised. These checks establish implementation consistency, not independent
scattering or nonlinear-yield accuracy for every material/interface.

Four-rank CUDA-aware and eight-rank host transport comparisons agree within
`2.71e-6` solid-velocity, `7.17e-7` pressure and `9.49e-8` flux relative L2
for the ordinary nonzero-signal controls. A tiny viscous fluid signal is
sensitive to `1e-6` PCG tolerance: its cross-layout flux difference is `10.207%`.
The matched `1e-10` control reduces that difference below `1.94e-11` relative
L2. Tight before/after viscous controls similarly resolve their loose-tolerance
differences (all field norms below `3.3e-10`). This does not qualify tiny-fluid
accuracy at the default tolerance, which is unchanged.

The wholly impermeable-normal control (`kz=0`, `kx,ky>0`) reports zero fixed
surface nodes and agrees with natural enforcement within `1.78e-6` relative
L2. Strong ignores the weak penalty factor, including zero. Extracted
production quadrature/algebra tests cover both tetrahedral partitions,
positive trace support, wet exclusion, nonmatching fragments, hanging closure
and projected storage; they supplement, rather than replace, the GPU runs.

### Independent planar reflection and refinement

The reference is an independently implemented one-dimensional inviscid Biot
modal solution, not a Hercules fine-mesh output or a digitized Wolf curve.
It solves the generalized two-field compressional eigenproblem, launches both
P modes from a mechanical force sheet, and imposes zero total normal traction
and zero pore pressure at the drained surface. Modal boundary residuals and
incident/reflected energy balance are checked. No Hercules element matrices,
fitted delay, gain or filtering enter the reference.

The domain is `2048 x 2048 x 1024 m`; the planar surface is at solver
`Z=64 m` for FEM and off-grid `Z=67.2 m` for VT. The source plane is at
`Z=192 m`, with a 12 Hz Ricker central frequency and time center `0.08 s`.
Uniform retained physical spacing is 8/4/2 m, timestep is fixed at
25 microseconds, gamma is `0.0025`, consistent PCG tolerance is `1e-8`, and
duration is `0.32 s`. Receiver depths are 32/64/192 m below the actual plane;
VT additionally has a 1 m receiver. Velocities/fluxes use their recorded
half-step timing. Inherited source-file sample spacings are 200/100/50
microseconds, interpolated onto the common step; source L2 differences versus
the finest history are only `0.00818%` and `0.00379%`, but are not zero.

The following full-record relative L2 errors are the maximum over the common
32/64/192 m receivers, in percent:

| Method | h (m) | Enforcement | Solid velocity (%) | Pressure (%) | Darcy flux (%) |
|---|---:|---|---:|---:|---:|
| FEM | 8 | natural | 1.17248 | 5.48335 | 4.87748 |
| FEM | 8 | strong | 1.17286 | 4.69897 | 3.59153 |
| FEM | 4 | natural | 0.29912 | 1.12917 | 0.81506 |
| FEM | 4 | strong | 0.29908 | 0.64555 | 0.60627 |
| FEM | 2 | natural | 0.09905 | 0.30174 | 0.18232 |
| FEM | 2 | strong | 0.09850 | 0.22082 | 0.15892 |
| VT | 8 | natural | 0.81274 | 5.18351 | 4.17414 |
| VT | 8 | strong | 1.11271 | 32.06482 | 21.56501 |
| VT | 4 | natural | 0.81501 | 1.63533 | 0.89086 |
| VT | 4 | strong | 0.83905 | 5.12224 | 4.41522 |
| VT | 2 | natural | 0.34156 | 0.56733 | 0.46539 |
| VT | 2 | strong | 0.38186 | 5.34811 | 4.13827 |

The deliberately shallow VT receiver is reported separately, without omitting
its errors from the original target audit:

| VT h (m) | Enforcement | Solid velocity at 1 m (%) | Pressure at 1 m (%) | Darcy flux at 1 m (%) |
|---:|---|---:|---:|---:|
| 8 | natural | 0.41666 | 171.83131 | 2.72278 |
| 8 | strong | 0.42504 | 100.00000 | 44.09762 |
| 4 | natural | 0.40169 | 13.54637 | 0.81299 |
| 4 | strong | 0.39541 | 77.78421 | 1.32712 |
| 2 | natural | 0.16049 | 4.36385 | 0.29919 |
| 2 | strong | 0.15606 | 75.07946 | 1.41722 |

Original fixture targets are 5% solid velocity and 10% pressure/flux at each
receiver. Errors above those thresholds remain accuracy gaps; successful
execution is not a target pass. For this fixed off-grid surface the 4 m and
2 m meshes both have their next node plane at `Z=68 m`, 0.8 m below the true
surface. The conservative VT zero-support band can therefore retain a large
near-surface relative error even when the deeper wavefield improves. This is
a reason to inspect cut placement and refinement, not to use zero pressure as
an accuracy certificate. Only the pressure/auxiliary space is restricted;
the mechanical surface is not translated. A secondary calculation at the
known lower node plane does not account for all the waveform discrepancy and
never replaces the true physical-surface reference above.

The continuum earliest lateral-edge influence is about `0.285 s`
(`1024/3590.329`); analysis also retains scores before `0.27 s`. The full
`0.32 s` window can include edge effects. This continuum estimate is not a
strict numerical causality bound or a universal timestep limit.

### Curved-terrain refinement and trace checks

The fixed 512 m cube uses the same 6 Hz source and five physical receivers,
gamma `0.0025`, PCG `1e-8`, 25-microsecond timestep and `0.4 s` duration.
The table reports maximum differences against the next finer mesh, in percent.
These are finite-domain **self-convergence differences, not analytical errors**.

| Method / enforcement | h comparison (m) | Solid velocity (%) | Pressure (%) | Darcy flux (%) |
|---|---|---:|---:|---:|
| FEM / natural | 4 vs 2 | 4.66996 | 7.55432 | 5.62238 |
| FEM / natural | 2 vs 1 | 2.27322 | 4.57959 | 2.70650 |
| FEM / strong | 4 vs 2 | 4.66717 | 9.15724 | 5.35867 |
| FEM / strong | 2 vs 1 | 2.27156 | 4.23648 | 2.60955 |
| VT / natural | 4 vs 2 | 0.40077 | 2.55815 | 1.05223 |
| VT / natural | 2 vs 1 | 0.20831 | 1.64209 | 0.37063 |
| VT / strong | 4 vs 2 | 0.43046 | 9.99677 | 2.44023 |
| VT / strong | 2 vs 1 | 0.22164 | 3.54683 | 1.13550 |

Exact supported FEM/VT planar trace stations and four retained curved-VT
surface stations record zero pressure with nonzero relative fluid flow.
A physical terrain-depth-zero point need not lie on a FEM stair-step face;
the corresponding pressure is not required to vanish. One curved VT point
lies in a discarded tetrahedron and is correctly rejected. Its raw failed
fixture is retained; the successful retained-point control does not validate
the discarded location. Curved geometry and shoreline support remain
approximations, and the hill data do not establish universal convergence order.

### Matched Wolf Figure 7 production control

Job `96528` repeats the flat inviscid 8 m case with the public strong option:
4096 m cube including 256 m CPML, 5 Hz source, 2 s duration,
200-microsecond timestep and gamma `0.0025`. Other settings match the earlier
private control `96316`. Maximum aggregate difference over recorded fields is
`0.00008733%` relative L2.

| Depth (m) | Gar6More3D displacement (%) | Derived velocity (%) | SeisSol pressure (%) | SeisSol flux (%) |
|---:|---:|---:|---:|---:|
| 0.5 | 1.68818 | 2.15074 | 1.57317 | 0.36422 |
| 500 | 3.77434 | 4.42265 | 0.66497 | 0.67225 |

The diagnostic surface pressure is exactly zero and its normal volume-flux
peak is `3.27677417e-06 m/s`. Pressure/flux
references are archived **SeisSol numerical data**, not analytical solutions.
Gar6More3D supplies displacement; velocity uses the established centered
differentiation, fixed scale and author time convention. No fitted correction
is applied. The original reference feature near 1.5 s and the earlier deep
velocity-reference discrepancy remain. This production run covers 2 s; the
older ten-second flat prototype evidence does not qualify all topographic
or coupled configurations. Single-precision accuracy is not established by
these double-precision runs.

### Provenance and limits

All new jobs use `cuda-615-SC`; no existing reservation was changed.
The coarse/fine/hill batches `96525/96526/96527` were stopped only after
their FEM cases completed; their VT cases run in separately submitted parallel
continuations `96584`–`96592`. Independent coarse VT hill controls can start
earlier on otherwise idle nodes. Partial attempts are not counted as passes.
The timed-out surface batch `96539` is retained, with its unfinished cases
completed by `96582`. Initial confounded gamma and parser-guard fixtures and
the discarded-VT-station failure remain in the raw audit with their corrective
controls. The extra shared/local-matching controls initially stopped before
simulation in `96606` because their inputs had not been staged; `96608`
repeats the same immutable fixtures after staging. The failed launch remains
recorded separately. Shared-constraint FEM/VT controls complete; all three
local-matching terrain attempts correctly hit the unchanged topography guard.
Flat local-matching controls for both acoustic formulations complete in `96617`.
Those submerged flat fixtures report zero dry-surface fixed nodes; they check
compatibility and inactive enforcement, not active local-matching drainage.
This feature does not make local-matching compatible with topography.
Unresolved execution/paired-comparison/MPI audit flags: zero.
Scientific accuracy limits are reported separately above, not waived by that
execution audit.

Immutable inputs, source snapshots, reference arrays, job states, station
records, analyses, figures and the PDF report remain local-only. Those investigation
helpers
are unavailable in a clean checkout. No production executable or unrelated
Horizon job was replaced. No general terrain/shoreline accuracy guarantee,
whole-system drainage balance, full CPML energy stability or universal
long-time stability follows from these tests.

## Geometry-dependent drained-surface defaults (2026-11-16)

Omitting `poroelastic_surface_pressure_enforcement` now selects:

| Active formulation / hydraulic surface / geometry | Enforcement default |
|---|---|
| `full_biot`, `drained`, topography off | `strong` |
| `full_biot`, `drained`, FEM topography | `strong` |
| `full_biot`, `drained`, VT topography | `natural` |
| Undrained surface, BrickUP or inactive poroelasticity | `natural` |

The hydraulic condition itself still defaults to `undrained`. This change does
not automatically drain a model. Explicit `natural`, `weak` or `strong` remains
authoritative, subject to the existing requirement that weak/strong need active
full Biot and a drained surface. There is no new `auto` input value. An irrelevant
`computation_method=vt` with topography disabled does not select the VT default.
To preserve the previous treatment of a flat/FEM drained input that omitted the
key, explicitly set `poroelastic_surface_pressure_enforcement=natural`.

The parser carries an internal unresolved marker through the existing integer
broadcast. Every rank resolves it after topography initialization using the
actual FEM/VT method. Startup reports the effective enforcement, geometry,
hydraulic BC and whether selection came from the default or explicit input.
This changes selection only: pressure/storage operators, strong/weak surface
geometry, stabilization, inertia, Darcy integration and MPI field schedules are
unchanged. In particular, VT strong retains the experimental approximation and
accuracy limitations documented above; it remains available explicitly.

An isolated clean double-precision Horizon build used `ml reset`, retained
GSL/PROJ paths and executable SHA-256 `afca6eb5ad04b5f7551b8f2ed4d5dbff303d078a0d0517f43165176067a6856e`. The comparison executable is the
preceding strong-surface build `550682fb135338750c9c3b159e80361ed1d8fbcfbecbff5107ffbdfd54801b0b`
(source parent `f7e3120`). No production executable was replaced.

Job `96716` runs twenty before/selected pairs on one node and four CUDA-aware
MPI ranks; job `96717` repeats four selected cases on two nodes and eight host
MPI ranks. Both use `cuda-615-SC`; no existing reservation was changed. All
**44 run records**, twenty paired comparisons and four transport comparisons
pass the finite-output, immutable-input/source/binary, station-count, selected-mode
and reported true-residual audits. Maximum relative L2 differences are:

| Comparison | Solid velocity | Pressure | Darcy flux |
|---|---:|---:|---:|
| Selected versus prior explicit choice | 3.3028765e-06 | 8.8453888e-07 | 1.5883182e-07 |
| Four-rank CUDA-aware versus eight-rank host | 2.1686048e-06 | 5.6575905e-07 | 1.5112098e-07 |

These are relative fractions, not percentages. The short fixtures exercise
flat/FEM/VT, lumped/consistent mass, both acoustic scalar formulations,
open/no-flow interfaces, elastic and selected nonlinear dispatch, viscous drag,
BrickUP drainage, inactive full-Biot selection, the unchanged undrained BC
default, case-insensitive geometry and explicit natural/weak/strong overrides.
The viscous pair retains a tight `1e-10` solver tolerance. They establish default
selection and regression consistency, not new continuum accuracy, long-time
stability, nonlinear-yield accuracy or a repeat of the full Wolf benchmarks.
The preceding strong/weak/natural scientific comparisons retain their original
binary, mesh, timestep and boundary-selection scope.

An extracted-production C++ parser/resolver check passes 240 combinations of
active/inactive poroelasticity, formulation, geometry, boundary condition,
explicit/omitted/invalid choices and reporting rank. It rejects premature
access to uninitialized geometry and preserves invalid-combination guards.
Both documentation build targets and their parser/rendered audits pass; desktop
and narrow-screen parameter links/defaults are checked in the browser.

Immutable fixtures, raw outputs, source manifests, scripts and detailed audit
JSON remain local-only;
these helpers are unavailable in a clean checkout. No scientific reference
output or existing benchmark input was changed to obtain agreement.

## VT natural-versus-weak accuracy study (2026-11-17)

**Keep natural as the VT default.** Weak beta `1` improves the unresolved
8 m shallow pressure trace but worsens deeper pressure and Darcy flux. At
4 m and 2 m, natural gives lower shallow and deeper pressure errors. A smaller
beta can improve one measurement while worsening another, and the ranking
changes with the terrain's position inside the cell. No default or solver
change is recommended from this study.

All sixteen new immutable cases passed completion, input/binary integrity,
finite output, station-count and true-PCG-residual audits. The independent
reference uses continuum two-field Biot modes and zero-traction/zero-pressure
reflection conditions; modal energy and boundary equations are checked. No
Hercules matrices, fitted gain, time shift or filtering enter this reference.
The numerical and reference fields are compared at their actual integer
pressure and half-step velocity times.

The main planar controls use an inviscid homogeneous material, consistent
mass, gamma `0.0025`, PCG tolerance `1e-8`, a `2048 x 2048 x 1024 m` domain,
surface `Z=67.2 m`, source sheet `Z=192 m`, 12 Hz central frequency and
`0.32 s` duration. The timestep is 25 microseconds at 8/4/2 m spacing. Source
histories are unchanged within each pair; inherited source sampling across
refinement differs by at most 0.0082% source L2 from the finest history.
The table gives full-record relative L2 errors; deeper maxima include the
32/64/192 m receivers, separately from the sensitive 1 m receiver. Absolute
errors and pre-0.27 s scores are also retained: lateral sheet-edge effects can
enter after roughly 0.285 s. The 1 m reference pressure is small but nonzero;
its large relative errors are not discarded.

| h (m) | Enforcement | p at 1 m | q at 1 m | Max deeper p | Max deeper q |
| ---: | --- | ---: | ---: | ---: | ---: |
| 8 | natural | 171.8234% | 2.7224% | 5.1841% | 4.1755% |
| 8 | weak, beta 1 | 128.2195% | 21.8995% | 17.8788% | 12.6701% |
| 4 | natural | 13.5464% | 0.8129% | 1.6354% | 0.8910% |
| 4 | weak, beta 1 | 26.1702% | 0.9364% | 1.9554% | 2.2601% |
| 2 | natural | 4.3638% | 0.2992% | 0.5673% | 0.4655% |
| 2 | weak, beta 1 | 16.2471% | 0.4061% | 1.2836% | 1.0425% |

### Penalty, cut placement and timestep

| h=4 m enforcement | p at 1 m | Max deeper p | Max deeper q |
| --- | ---: | ---: | ---: |
| natural | 13.5464% | 1.6354% | 0.8910% |
| weak beta 0.25 | 4.3413% | 1.5102% | 1.2562% |
| weak beta 1 | 26.1702% | 1.9554% | 2.2601% |
| weak beta 4 | 62.1226% | 3.9681% | 3.7698% |
| weak beta 16 | 75.2977% | 4.9724% | 4.3609% |

Beta `0.25` improves shallow pressure in this one cut placement, while deeper
flux is less accurate than natural. Increasing beta is not reliably beneficial.
At surface `Z=66 m`, beta `1` reduces shallow pressure error from 74.5825% to
64.9219% but increases maximum deeper pressure error from 0.5161% to 6.2485%.
At `Z=67.8 m`, it improves deeper pressure/flux from 0.6938%/1.1136% to
0.5792%/0.8090%. These positions prevent choosing a universal beta from a
single surface cut.

The native `min_enclosed_volume_perc=0.1` remains fixed. Extracted-production
geometry checks give represented surface fractions of 100%, 64% and 90.25%
at within-cell cut coordinates 0.5, 0.8 and 0.95. With cutoff zero, the host
geometry check gives the complete plane, but no cutoff-zero GPU propagation
run is claimed. Weak acts on retained tetrahedra; discarded material and
native volume quadrature already affect natural runs. Surface area alone
does not predict waveform error.

The h=4 m beta-1 half-step control changes matching-time fields by at most
0.002440% relative L2. Fresh h=8 m natural and h=4 m beta-zero runs
reproduce their archived natural controls within 0.000032% aggregate
relative L2. This verifies reuse and beta-zero recovery at the tested settings.

### Curved terrain

The 512 m-cube hill uses a 6 Hz source, five receivers, 0.4 s duration and
25-microsecond timestep. These are differences against the next finer mesh,
not analytical errors, and do not determine which closure is closer to the
exact curved-terrain solution.

| Enforcement | Mesh pair (m) | Max v difference | Max p difference | Max q difference |
| --- | --- | ---: | ---: | ---: |
| natural | 4 vs 2 | 0.4008% | 2.5581% | 1.0522% |
| natural | 2 vs 1 | 0.2083% | 1.6421% | 0.3706% |
| weak | 4 vs 2 | 0.4035% | 3.6932% | 1.1696% |
| weak | 2 vs 1 | 0.2091% | 1.1574% | 0.4086% |

Both closures show smaller differences with refinement. Outer-boundary
approximations remain common to these hill comparisons. This study does not
establish viscous surface accuracy, general coastal/CPML long-time stability,
a universal timestep limit, or a global optimum penalty. Strong is a separate
experimental VT support restriction and is not tested in this comparison.

Source commit `754934f`; executable SHA-256
`afca6eb5ad04b5f7551b8f2ed4d5dbff303d078a0d0517f43165176067a6856e`.
Jobs `96736`–`96742` completed on `cuda-615-SC`; no existing reservation was
changed. `IO_testing` did not authorize the account at submission. Inputs,
raw records, manifests, derived metrics, waveform/refinement/penalty/FAS figures
and the reviewed PDF remain local-only, unavailable in a clean checkout.
Existing references and benchmark inputs were not modified to obtain agreement.

## Fixed buried groundwater for full Biot (2026-11-17)

This extension supersedes the earlier full-Biot exclusion of buried/spatial
groundwater grids. It uses `surface_groundwatertable=yes` and
`groundwaterlevel_file`, with the established ground-relative, bilinearly
interpolated depth grid. It adds no new input parameter or default change.
Full Biot supports flat/FEM/VT **whole-element saturation** and a hydraulic
condition at the resulting mesh-face porous/dry-cap boundary. It does not add
subcell saturated volumes, unsaturated/capillary flow, a moving water table,
or automatic conversion of saturated material properties to dry values.

The dry cap keeps its supplied single-phase solid properties and shares solid
motion with the porous region. Owner-confirmed above-groundwater metadata
separates that interface from other elastic/material/depth-limit contacts.
The latter remain sealed. The groundwater contact follows
`poroelastic_surface_pressure_bc`: `undrained` (default) constrains normal
relative flow; `drained` uses `natural`, `weak` or `strong`. Omission still
selects strong for drained flat/FEM and natural for VT. Native VT terrain-cell
contacts integrate retained P1 face triangles clipped to terrain and neighbor
subfaces. All new work is setup-time classification/quadrature; existing
inertia, storage, Darcy, pressure, GPU/MPI and output paths are retained.

Groundwater depths must be finite and nonnegative for full Biot. The grid covers
the complete horizontal footprint, with positive spacing that divides both
horizontal lengths, then north-major/east-minor values. Zero depth reproduces
surface saturation, including retained terrain cells whose cube centers lie
above the continuous surface. Positive buried depths use the actual cube center;
VT does not substitute the centroid of its retained solid volume. Wet explicit
seabed columns retain water-surface saturation. The porous depth limit remains
relative to ground/bed. Other full-Biot/CPML restrictions are unchanged; reduced
BrickUP porous CPML still requires surface saturation.

### Independent buried-interface comparison

These are additional interface tests, not reproductions of Wolf's free-surface
Figure 7. The model has a dry elastic cap above a porous half-space. The
continuum reference solves the two-field Biot generalized eigenproblem and
matches **solid velocity and total traction** across the interface, with either
zero pore pressure or zero normal relative flow. It is not an acoustic open-pore
condition: solid velocity, rather than solid velocity plus Darcy flux, is
continuous into the dry solid. Boundary residuals and modal energy balance
are checked independently of Hercules's matrices.

The qualified reference domain is a **2048 m cube** with uniform flat meshes,
water table at solver
`Z=512 m`, and a vertical mechanical force sheet at `Z=640 m`, 128 m below the
table. The earlier `2048 x 2048 x 1024 m` models are retained as bottom-distance
controls: moving the source deeper made their old bottom capable of returning
within the comparison window. The doubled-depth cases also extend the porous
depth limit to 2048 m, so an artificial material cutoff does not remain at the
old bottom. Source, receivers, spacing, timestep and physical material are
unchanged. No CPML is used in these planar controls; the existing mechanical
boundary treatment is retained. Lateral sheet-edge effects can enter the full
record after about 0.285 s. Primary scores use `t<0.27 s`; full 0.32 s scores
and matched boundary-distance differences are retained in the local analysis.
The window does not contain every late slow-mode return at the deepest probe.

The source has 12 Hz central frequency, Ricker time shift 0.08 s and force per
area `2*rho*Vp_undrained*0.001`. Material values are drained Vp 3296.268286 m/s,
Vs 1961.161351 m/s, mixture rho 1664 kg/m3, pore-fluid rho 1040 kg/m3,
porosity 0.4, alpha 0.522666667, M 6.019261637 GPa, tortuosity 2, viscosity
zero and positive permeability 6e-9 m2. The reference P speeds are 3590.329282
and 1055.976329 m/s. The dry cap uses the same supplied density/Vp/Vs as a
single-phase elastic material for this synthetic comparison. Consistent mass,
gamma 0.0025, PCG tolerance 1e-8 and a 25-microsecond timestep are fixed.
No fitted gain, time shift or output filtering is used. Pressure is evaluated
at integer time and velocities at `t-dt/2`.

Receivers are at depths -128 m (dry cap), 0.5, 32, 64 and 192 m relative to the
water table. The table reports relative L2 errors before 0.27 s. Bulk maxima
cover only the three 32/64/192 m receivers; the shallow pressure and transmitted
dry-cap velocity are shown separately. Targets of 5% velocity and 10%
pressure/flux apply to the bulk receivers and do not conceal shallow errors.

| Geometry / BC / enforcement | h (m) | Bulk v | Bulk p | Bulk q | p at 0.5 m | Dry-cap v |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| flat natural | 4 | 0.12665% | 0.94239% | 1.31041% | 21.44203% | 0.12534% |
| flat strong | 4 | 0.12736% | 0.62794% | 0.77124% | 4.05221% | 0.12439% |
| fem strong | 8 | 0.46673% | 4.37707% | 4.05159% | 10.51185% | 0.49937% |
| flat natural | 8 | 0.46198% | 5.96574% | 5.46499% | 100.35861% | 0.50242% |
| flat strong | 8 | 0.46673% | 4.37708% | 4.05159% | 10.51185% | 0.49347% |
| flat undrained | 8 | 0.45757% | 6.03397% | 5.87193% | 7.91692% | 0.44152% |
| flat weak | 8 | 0.46349% | 3.47628% | 4.00271% | 17.79691% | 0.49794% |
| vt natural | 8 | 0.46198% | 5.96574% | 5.46499% | 100.35860% | 0.51513% |

A small nodal trace is not sufficient evidence of interface accuracy. Strong
fixes the mesh trace, while pressure at a point 0.5 m below it is nonzero in the
continuum solution. Off-grid or sloping water tables retain a stair-step position
error, and a retained VT terrain trace can add the previously documented strong
support approximation. Curved/coastal execution tests do not supply an exact
curved groundwater reference.

### Build, regression and transport

The isolated clean Horizon build passed with the default modules (`ml reset`).
Executable SHA-256 is `454882d865a51647a37f7b5fe512983d88b59787564007bf24f9cbf94e3be97a`;
the source manifest pins 106 files. The before binary is the preceding
geometry-default build (`afca6eb5ad04b5f7551b8f2ed4d5dbff303d078a0d0517f43165176067a6856e`).
Local extracted-production tests passed 42 reader/classifier UBSan checks and
77 groundwater-face quadrature checks, plus the existing exposed-face and
independent projected mass/flux/strong-support algebra checks. No local CUDA
build was used as a substitute for the target build.

The final audit accepts **79 completed successful executions**
(72 unique fixtures), 19 matched comparisons and
7 one-node CUDA-aware versus two-node host-MPI comparisons.
Five further fixtures intentionally reject negative depths, nonfinite spacing,
incomplete data, irregular spacing and an all-dry domain. Their enclosing Slurm
jobs report FAILED because rejection is expected; each exact diagnostic and
nonzero exit is checked separately. No partial output is counted as a pass.

Coverage includes both mass treatments; all hydraulic surface closures; FEM/VT
terrain with varying groundwater depth; mixed porous/elastic/CPML regions;
pressure and velocity-potential acoustics, open/no-flow coastal contacts;
nonlinear and BKT/BKT2/BKT3/BKT3F interiors; unchanged standalone elastic,
Rayleigh and acoustic paths; and compact, batched, owner-row and element mass
operators. The latter paths are confirmed selected in the monitor before
comparison. This does not independently validate nonlinear yielding,
hysteretic response, every attenuation law, viscous interface accuracy or
all geometry/constraint combinations.

- Before/after: maximum aggregate relative L2 differences v 0.000198134%, p 7.05385e-05%, q 1.20295e-05%.
- Zero-depth grid versus surface saturation: maximum aggregate relative L2 differences v 0.000283674%, p 7.57729e-05%, q 1.77001e-05%.
- Mass-path and balance-diagnostic controls: maximum aggregate relative L2 differences v 9.65225e-06%, p 1.91278e-06%, q 3.82262e-07%.
- MPI/backend controls: maximum aggregate relative L2 differences v 1.57777e-05%, p 1.20773e-05%, q 7.71215e-06%.

The weak boundary balance residual/scale L2 is `5.18722724e-11`;
its integrated penalty dissipation is nonnegative (1.7931817 J).
This is a pressure-update balance, not a whole-system CPML energy proof.
Dry receiver checks also cover a point physically above the table but inside
a whole-cell porous owner: pressure stays undefined and no fluid companion
is emitted. The no-station/no-acoustic CPML case verifies that the grid is
loaded before the reference-speed/material scan.

### Ten-second CPML controls

The 256 m cube retains a 64 m-deep table, 32 m CPML, 16 Hz point-source rate,
consistent mass and gamma 0.0025. Uniform h=8 m and dt=100 microseconds are
verified from the actual logs. Ratios below compare each receiver's last
0.5 s peak to its first 0.5 s peak, then take the maximum over receivers.
These are boundedness controls, not a slow-wave accuracy or universal stability
claim; no enlarged-domain absorber error is inferred from them.

| Enforcement | Late/early v | Late/early p | Late/early q |
| --- | ---: | ---: | ---: |
| strong | 5.39354e-05 | 4.66867e-07 | 1.90517e-06 |
| weak | 5.12721e-05 | 4.95872e-07 | 2.1499e-06 |

All jobs used `cuda-615-SC`; no existing reservation was changed. Jobs: `96806`, `96807`, `96808`, `96809`, `96810`, `96811`, `96814`, `96815`, `96817`, `96821`, `96822`, `96823`, `96824`.
The account was not authorized for `IO_testing` at submission. Fixtures, source
manifests, independent reference, raw station outputs and detailed audit/plot
helpers remain local-only,
unavailable in a clean checkout. No existing reference or benchmark input
was changed to obtain agreement.


## BrickUP fixed groundwater through CPML — 2026-11-17

This extends `surface_groundwatertable=yes` / `groundwaterlevel_file` to
`poroelastic_formulation=brick_up` with regular porous CPML. It supersedes the
surface-groundwater-only restriction in earlier dated entries. No parameter or
default changes are introduced. This is reduced-Biot implementation/absorber
validation, not a Wolf full-Biot reproduction.

The existing center-based saturation mask continues through the layer. Dry cap
cells use solid CPML and saturated cells retain pressure storage, coupling and
Darcy diffusion. BrickUP keeps its existing locally horizontal drained-plane
compliance lift, including the original nodal-volume-weighted 2x2 quadrature.
The lift's direct coefficient and memory tail now use the same Cartesian J
filter as pressure storage. This also corrects surface-saturated drained
BrickUP CPML. Without this correction the boundary/storage ratio changes
through the stretch. Full Biot retains its distinct mesh-face hydraulic law.
No pressure histories, global solve, communication, allocation or transfer is
added per step; eight coefficient doubles per porous CPML cell are populated
at setup. Full-Biot kernels compile out the reduced drainage correction.

Native Cartesian stiffness, no checkpoint/restart, no geostatic ramp, and
regular VT-buffer requirements remain. Full-Biot open-pore constraint guards
are untouched. Whole-cell saturation is still an approximation; no moving
table, unsaturated/capillary flow, groundwater subcell volume or initial
hydrostatic pressure is introduced. Properties above the table remain the
supplied single-phase solid values. Nonlinear effective-skeleton forces and
attenuation retain their existing eligibility, with nonlinear/BKT updates
excluded from CPML cells.

### Provenance and checks

The clean Horizon build used default modules (`ml reset`), existing GSL/PROJ
libraries, and 106 frozen production/build files. Binary SHA256:
`37efcfa63ae5490b0e98c4c7dc812a096be76740a52f6f1c56621877b0731974`.
The unchanged baseline was full-Biot groundwater commit `990e21b`, binary
`454882d865a51647a37f7b5fe512983d88b59787564007bf24f9cbf94e3be97a`.
All production sources matched the frozen build manifest at final audit.

Local host algebra checks exercised 819,200 storage-plus-drainage recurrence
balances against an independently evaluated complete J transfer. Scaled error
was below 7.9e-15. Cases include no damping, non-unit direct stretches, side/
edge/corner profiles, repeated and zero poles, and zero penalty. Extracted
production setup passed 8,000 nodal quadrature checks in each of double and
single precision; drainage weights matched independent integrals within
5.8e-16 relative error. Non-CPML assembly agrees with the earlier routine
within roundoff. Existing CPML coupling/storage/inertia algebra also passed.
These host checks do not replace target GPU/MPI propagation checks.

Horizon jobs **97112–97116, 97118–97119, 97125, 97127–97128, 97130–97131,
97134** all finished on `cuda-615-SC`, using one, two or four nodes. IO_testing
did not list this account, so no jobs used it. The suite audited **65 unique
successful fixtures / 68 successful executions**, **12 matched regression
pairs**, **three host/CUDA-aware MPI comparisons**, and **five explicit
rejections**. No requested case remains pending and no unresolved execution or
regression failure remains. Absorber accuracy limitations are recorded below. Slurm FAILED for 97114/97130/97134 includes the
intentional rejection cases and, for 97114, the invalid fixtures below.

Eight initial coastal fixtures from an older investigation were rejected
before solver initialization because their wet surface slopes violate the
current geometry contract. These failed attempts remain archived and are not
passes. Current compatible shoreline fixtures replaced them. An additional
eight short-source coastal controls retained an existing coupled explosion;
eight full 0.3 s force-pulse controls provide the stronger excited-wave checks.
The checkpoint test encountered the parser's earlier absorber-history guard;
the first geostatic test encountered the earlier sealed nonlinear-ramp guard.
A permeable follow-up reached and confirmed the specific CPML geostatic guard.

Execution coverage includes flat/FEM/compatible VT buffers; mesh-aligned,
inside-cell and varying tables; drained/undrained, zero penalty, sealed and
anisotropic permeability; regular CPML profiles; dry/porous solid contacts;
acoustic open/no-flow with pressure and velocity-potential formulations;
nonlinear and Rayleigh-damped interiors; and unchanged elastic, viscoelastic,
acoustic, acoustic–elastic and full-Biot groundwater controls. The terrain and
coupled cases are execution/regression evidence, not a continuum accuracy
qualification of arbitrary slopes, shorelines or constitutive laws.

Zero-depth grids reproduce surface saturation; zero penalty reproduces the
undrained pressure update. Maximum aggregate relative L2 differences among
all 12 matched pairs were 7.20e-7 velocity, 9.50e-7 pressure and 7.04e-7 full-Biot
flux (the maxima are in the unchanged damped full-Biot control). Two-node host
versus one-node CUDA-aware results differed by at most 1.44e-10 velocity and
6.48e-12 pressure in buried drained/undrained and varying-table cases. Recorded
fields, station pressure eligibility, source/input hashes, runtime settings,
completion and binary provenance were audited. Full-Biot true residuals
retained their configured tolerance.

### Absorption and near-table pressure

The finite-Darcy case uses the earlier homogeneous reduced-Biot material:
mixture density 1664 kg/m3, drained Vp 3296.268286 m/s, Vs 1961.161351 m/s,
porosity 0.4, alpha 0.522666667, fluid bulk modulus 2.5e9 Pa, viscosity 1e-3 Pa s,
and permeability 6e-9 m2 in each direction. These are benchmark values, not
new defaults. A 25 Hz central-frequency Ricker point force acts 256 m below
ground. Uniform 8 m cells, dt=100 microseconds, and a 0.5 s record are retained
between matched small/large domains. The reference is a 2048 m cube with the
same ground/table, source and translated receivers, and sufficiently distant
side/bottom truncations for that record. The 16 m near-table control has its
own same-spacing large-domain reference. No filtering, gain fit or time shift
is applied. Values below are worst per-receiver relative L2 differences.

The buried table is at 67 m, inside an 8 m cell; a separate 64 m aligned table
was also tested. Bulk receivers include one at the compact model's inner CPML
boundary; exclude absorbing-layer fields from physical predictions. The
near-table set is horizontally 80 m from the source at table-relative depths
-8, 0.5, 8 and 32 m; the dry probe has no defined pore pressure. The wider
1024 m control translates source/receivers horizontally and preserves their
absolute depths, moving the lateral/bottom boundaries outward.

| Setup | Solid velocity | Pore pressure |
|---|---:|---:|
| 512 m cube, 128 m CPML, default profile | 22.01959% | 1.21934% |
| 512 m cube, 128 m CPML, comparison profile | 0.08379% | 0.07587% |
| 1024 m cube, 256 m CPML, default profile | 0.16713% | 0.07996% |
| Near-table probes, 8 m, compact default profile | 11.68406% | 2.51392% |
| Near-table probes, 16 m, compact default profile | 28.38585% | 24.76209% |
| Near-table probes, 8 m, comparison profile | 0.03811% | 0.06821% |
| Near-table probes, 8 m, wider default profile | 0.05545% | 0.16761% |

The default profile is R=0.01, order=2, kappa_max=5, alpha_max=0.5 Hz and
alpha_floor_fraction=0.1, with no multiaxial cross damping. The comparison
profile uses R=1e-4, order=3, kappa_max=1, alpha_max=6.25 Hz and zero alpha
floor, with all source/material/mesh/timestep inputs held fixed. This is an
explicit test profile, not a change to production defaults. The wider-buffer
control retains every default profile parameter.

The compact default-profile velocity error already exists in the surface
case: 22.03323% before and 22.03335% after the drainage correction, versus
22.01959% with the buried table. Its worst strictly interior bulk-receiver
velocity error is about 11.51%; the largest figure includes the inner-boundary
probe. These results do not qualify the compact default absorber as accurate.
The comparison profile and wider default buffer bring the buried case below
0.17% velocity and 0.08% pressure at bulk receivers; near-table maxima are
0.056% velocity and 0.168% pressure. They establish small numerical boundary
error for these tested configurations, not absolute continuum or drainage-law
accuracy. Near-table default-profile error decreases from 16 m to 8 m but
remains material; refinement and absorber distance/profile controls must be
assessed together.

Two 10 s buried-table controls retain the compact default profile, one drained
and one undrained. All outputs remain finite, with an overall decaying coda
and no sustained late growth. Worst last-0.5-s / first-0.5-s peaks are 2.066%
velocity and 0.03441% pressure for drained, and 2.075% / 0.03449% for undrained.
The coda is not monotone and the velocity tail is still appreciable; boundedness
does not cure the compact absorber's measured reflection error.

Raw inputs, diagnostics, histories and helper scripts are local-only, unavailable in
a clean
checkout. Their README records preparation, immutable build/run manifests,
analysis and overwrite behavior. Documentation lint, TypeScript, both site
builds/rendered audits, desktop rendering and 390 px parameter-link/layout
checks passed. No website publication, merge or tag was performed.


## Full-Biot four-curve groundwater CPML comparison — 2026-11-17

The maintainer requested the preceding four-curve BrickUP comparison with
`poroelastic_formulation=full_biot`. This study uses the unchanged solver at
`c9f9394`, retaining the same source history, uniform 8 m mesh, 100 microsecond
timestep, 0.5 s record, material and relative source/receiver locations.
The force Ricker central frequency is 25 Hz; `simulation_wave_max_freq_hz`
remains 24.271798903240597 Hz as the existing uniform-mesh selector.
This is a same-grid numerical boundary comparison, not a Wolf reproduction
or a continuum-accuracy/slow-wave-resolution qualification.

### Four runs and current defaults

| Figure color | Total box dimensions | Boundary |
|---|---|---|
| Black | 2048 x 2048 x 2048 m | Distant dashpots; same-grid full-Biot reference |
| Red | 512 x 512 x 512 m | 128 m side/bottom CPML, default profile |
| Blue | 1024 x 1024 x 1024 m | 256 m side/bottom CPML, default profile |
| Orange | 512 x 512 x 512 m | 128 m side/bottom CPML, comparison profile |

Box dimensions include the absorber. The reference's inherited thickness key
is inactive with `abc_type=dashpot`. The blue case changes both boundary
location and thickness; the orange case changes only profile settings.
The free ground surface and absolute source/receiver depths are unchanged.
Source X/Y and receiver X/Y translate together as the box widens.

Default profile: reflection target 0.01, order 2, kappa_max 5, alpha_max 0.5 Hz,
alpha_floor_fraction 0.1; comparison profile: 1e-4, 3, 1, 6.25 Hz, 0 respectively.
Both have zero multiaxial cross damping and no CPML material damping.
The requested reflection target is not the measured error.

The inputs deliberately omit mass treatment, stabilization, drained enforcement
and solver tolerance. Runtime logs verify **consistent mass, gamma 0.0025,
strong drained enforcement for flat geometry, and PCG tolerance 1e-6**.
The production optimized mass controls remain in use: element-inverse
preconditioning, quadratic predictor with cached products, fused local kernels,
combined solid/fluid packets, sequential solves, and full vectors.
Physical porous cells use consistent mass/storage; CPML and dry cells retain
their existing lumped contributions. No implementation or default was changed.

Material values match the preceding reduced case: mixture density 1664 kg/m3,
drained Vp 3296.268286 m/s, Vs 1961.161351 m/s, porosity 0.4, alpha 0.522666667,
fluid density 1040 kg/m3, fluid bulk modulus 2.5e9 Pa, viscosity 1e-3 Pa s,
permeability 6e-9 m2 in all directions, tortuosity 2. The viscosity is nonzero;
this is not the inviscid limit. Independent fluid inertia is active. Runtime
logs confirm two SDIRK inertia solves and two scalar storage solves per step.

The fixed physical groundwater table is 67 m below flat ground. The source is
256 m deep. Four receivers, 80 m horizontally from the source, lie at table-relative
depths -8, 0.5, 8, 32 m. The dry receiver has undefined pore pressure and no fluid
companion. With 8 m center-based saturation, the first porous mesh row starts at
64 m. Full-Biot strong drainage fixes pressure at that discrete dry-cap mesh
boundary; it does not impose a subcell boundary at 67 m. This differs from
BrickUP's cut-plane compliance, so the two formulations are not expected to
produce the same waveforms. Each is compared with its own larger-domain run.

### Completion and numerical differences

Horizon jobs 97159 (two compact cases, one node/four ranks), 97160 (wider case, two
nodes/eight ranks), 97161 (reference, four nodes/16 ranks) completed on cuda-615-SC.
Default modules were reset with `ml reset`; required GSL/PROJ paths and binary
libraries were verified. IO_testing did not list this account. Existing
reservations were not changed. All four runs completed 5000 steps; solid and
fluid histories were finite, with expected pressure/fluid eligibility.
Exact input/source/runtime/binary provenance was audited. The identical source
history SHA256 is
`4a4b97038b4cb89bc1a913c9f4e89bfe3e3e6f17e6391376953bc5020bf30a98`.
Each run reported 10000 inertia and 10000 storage solves, preserving fresh
true-residual acceptance. Maximum reported true relative residual was
`9.99967967e-07`, below 1e-6.

At the plotted receiver, 0.5 m below the physical table, the following unfiltered
relative L2 differences use the corresponding full-record reference norm.
Velocity and flux columns here are the plotted vertical components.

| Candidate | Vertical solid velocity | Pore pressure | Vertical Darcy flux |
|---|---:|---:|---:|
| 512 m / default | 6.09457% | 3.52128% | 5.40005% |
| 512 m / comparison profile | 2.89657% | 3.06164% | 2.72994% |
| 1024 m / default | 2.30477% | 0.96289% | 1.97311% |

Worst per-receiver differences across all available fields (four velocity
receivers, three pressure/flux receivers) use vector norms for velocity/flux:

| Candidate | Solid velocity vector | Pore pressure | Darcy-flux vector |
|---|---:|---:|---:|
| 512 m / default | 5.84308% | 3.63317% | 5.96491% |
| 512 m / comparison profile | 2.82521% | 3.23412% | 2.87084% |
| 1024 m / default | 2.46215% | 0.96289% | 2.40179% |

At the displayed receiver the wider default case has the smallest full-record
velocity/pressure/flux differences. The modified compact profile reduces the
late velocity coda relative to the compact default case, but retains roughly
3% pressure error. The profile change therefore does not reproduce the very
small errors seen in the BrickUP comparison. No single profile parameter is
isolated: the orange case changes five together. These four tests also do not
isolate the effects of consistent-to-lumped mass transition, stabilization or
hydraulic closure. Further controls would be needed to attribute the remaining
error to any one mechanism.

The primary figure shows velocity/pressure above and differences divided by
the full-record reference peak below. A companion figure shows vertical
relative-fluid flux q=d[porosity*(Uf-u)]/dt. No fitted time shift, gain,
filtering or resampling is used. Velocity/flux have the same half-step
stagger in all four records. Their cross-domain differences therefore use
identical recorded timestamps.

The reference keeps the same physical free surface; that physical reflection
is included in all curves. Its continuum fast-wave estimate gives an earliest
lateral truncation return beyond 0.54 s for these receivers, outside the 0.5 s
record. This does not remove common discretization, boundary representation,
stabilization or tolerance error. This study does not establish 10 s stability,
arbitrary terrain/coastal accuracy, or a universal CPML parameter choice.
No production default was changed on the basis of these four runs.

Binary SHA256:
`37efcfa63ae5490b0e98c4c7dc812a096be76740a52f6f1c56621877b0731974`.
The reused clean build matches all 106 frozen production/build source hashes.
Inputs, job manifests, raw records, plots and audit helpers remain local-only; they
are unavailable
in a clean checkout. The local README documents dependencies, inputs,
units, commands, overwrite behavior and limitations. No merge, tag or website
publication was performed.

## Surface-saturated BrickUP and full-Biot CPML comparison — 2026-11-17

The maintainer requested the same four-curve comparison for both formulations
without a buried groundwater table. All eight cases set
`surface_groundwatertable=no`, omit `groundwaterlevel_file`, and use porosity
0.4 from flat ground Z=0 downward. There is no dry cap or acoustic water layer.
The unchanged solver/binary at `c9f9394` was reused after matching all 106 frozen
production/build source hashes. No new implementation or default was required.

### Matched setup and four curves

| Figure color | Total model dimensions | External boundary |
|---|---|---|
| Black | 2048 x 2048 x 2048 m | Distant dashpots; separate reference for each formulation |
| Red | 512 x 512 x 512 m | 128 m side/bottom CPML, default profile |
| Blue | 1024 x 1024 x 1024 m | 256 m side/bottom CPML, default profile |
| Orange | 512 x 512 x 512 m | 128 m side/bottom CPML, comparison profile |

Dimensions include CPML. There is no top CPML: the physical ground surface
remains free and drained in every run. The wider case changes both boundary
distance and thickness; the orange case changes only the profile. Default
profile values are reflection target 0.01, polynomial order 2, kappa_max 5,
alpha_max 0.5 Hz, alpha_floor_fraction 0.1. The comparison profile uses 1e-4,
3, 1, 6.25 Hz, 0 respectively. Both have zero multiaxial cross damping and no
CPML material damping. The reflection target is not a measured error.

The source is 256 m below ground, with the same byte-identical 25 Hz Ricker
force history as the buried-table comparisons. Its three force amplitudes
are (1e9, 4e8, -2e8) N in solver X/Y/Z. The central frequency is 25 Hz;
`simulation_wave_max_freq_hz=24.271798903240597` remains the existing uniform
8 m mesh selector. All runs use dt=100 microseconds, 5000 steps and a 0.5 s
simulation. Station records contain 2500 samples at 200 microsecond spacing,
from 0 through 0.4998 s. Source/receiver X/Y translate together with box size.
Receivers are 80 m horizontally from the source at ground-relative depths
0, 0.5, 8 and 32 m. The Z=0 receiver replaces the prior dry-cap probe.

Material values remain mixture density 1664 kg/m3, drained Vp 3296.268286 m/s,
Vs 1961.161351 m/s, Biot alpha 0.522666667, fluid density 1040 kg/m3, fluid
bulk modulus 2.5e9 Pa, viscosity 1e-3 Pa s and permeability 6e-9 m2 in each
direction. Full Biot uses tortuosity 2 and independent fluid inertia. These
are finite-viscosity cases, not the inviscid limit.

BrickUP retains its lumped update and existing drained compliance factor 1;
full-Biot stabilization gamma does not apply to that formulation. Full Biot
uses the omitted-input defaults: consistent mass/storage in physical porous
cells, strong drained enforcement for flat geometry, gamma 0.0025, PCG
tolerance 1e-6 and the existing optimized mass controls. CPML retains its
existing lumped contributions. Each formulation is compared with its own
reference, not with the other formulation or a continuum analytical solution.

### Completion and audited results

All six Horizon jobs completed with exit code 0 under `cuda-615-SC`:

| Job | Cases | Nodes / MPI ranks | Job elapsed |
|---|---|---:|---:|
| 97169 | BrickUP compact default and comparison profile | 1 / 4 | 37 s |
| 97170 | BrickUP wider default | 2 / 8 | 33 s |
| 97171 | BrickUP reference | 4 / 16 | 39 s |
| 97172 | Full-Biot compact default and comparison profile | 1 / 4 | 102 s |
| 97175 | Full-Biot wider default | 2 / 8 | 98 s |
| 97176 | Full-Biot reference | 4 / 16 | 129 s |

The two compact cases execute sequentially within their respective jobs;
these job elapsed times are not comparable per-step solver benchmarks.
The launcher used `ml reset` and verified the binary's required libraries.
No reservation was altered. IO_testing did not list the account.

All eight runs passed exact input/source/binary/runtime provenance checks,
completed 5000 steps, and produced finite pressure/solid histories at all four
porous receivers. Full Biot additionally produced all four fluid companions;
BrickUP correctly produced none. Runtime logs verify saturation from the local
ground surface and the selected formulation/boundary defaults. Every full-Biot
run reported 10000 inertia and 10000 storage solves; the largest reported
true relative residual was 9.99981967e-7, below 1e-6.

At the plotted receiver, **0.5 m below ground**, relative L2 differences over
the full record are:

| Formulation | Candidate | Vertical solid velocity | Pore pressure | Vertical Darcy flux |
|---|---|---:|---:|---:|
| BrickUP | 512 m / default | 9.01869% | 3.46398% | Not an independent field |
| BrickUP | 512 m / comparison profile | 0.07897% | 0.04064% | Not an independent field |
| BrickUP | 1024 m / default | 0.03921% | 0.10001% | Not an independent field |
| Full Biot | 512 m / default | 5.68173% | 3.52853% | 4.38278% |
| Full Biot | 512 m / comparison profile | 3.62676% | 2.61480% | 3.10767% |
| Full Biot | 1024 m / default | 2.06668% | 1.27687% | 1.45446% |

Worst differences across the four receivers use vector norms for solid
velocity and Darcy flux:

| Formulation | Candidate | Solid velocity vector | Pore pressure | Darcy-flux vector |
|---|---|---:|---:|---:|
| BrickUP | 512 m / default | 10.21738% | 4.15158% | Not an independent field |
| BrickUP | 512 m / comparison profile | 0.05238% | 0.06102% | Not an independent field |
| BrickUP | 1024 m / default | 0.04642% | 0.10065% | Not an independent field |
| Full Biot | 512 m / default | 6.18786% | 3.52853% | 6.30084% |
| Full Biot | 512 m / comparison profile | 3.23172% | 3.20594% | 3.30479% |
| Full Biot | 1024 m / default | 2.23097% | 1.27687% | 2.33831% |

Full-Biot pressure is identically zero at the surface probe in all four runs.
Its relative surface-pressure error is therefore undefined (zero reference
norm), and is excluded from relative-pressure maxima; its absolute difference
is exactly zero. BrickUP's weak compliance correction leaves a nonzero surface
pressure: the reference peak is 165.2541 Pa. Consequently the very small
BrickUP CPML differences do not establish exact drained-boundary accuracy;
all four BrickUP runs share that boundary approximation.

Both formulations reproduce the earlier qualitative trend: the compact
default profile leaves a larger late coda; the modified compact profile
substantially reduces it. BrickUP then closely matches its reference, whereas
full Biot retains percent-level differences. The wider default case has the
smallest full-record full-Biot differences at the displayed receiver, while
the modified compact case has less late coda there. These controls do not
isolate the cause of the remaining full-Biot error (for example the
consistent-to-lumped transition, stabilization or hydraulic treatment).

The two principal figures show vertical solid velocity and pore pressure,
with difference/reference-peak panels below. A companion figure shows full-Biot
vertical Darcy flux. Histories are unfiltered and use identical recorded
timestamps, without fitted shifts, gains or resampling. The reference keeps
the same physical free-surface reflection. A continuum fast-wave estimate
places its earliest lateral truncation return beyond 0.55 s for full Biot,
outside the record. Common mesh dispersion, surface approximation and
stabilization/tolerance errors are not removed by this numerical comparison.
These 0.5 s runs do not establish long-time stability or general terrain,
elastic/acoustic-interface, viscoelastic or nonlinear coupling accuracy.
No solver changes or default-profile changes were made.

Binary SHA256:
`37efcfa63ae5490b0e98c4c7dc812a096be76740a52f6f1c56621877b0731974`.
Source-history SHA256:
`4a4b97038b4cb89bc1a913c9f4e89bfe3e3e6f17e6391376953bc5020bf30a98`.
Inputs, manifests, raw output, analysis and PNG/PDF figures remain local-only
in the private validation archive, unavailable in a
clean checkout. Its README records dependencies, commands, units, overwrite
behavior and limits. Documentation lint, TypeScript, both site builds and their
rendered/parameter audits passed. Merge/tag remains on hold; no website
publication occurred.

## Consistent full-Biot CPML extension (2026-11-17)

This change extends `poroelastic_mass_matrix=consistent` from physical porous
cells into regular full-Biot CPML cells. It removes the previous
consistent-interior/lumped-layer transition. No new parameter is needed; the
existing conditional full-Biot default selects it. Explicit lumped full Biot,
BrickUP, elastic and scalar acoustic CPML retain their schemes. The comparisons
in the preceding sections describe the earlier implementation.

### Discrete operators and retained optimizations

For element-constant CPML coefficients, the exact Q1 shape mass integrates the
complete stretched solid/fluid inertial block, including cross inertia, the
polynomial velocity/restoring terms and convolution memories. The same Q1
integration pairs leading pressure storage with its convolution tail. Spatial
integration of the memories uses a separable three-axis shape-mass application,
without new history variables or field exchanges. Physical solid mass damping
and Darcy drag retain their existing lumped operators and matching lumped CPML
histories. Leading Darcy drag and physical-interior solid mass damping retain
the existing SDIRK2 stages. Optional solid mass damping inside CPML keeps its
lumped backward difference and matching histories; base CPML polynomial
damping keeps its explicit time difference with consistent spatial integration.

The physical domain and layer use the same projected distributed PCG, default
1e-6 tolerance, element-inverse preconditioner, quadratic predictors, cached
products, fresh residual verification, fused kernels and combined solid/fluid
CUDA-aware or host MPI packets. Uniform Q1 assembled rows now have separate
inertia/storage coefficients where they touch CPML; physical-only rows retain
byte integer weights. Weighted rows add 432 bytes per local row. General
material/mesh and native VT physical cells use the corresponding general element
operator. Optional compact vectors, solve batching, single-reduction CG and
aggregate packing remain available. Experimental owner-complete one-round rows
explicitly fall back to the existing two-exchange operator for weighted porous
CPML. They remain off by default; this extension does not claim that experiment
has been generalized to weighted rows.

The initialization subtracts the complete stretched lumped contribution before
inserting the consistent operator, so shared physical/CPML and porous/nonporous
nodes retain the correct sum. Elastic CPML skips its former solid inertia
history force on these cells; the coupled kernel supplies it once, together
with fluid and cross terms. Its stiffness contribution is unchanged. Strong,
weak and natural hydraulic boundary operators, source normalization and station
quantities retain their existing contracts.

Existing restrictions remain: CPML needs a regular Cartesian buffer, native VT
cut cells must stay outside it, and full-Biot open-pore CPML coupling requires
field-specific constraints. Nonlinear elastoplastic/plastic and hysteretic
constitutive updates and BKT-family attenuation are excluded inside the linear
CPML buffer; applicable physical interiors remain supported. Fixed groundwater
classification remains whole-cell, not a subcell interface. Checkpoint/restart,
geostatic ramp, fixed-base and DRM restrictions are unchanged. This does not
relax the earlier shared-open-pore stability guard.

### Matched absorption comparisons

These are 25 Hz force-source CPML controls, not new Wolf analytical comparisons.
They repeat the preceding surface-saturated and buried-table tests at h=8 m,
dt=100 microseconds, duration 0.5 s, gamma=0.0025 and PCG tolerance 1e-6. Total
compact dimensions are 512 m cubed with 128 m side/bottom CPML; wider dimensions
are 1024 m cubed with 256 m CPML. Each reference is 2048 m cubed with distant
dashpots. The physical ground is free, with no top CPML. The buried table is
67 m below ground, using the existing whole-cell saturation geometry. Source
history and relative source/receiver geometry are unchanged.

The default profile uses reflection target 0.01, order 2, kappa_max 5,
alpha_max 0.5 Hz and alpha_floor_fraction 0.1. The modified profile uses 1e-4,
3, 1, 6.25 Hz and 0. Both have zero multiaxial cross damping. These inputs are
unchanged from the matched earlier tests. No profile default was changed.

The later [alpha-only follow-up](CPML_ALPHA_ONLY_VALIDATION.md) changes only
alpha while holding the other controls fixed. Increasing alpha alone does not
produce the improvement below; reducing it gives only a small gain in those
matched controls. The multi-parameter results here do not isolate alpha's effect.

At the receiver 0.5 m below the hydraulic surface, full-record relative L2
errors against the matching larger-domain reference are:

| Hydraulic surface | CPML case | Previous velocity Z | Consistent velocity Z | Previous pressure | Consistent pressure |
|---|---|---:|---:|---:|---:|
| Ground surface | Compact, default profile | 5.68173% | 0.88480% | 3.52853% | 1.21231% |
| Ground surface | Compact, modified profile | 3.62676% | 0.05809% | 2.61480% | 0.02376% |
| Ground surface | Wider, default profile | 2.06668% | 0.23227% | 1.27687% | 0.27985% |
| Buried table | Compact, default profile | 6.09457% | 3.25906% | 3.52128% | 1.37452% |
| Buried table | Compact, modified profile | 2.89657% | 0.04829% | 3.06164% | 0.04436% |
| Buried table | Wider, default profile | 2.30477% | 0.21288% | 0.96289% | 0.16191% |

The viscous material retains eta=1e-3 Pa s and permeability 6e-9 m2 in all three
directions. Across every receiver and the tested solid-velocity, pressure and
Darcy-flux vector/vertical fields, the modified-profile maximum relative L2
errors are 0.07220% for saturation from the ground and 0.08294% for the buried
table. Zero surface-pressure reference norms are excluded from relative errors;
the actual strong surface trace is exactly zero. The repeated non-CPML surface
reference differs from the previous reference by at most 6.55e-8 relative L2.

Additional ground-surface controls set only viscosity to zero and use a matching
inviscid larger-domain reference. The new consistent compact/default case has
1.85358% vertical-velocity and 10.86614% pressure error at the shallow receiver;
the modified profile reduces these to 0.09441% and 0.35247%. The worst modified
inviscid relative error over all tested fields/receivers is 0.52922%. Consistent
mass improves matching but does not make every CPML profile equally effective.
The inviscid dashpot reference still uses two inertia stages because its distant
physical boundary damping is active; the inviscid CPML controls use one.

All differences use unfiltered histories with identical timestamps, without
resampling, fitted shifts or gains. These numerical references isolate boundary
and discretization-transition differences at fixed resolution; common physical
mesh dispersion, stabilization and boundary-approximation errors remain. The
source spectrum and 8 m resolution are not a claim of continuum accuracy.

### Long-time controls

The viscous and inviscid surface-saturated compact/default cases completed
10 s (100000 steps) with finite solid, fluid and pressure records. Half-second
receiver envelopes decay overall; the inviscid tail fluctuates rather than
decaying monotonically. In the last half-second, the largest normalized peak
across all tested receivers and fields is 1.91e-11 for the viscous case and
3.33e-4 (0.0333%) for the inviscid case, relative to each field's initial
half-second peak. No growing receiver envelope was observed. These measurements
do not establish unconditional stability, later-time behavior or long-time
stability of the coastal acoustic interface.

### Reproducibility and scope

The frozen build contains 106 source/build files, checked by SHA256 before a
clean Horizon CUDA/MPI build using `ml reset`; dynamic library resolution passed.
Candidate binary SHA256:
`f97c64e81012b098a29d7eeb76feb7092471515ab22a0d9d1e85f3c4accf2955`.
Old regression binary SHA256:
`37efcfa63ae5490b0e98c4c7dc812a096be76740a52f6f1c56621877b0731974`.
Source-history SHA256 for the matched surface/buried force controls:
`4a4b97038b4cb89bc1a913c9f4e89bfe3e3e6f17e6391376953bc5020bf30a98`.

An independent host check uses eight-point Gauss quadrature and separately
evolved histories to check the production Q1 inertia/storage updates for zero,
one, two and three stretched directions, zero/repeated poles, non-unit kappa,
physical solid damping and zero/nonzero Darcy drag. Float and double histories
passed, as did constant row sums, with undefined-behavior instrumentation.
A clean GPU build and real MPI runs provide separate evidence beyond that check.

All investigation inputs, helpers, manifests, raw outputs, analysis and PNG/PDF
figures remain local-only.
The local README documents requirements, commands, units, overwrite behavior and
limitations. The three figures are `consistent_cpml_surface`,
`consistent_cpml_buried` and `consistent_cpml_long_time`, each as PNG and PDF.

Historical 32-node performance gains and the weight-1 selection measured
consistent interiors with lumped CPML. They are not timings or a partition-weight
calibration of this extension. No universal speedup, optimal CPML profile or
production partition weight is inferred here. Default performance switches and
the global CPML partition weight remain unchanged.

### Completed integration and regression checks

There are 57 successful executions of 52 distinct accepted fixtures,
plus one documented input rejection before time stepping. Every accepted run
passed exact input/source/binary/runtime provenance, step and sample counts,
finite motion/fluid/eligible pressure fields, exactly zero relative-fluid flux
in sealed directions, and fresh PCG-residual checks.
The largest reported true relative residual was 9.99999038e-07, below 1e-6.
The parser truncates floating step counts: for example 0.3 s / 25 microseconds
produces 11999 steps, and 0.3 s / 100 microseconds produces 2999. The audit
matches that existing behavior, without changing input times or solver code.

All 18 submitted allocations used `cuda-615-SC`, with four GPU MPI ranks per
node, `ml reset` and verified binary libraries. No reservation was changed.
The jobs and their case groups are recorded below. Job 97246 has a failed job
exit because its first BKT2 fixture lacked two mandatory Q-control flags;
parsing stopped before stepping. Its other two cases completed normally.
A fresh immutable fixture with `use_parametric_q_factor=no` and
`use_infinite_qk=no` passed in 97260. This is an input-preparation correction,
not a solver fix or a successful run of the rejected fixture.

| Job | Cases | Nodes / ranks | MPI backend |
|---|---|---:|---|
| 97222 | `surface_compact_default`, `surface_compact_profile`, `surface_inviscid` | 1 / 4 | `cuda_aware` |
| 97223 | `surface_reference` | 4 / 16 | `cuda_aware` |
| 97224 | `surface_wide_default`, `buried_wide_default` | 2 / 8 | `cuda_aware` |
| 97225 | `brick_up_lumped_old`, `brick_up_lumped_new`, `full_biot_lumped_old`, `full_biot_lumped_new`, `elastic_old`, `elastic_new`, `acoustic_old`, `acoustic_new`, `acoustic_elastic_old`, `acoustic_elastic_new`, `viscoelastic_old`, `viscoelastic_new`, `brick_nonlinear_old`, `brick_nonlinear_new` | 1 / 4 | `cuda_aware` |
| 97226 | `surface_elements`, `surface_unfused`, `surface_compact`, `surface_batched`, `surface_single`, `surface_one_round` | 1 / 4 | `cuda_aware` |
| 97227 | `buried_compact_default`, `buried_compact_profile`, `full_varying`, `surface_sealed`, `surface_undrained`, `full_rayleigh_cpml`, `full_viscoelastic`, `full_nonlinear` | 1 / 4 | `cuda_aware` |
| 97228 | `fem_noflow`, `fem_open_p`, `fem_open_phi`, `fem_open_point`, `fem_weak` | 1 / 4 | `cuda_aware` |
| 97229 | `vt_noflow`, `vt_open_p`, `vt_open_phi`, `vt_weak` | 1 / 4 | `cuda_aware` |
| 97231 | `long_viscous` | 1 / 4 | `cuda_aware` |
| 97232 | `long_inviscid` | 1 / 4 | `cuda_aware` |
| 97233 | `surface_compact_default`, `fem_open_p` | 1 / 4 | `host` |
| 97234 | `surface_compact_default`, `vt_open_p` | 2 / 8 | `cuda_aware` |
| 97242 | `inviscid_reference` | 4 / 16 | `cuda_aware` |
| 97243 | `inviscid_profile`, `surface_pack` | 1 / 4 | `cuda_aware` |
| 97246 | `full_bkt2`, `surface_anisotropic`, `full_heterogeneous` | 1 / 4 | `cuda_aware` |
| 97260 | `full_bkt2_valid` | 1 / 4 | `cuda_aware` |
| 97263 | `vt_native_cut` | 1 / 4 | `cuda_aware` |
| 97264 | `vt_native_cut` | 2 / 8 | `host` |

The default row operator versus the general element operator, disabled
fusion/cache/combined packets, compact vectors, batched inertia/stabilization,
single-reduction CG, aggregate packing and the explicit one-round fallback
agree with the default case to at most 1.03623e-08 relative L2
across every tested field/receiver. Default versus general-element alone is
1.45101e-9. The matched old/new BrickUP lumped, full-Biot lumped, elastic,
acoustic, acoustic-elastic, Rayleigh and nonlinear controls differ by at most
2.29451e-10 relative L2. These small differences
include GPU reduction/atomic ordering and output rounding; bit identity is not
claimed.

Changing host/CUDA-aware MPI and/or four/eight ranks in the flat, FEM open-pore,
VT coastal and native-cut controls gives at most 1.10169e-06
relative L2 difference across tested fields. The ordinary coastal VT controls
contain zero retained cut cells; they exercise VT interface and constraint
paths but do not alone qualify native P1 integration. The additional native
case reports 268 retained cuts and 534 owned surface constraints, outside the
regular CPML buffer. Its 0.2 s, 5-microsecond undrained/no-flow acoustic-coupled
Biot-explosion record has an excited source peak at 0.04 s. Both communication
backends and rank layouts complete it with matching observables.

Other completed controls cover variable groundwater, partial/directional and
fully sealed permeability, drained/undrained hydraulics, FEM/VT weak enforcement,
pressure/velocity-potential acoustic formulations, open-pore/no-flow coupling,
FEM point projection and face-mass projection, a two-layer material forcing the
general mass operator, Rayleigh damping inside/outside CPML, BKT2 attenuation
in the physical interior and nonlinear solid interiors. Coastal records cover
0.3 s at 25 microseconds (FEM) or 5 microseconds (VT), retaining their 1e-8
PCG tolerance; the native-cut control explicitly uses 1e-6. The synthetic material
values are diagnostic, not production calibrations. Execution/finite-field and
MPI agreement checks do not establish analytical interface accuracy, nonlinear
constitutive calibration or long-time coastal stability. No claim is made that
these tests cover every combination of physics and boundary settings.

For context only, one-node/four-rank compact viscous Solver timers were
29.820 s with defaults,
29.801 s with the general element operator,
42.297 s with fusion/cache/combined packets disabled,
and 26.475 s with solve batching enabled.
These small runs used different allocations for the default and option group;
they establish neither a production speedup nor a new default choice. Fresh
32-node whole-step timing and partition-weight calibration remain outstanding.

Documentation lint, TypeScript, both site production builds and rendered/parameter
audits passed. A browser check verified the CPML prose, canonical mass-parameter
link and narrow-screen layout. No website publication, merge or tag was performed.

## Consistent full-Biot CPML production calibration (2026-11-17)

**Completed:** all six production records, four instrumentation controls and
two full-record profiles passed. Weight `1` is the fastest tested choice for this
homogeneous inviscid case. The conditional full-step GPU service-cost ratio is
about **1.87 CPML/interior**, not a replacement partition weight or global default.

The earlier 32-node weight selection used consistent physical cells and lumped
CPML. The current implementation also applies consistent Q1 inertia/storage and
matching memory terms inside porous CPML. Existing small-case logs do provide
useful complete-run timings, but every case below used weight 4.38, with detailed
GPU profiling disabled. They cannot determine the best production weight or
separate the complete CPML/interior computational cost.

| Matched viscous case | Nodes / GPU ranks | Old hybrid Solver (s) | Consistent CPML Solver (s) | Old/new ratio |
|---|---:|---:|---:|---:|
| 512 m cube, 128 m CPML, default profile | 1 / 4 | 34.571853 | 29.820085 | 1.15935 |
| 512 m cube, 128 m CPML, modified profile | 1 / 4 | 35.120696 | 29.381748 | 1.19532 |
| 1024 m cube, 256 m CPML, default profile | 2 / 8 | 70.077036 | 66.386666 | 1.05559 |
| 2048 m distant reference, dashpots | 4 / 16 | 95.960770 | 95.373759 | 1.00615 |

These are maximum-rank Solver timers for complete 0.5 s / 5000-step records,
with uniform 8 m cells and a 100-microsecond timestep. Parameter values and
source files match within each pair. The executions used different allocations;
no timing repeat was performed for these pairs. The distant reference contains
no CPML. None of these ratios is a 32-node production speedup or a material cost
ratio. Their numerical/reference scope remains that of the preceding surface
and buried-groundwater validation records.

### Matched full-size comparison

Job **97288** completed six Figure 5 records sequentially on the same
32 nodes / 128 GPU MPI ranks in `cuda-615-SC`: old hybrid at weight 1, consistent
CPML at weights 1, 2, 4.38 and 1.5, then a repeat at weight 1. Every case regenerates
its own mesh. The 4096 m cube uses uniform 4 m cells, a 256 m CPML at four lateral
faces and the bottom, the default absorption profile, the exact author 16 Hz
source-rate samples, the same eight receivers, inviscid material, timestep
200 microseconds, stabilization 0.0025 and mass tolerance 1e-6. Each record has
7500 steps over 1.5 s. Current default mass-solver performance switches are fixed,
and profiling is disabled. There are 1,073,741,824 cells, of which 303,038,464
are CPML. This is a homogeneous flat model without acoustic or nonlinear cells;
it does not calibrate mixed physics, topography, adaptive meshes or viscous
production runs.

| Completed case | Weight | Maximum-rank Solver (s) | Solver (ms/step) | Mean inertia / storage iterations |
|---|---:|---:|---:|---:|
| Old hybrid CPML | 1 | 817.262068 | 108.968276 | 3.979 / 3.695 |
| Consistent CPML | 1 | 773.032962 | 103.071062 | 3.029 / 2.476 |
| Consistent CPML | 2 | 821.149312 | 109.486575 | 3.029 / 2.476 |
| Consistent CPML | 4.38 | 1000.411781 | 133.388237 | 3.029 / 2.476 |
| Consistent CPML | 1.5 | 794.480000 | 105.930667 | 3.029 / 2.476 |
| Consistent CPML, repeat | 1 | 773.923964 | 103.189862 | 3.029 / 2.476 |

At weight 1, the new implementation is **1.057215 times faster** in complete
Solver time, a 5.41% reduction. The total program wall timers are 893.930269 and
850.448916 s, respectively; these include setup and differ from the full Slurm
allocation, which contains several cases. This comparison keeps timestep,
stabilization, source, mesh and tolerance fixed. It must not be substituted into
or multiplied by the historical speedup against the greater-than-seven-hour
50-microsecond / 1e-10 baseline.

The two weight-1 runs average **773.478463 s**, with a range of 0.891002 s
(0.115% of the mean). Relative to that mean, weights 1.5, 2 and 4.38 take
2.715%, 6.163% and 29.339% longer. Weight `1` is therefore the preferred tested
choice; the gap to the next tested choice exceeds the observed repeat range.
Only weight 1 was repeated, so this is not a statistical confidence interval.
The matched old/new ratio using the new mean is 1.056606, or 5.36% less Solver
time. Mean iteration counts are unchanged across the new weights. Their worst
recorded-field relative L2 difference from the first weight-1 run is 7.01e-12;
reference errors agree at the precision shown below.

Total program wall timers for weights 2, 4.38, 1.5 and the weight-1 repeat are
954.564117, 1167.542429, 919.266643 and 851.119433 s. These include mesh/setup
work and are distinct from Solver time. The complete six-case allocation took
1:35:08. No reservation was changed.

All six production cases passed binary/input/source/runtime provenance, full finite
motion/fluid/pressure histories and all 7500 inertia / 15000 storage solves with
fresh true residuals at or below 1e-6. The largest new-versus-old relative L2
change across recorded fields is 0.5496%. This compares different CPML spatial
mass treatments; it is not a roundoff-equivalence claim.

| Worst receiver reference L2 error | Old hybrid | Consistent CPML |
|---|---:|---:|
| Solid velocity, Karpfinger | 8.076834% | 8.080028% |
| Darcy flux, Karpfinger | 7.456355% | 7.456520% |
| Pressure, secondary SeisSol numerical reference | 7.641824% | 7.642259% |

Reference normalization and timing are fixed: Karpfinger scale 7.025e12,
velocity/flux at native time minus half a timestep, pressure at native time,
and the reference window 0.001–1.499 s. No gain or time shift is fitted. The
analysis reproduces the earlier job 95827 scores before accepting new results.
Pressure is a secondary numerical comparison, not an analytical pressure
solution. Similar reference errors do not establish spatial convergence or
unconditional long-time stability.

### Full-record GPU cost measurements

Job **97293** completed four matched 0.5 s viscous pulse controls on one node:
production executable, private profiler with events disabled, original-dispatch
event profiling, and separated physical/CPML traction and continuity profiling.
Their Solver timers were 29.584049, 29.565136, 55.840556 and 56.176625 s.
Every control passed provenance, complete finite fields and true-residual checks.
Worst waveform relative L2 difference from production was 7.79e-9; worst FAS
relative L2 difference was 5.57e-9. The small fixture qualifies the instrumentation;
its four ranks have nearly identical material proportions, so their composition
fits cannot provide a reliable CPML/interior cost ratio. Profiling overhead is
also visible; only unprofiled runs select the weight.

Job **97306** completed both profiling modes for the full 1.5 s Figure 5 record
on 32 nodes in `cuda-615-SC`, after successful completion of 97288 and 97293.
Original and split dispatch took 832.582530 and 832.430078 s Solver time; the
allocation took 31:15. Both passed provenance, finite complete histories and
all fresh true residuals at 1e-6. Worst relative waveform/FAS changes from the
unprofiled weight-1 run were 8.64e-12/8.88e-12. All three allocations exited 0.
The private profiler separately measures force/history kernels,
consistent operators, common nodal operations and packet helpers. It now counts
CPML cells in the consistent operator and records uniform versus CPML-adjacent
weighted rows; the earlier physical-only cost formula would omit that work.
Nested row timers are children of the operator timer and are not added twice.
CPU exchange/reduction waits overlap GPU work and are excluded from GPU service
costs. CPML face/edge/corner counts are collected at setup. Fixed event pools
are consumed at existing reduction waits, without additional step synchronization.
The model fits all 128 ranks' GPU service time over the complete record:

```math
T_r^{GPU}=c_{physical}N_{physical,r}+c_{CPML}N_{CPML,r}.
```

The reported coefficients divide by 7500 steps and convert to ns/cell/step.
They include consistent operators, common nodal work, stiffness, pressure-force,
continuity/history updates, stabilization and source work. Packet/constraint
helpers and output are measured separately; CPU waits are not added. Common
nodal and launch costs are amortized over cells in this fixed allocation, so
these coefficients are conditional costs, not intrinsic material constants.

| Full-record fit | Interior (ns/cell/step) | CPML (ns/cell/step) | CPML/interior | Relative RMS |
|---|---:|---:|---:|---:|
| Original production dispatch | 6.46664 | 12.06094 | **1.86510** | 0.7276% |
| Separated physical/CPML dispatch | 6.48394 | 12.11051 | **1.86777** | 0.6723% |

Both normalized design matrices have condition number 2.478. Separating the
launches changes aggregate measured GPU service by 0.328% and the per-rank
service vector by 0.352% relative L2; the two estimates agree closely. A richer
physical/CPML-face/edge/corner fit reduces relative RMS to 0.242%/0.215%, with
condition number 9.815 and composition-weighted CPML/interior ratios
1.82972/1.83474. The useful conclusion is **about 1.8–1.9 for this complete
inviscid workload**, not an exact universal ratio. The original-dispatch
four-region coefficients are 6.51896, 11.75692, 13.40535 and 13.36145
ns/cell/step; the small corner population does not establish a universal ordering
of edge versus corner costs.

The split-dispatch direct element kernels, covering stiffness and Biot
traction/continuity histories but excluding the mass/storage solves, stabilization
and common nodal work, cost 0.70077 ns per physical cell per step versus
5.09241 ns per CPML cell per step: **7.27 times**. This partial-kernel ratio is
not the full-step ratio above. The original-dispatch GPU service averages are:

| Component | Mean GPU ms/rank/step |
|---|---:|
| Consistent operators, including row children once | 26.9874 |
| Common nodal/fused work | 20.1790 |
| Physical stiffness | 0.9431 |
| Elastic CPML stiffness/history part | 2.4399 |
| Mixed Biot traction and continuity/history | 12.6684 |
| Stabilization spatial application | 4.2727 |
| Source work | 0.000085 |
| **Physics service total** | **67.4906** |
| Packet/constraint helpers, separate | 6.3056 |

These are rank means of GPU event intervals, not sums of per-stage rank maxima
or a decomposition of elapsed Solver time. Consistent operators plus common
nodal work account for 69.89% of the fitted GPU service. Directly timed
CPML-adjacent weighted mass rows take 0.19215 ns/node/application for inertia
and 0.13106 for storage, versus 0.16176 and 0.08889 for purely physical rows:
ratios 1.188 and 1.474. Row timings are diagnostic children of the consistent
operator total, not additional work. They include both operator and
preconditioner applications and are not standalone solve costs.

A summed GPU service ratio is not automatically the best synchronous partition
weight. Mass/storage solves, force stages and communication have separate rank
maxima; increasing the weight redistributes all ordinary work as well as CPML.
The measured complete-loop comparisons, including the weight-1 repeat, select
`1` for this case even though the fitted aggregate GPU cost ratio is about 1.87.
Mixed acoustic, elastic, reduced-Biot, attenuating or nonlinear models, viscous
production, topography and adaptive meshes require their own complete-run timing
and physics checks. The unprofiled comparisons, not event-instrumented timings,
are the performance evidence.

The new production binary has SHA256
`f97c64e81012b098a29d7eeb76feb7092471515ab22a0d9d1e85f3c4accf2955`
(source in `11ac0e1`); the old hybrid is
`37efcfa63ae5490b0e98c4c7dc812a096be76740a52f6f1c56621877b0731974`
(source in `c9f9394`). The private diagnostic binary is
`6868460dac66907b296078bf21f32bcef67ce287bceea47b2f4009c999177427`.
Its isolated clean CUDA/MPI build used `ml reset`, a 107-file source manifest
and verified dynamic-library resolution. Production solver source is unchanged
by this measurement study. No existing reservation was altered.

Inputs, diagnostic helpers, frozen source, raw results, audits and figures remain
ignored and local-only.
The local README documents interfaces, dependencies, units, overwrite
behavior and limitations. `production_timings`, `production_waveforms` and
`production_cost` figures are retained there as PNG/PDF; they are not assets of
a clean checkout. The global `cpml_partition_weight=4.38` and all physical and
solver defaults remain unchanged. This homogeneous calibration cannot establish
a mixed-material global default.

The frozen dependency install, lint, TypeScript check, both production builds
and rendered/parameter audits passed. Browser checks verified the calibration
link card, canonical parameter links and narrow-screen layout. Production solver
source is unchanged; the isolated diagnostic CUDA/MPI build and the twelve
completed Horizon records supply the runtime evidence. No website publication,
merge or tag was performed.


## Material-dependent CPML partition defaults (2026-11-17)

At the maintainer's request, omitted `cpml_partition_weight` now assigns cost
`1` to saturated poroelastic CPML and `4.38` to elastic/acoustic CPML, including
dry caps. Ordinary cells keep cost `1`. Both `brick_up` and `full_biot` use this
rule, with any supported mass option. An explicitly supplied number still
sets one uniform weight for **all** CPML cells; explicit `1` retains ordinary
equal-leaf-count partitioning. Omit the key to select material-dependent defaults.
The earlier statements that the global default remained `4.38` describe those
historical builds; this section records the subsequent interface change.

Partitioning and solver enrollment now share the same completed-cell material
and saturation classifier, including physical center coordinates, mesh origins,
buried groundwater, terrain and the porous depth limit. The groundwater grid
is loaded before completed-mesh partitioning. The explicit/omitted distinction
is broadcast to every solver rank. Construction and mesh-only partitions remain
ordinary. Boundary/interface operators, time integration, pressure/stabilization,
GPU kernels, MPI exchange mechanisms and station quantities are unchanged.
The regular porous CPML buffer and existing coupling/nonlinear restrictions remain.

A clean Horizon build used `ml reset` and the default toolchain; its binary SHA256
was `b71b5d1628c57f4905e76bb16d11001c4ba854820b939755721c7ceba21109a1`. The local frozen snapshot contains 106 production/build files.
Jobs **97464**, **97465**, **97466** and **97470** used `cuda-615-SC`, each with one
node and four GPU ranks. The controls retained existing pulse sources, material
values, mesh, timestep and CPML profiles, shortening duration to **0.2 s**.
Jobs 97464/97465 exercised CUDA-aware MPI and direct/equal-rank completed-mesh
partitioning; 97466/97470 exercised host MPI and two-phase saved-mesh import.

Each pair compares omission against an explicit uniform weight: `4.38` for the
nonporous controls and `1` for porous controls. Thus the mixed-material pairs
intentionally have different ownership. Provenance, completion, step/sample counts,
station material/pressure availability, finite eligible fields, timestamps and
fresh PCG residuals were audited. Owner ranks and local element identifiers
may change. The table gives the largest native station-field relative L2
**difference**, not analytical error or speedup; fields include displacement,
velocity, pressure and eligible full-Biot fluid histories.

| Job | Control pair | Largest field difference |
|---|---|---:|
| 97464 | `surface_compact_default` | 1.0161254e-06% |
| 97464 | `full_biot_lumped_new` | 1.4042611e-11% |
| 97464 | `brick_up_lumped_new` | 0% |
| 97464 | `elastic_new` | 1.1316178e-11% |
| 97464 | `acoustic_new` | 1.0759313e-11% |
| 97464 | `acoustic_elastic_new` | 1.1716e-11% |
| 97465 | `buried_compact_default` | 0.012510413% |
| 97465 | `full_varying` | 0.00039041915% |
| 97465 | `fem_open_p` | 1.0706315e-05% |
| 97465 | `vt_open_p` | 3.2272908e-05% |
| 97465 | `full_bkt2_valid` | 2.1203755e-11% |
| 97465 | `full_nonlinear` | 1.2725714e-11% |
| 97466 | `buried_compact_default` | 0.023611129% |
| 97466 | `fem_open_p` | 9.390847e-06% |
| 97470 | `buried_tight` | 4.5425424e-05% |

The initial pair threshold was **0.02% L2**. One observation exceeded it:
near-groundwater pressure at station 1 in host/import job 97466 differed by
**0.0236111%**, with maximum absolute difference **0.1356 Pa** on a reference
peak **893.0431 Pa**. This exceedance is retained in the audit. Job 97470 repeated
that pair with only `poroelastic_mass_tolerance=1e-8`; its largest field difference
fell to **0.0000454254%**, approximately **520 times smaller**. This supports
iterative stopping sensitivity under repartitioning; it is not bitwise invariance
or a general error bound for the unchanged default tolerance `1e-6`.

Extracted production C++ classifier/callback tests with UBSan covered both
formulations, ordinary/acoustic/air cells, dry caps, depth limits, terrain,
origin offsets and explicit weights `1`, `2` and `4.38`. The full clean solver
build, documentation lint/type checks, both documentation production builds
and rendered/parser audits passed. Browser checks confirmed the parameter row,
links and the narrow-screen card layout.

These controls cover standalone elastic, acoustic and both porous formulations,
acoustic–elastic, porous–elastic dry caps and FEM/VT porous–acoustic coupling,
plus BKT2 and nonlinear elastoplastic/plastic solid interiors. They establish
implementation and partition-invariance evidence at the stated tolerances,
not optimal weights for those models. The earlier complete Figure 5 timing
comparison remains the performance evidence for choosing porous weight `1`.
Mixed-model optimal weights, adaptive production performance, reduced-rank
unified handoff and every supported constitutive/interface combination are not
newly calibrated by these short controls.

Investigation helpers, frozen inputs, logs and `analysis.json` are local-only
in the private validation archive; they are
ignored and unavailable in a clean checkout. Existing published benchmark
inputs and historical outputs were not rewritten.
