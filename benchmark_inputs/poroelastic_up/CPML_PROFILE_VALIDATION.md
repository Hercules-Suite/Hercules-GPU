# Poroelastic CPML reflection, kappa and order study — 2026-11-17

The user clarified that “gamma” in this request means `cpml_kappa_max`, **not**
`poroelastic_pressure_stabilization`. Full-Biot pressure stabilization remains
0.0025; this control is inactive for BrickUP.
This follows the [alpha-only comparison](CPML_ALPHA_ONLY_VALIDATION.md), which
found little benefit from lowering alpha and larger errors when increasing it.
Here alpha stays at 0.5 Hz and its floor fraction stays at 0.1 throughout.

Changing kappa and the reflection target has a substantial effect in these
controls. Polynomial order has a smaller, case-dependent effect. The best
individual settings do not necessarily give the best combination. All results
below are measured numerical-reference differences for these finite meshes,
not continuum analytical errors or guaranteed reflection coefficients.
The tables and figures use the v14.0 baseline. To reproduce a curve
labelled “default” with a newer executable, set `cpml_kappa_max=5`
explicitly. The current default of `3` corresponds to the kappa-only
sweep row at reflection target `0.01`, order `2` and alpha `0.5` Hz.

**Scope update, 2026-11-18:** the subsequent 4 Hz Izmir acoustic–elastic run
98047 diverges with `cpml_reflection_target=1e-4`, `cpml_kappa_max=2` and
order `2`, even with 20 layers. Its dominant observed surface growth is in
the land portion of the eastern CPML. The homogeneous porous optima below
must not be promoted to general CPML defaults. See the
[4 Hz investigation](#4-hz-izmir-counterexample-job-98047) below.
The [homogeneous elastic BKT2 control](#homogeneous-elastic-bkt2-control)
compares the same profiles in a flat, fully solid model.

## Controlled experiment

The five families, executable, physical inputs and references are the same as
the preceding alpha-only experiment: BrickUP viscous surface saturation and
buried groundwater, consistent full-Biot viscous surface saturation and buried
groundwater, and consistent full-Biot inviscid surface saturation. The buried
table is 67 m below the flat physical ground. Material constants, source and
station coordinates are recorded in that linked validation record.

Each compact model is a **512 m cube including 128 m side/bottom CPML**, with
no top CPML, uniform 8 m cells, 100-microsecond timestep, and a 25 Hz
central-frequency force source. The short comparison spans 0–0.5 s. The same
four receivers are sampled every two steps. Matching **2048 m cube numerical
references** use the same grid, source, relative station positions and binary.
Their distant-boundary continuum return estimate is later than the comparison
window. The physical top is free and drained. BrickUP retains natural/weak
pressure enforcement and lumped mass; full Biot retains strong enforcement,
consistent physical/CPML mass and 1e-6 PCG tolerance with production optimized
controls. The historical effective partition weight is explicitly fixed at
4.38, including in the references; this study does not calibrate partition cost.

The initial sweep changes **one input at a time** around the default profile:

| Input | Tested values | Baseline |
|---|---|---|
| `cpml_reflection_target` | 0.1, 0.01, 0.001, 1e-4, 1e-6, 1e-8 | 0.01 |
| `cpml_kappa_max` | 1, 1.5, 2, 3, 5, 8, 10 | 5 |
| `cpml_poly_order` | 1, 2, 3, 4 | 2 |

Orders 1–4 are the entire accepted range; kappa must be at least 1 and the
reflection target strictly between zero and one. With a shared fresh baseline,
this gives 15 runs per family, or 75 new runs. Exact parameter comparisons
allow only the selected line to change, and all source/material auxiliary
files are byte-identical within a family. Previously implicit defaults are
written explicitly in all new baselines and variants. The five unchanged
reference runs are reused and their input/binary/runtime provenance is rechecked.

For each field, error is the L2 norm of the unfiltered waveform difference
divided by the reference norm, maximized over receivers. Solid velocity and
Darcy flux use all three vector components. The selection score is the largest
of the velocity-vector, pressure and flux-vector errors, wherever defined.
Dry pressure is undefined and excluded; zero strong surface pressure is checked
in absolute terms and excluded from relative errors. No shifts, gains, filters
or resampling are fitted. A change in the maximizing field/receiver can change
the score without every field improving equally.

## One-parameter results

Every number in the following table is a worst-field error in percent of the
matching reference norm. Each column changes only its named parameter from
the baseline; the best tested value is in parentheses.

| Family | Baseline | Best reflection target | Best kappa | Best order |
|---|---:|---:|---:|---:|
| BrickUP, surface | 10.21738% | 0.76897% (1e-8) | 0.30085% (2) | 5.57503% (4) |
| BrickUP, buried | 11.68406% | 1.14645% (1e-8) | 0.37144% (2) | 5.88762% (4) |
| Full Biot, viscous surface | 1.21231% | 0.09240% (1e-6) | 0.84821% (1) | 1.21231% (2, unchanged) |
| Full Biot, viscous buried | 1.90392% | 0.27224% (1e-8) | 0.77872% (2) | 1.46412% (3) |
| Full Biot, inviscid surface | 12.61246% | 4.43688% (1e-8) | 0.77306% (1) | 6.39670% (4) |

Kappa 1–2 gives large improvements for BrickUP and inviscid full Biot; kappa
8–10 is worse than the baseline in every family. Reducing the reflection
target generally helps, but it is not monotonic in every case: viscous full
Biot at the surface has score 0.09240% at 1e-6 and 0.09899% at 1e-8.
Increasing order helps both BrickUP cases and inviscid full Biot, but order 2
remains best in the viscous surface full-Biot case. Therefore “stronger damping”
or “higher order” is not a reliable standalone tuning rule.

## Combined profiles

After completing the individual sweeps, 25 additional short runs test pairs
and triples of the independently best settings, plus the moderate recipe
$R=10^{-4}$, $m=3$, $\kappa_{\max}=1$ and its best-individual-kappa variant.
Identical profiles are deduplicated. These runs deliberately change several
inputs together; they are separate from the one-parameter comparisons above.
Alpha, its floor fraction, pressure stabilization, timestep and thickness
remain fixed. No waveform-based adjustments are fitted.

The lowest tested score in every family uses **reflection target 1e-4 and
order 3**. The preferred kappa differs:

| Family | Kappa | Solid velocity error | Pressure error | Darcy-flux error | Worst-field score |
|---|---:|---:|---:|---:|---:|
| BrickUP, surface | 2 | 0.07091% | 0.08114% | — | 0.08114% |
| BrickUP, buried | 2 | 0.05300% | 0.08778% | — | 0.08778% |
| Full Biot, viscous surface | 1 | 0.04923% | 0.04871% | 0.05097% | 0.05097% |
| Full Biot, viscous buried | 2 | 0.03544% | 0.08150% | 0.03338% | 0.08150% |
| Full Biot, inviscid surface | 1 | 0.07582% | 0.45728% | 0.17983% | 0.45728% |

The late-window check uses the error norm over 0.25–0.5 s divided by the
**full-record** reference norm. It also improves: inviscid full-Biot pressure
changes from 12.60706% to 0.45645%; BrickUP surface velocity changes from
10.21631% to 0.03496%. Source data, amplitude and timing are identical, and
the histories are not shifted or filtered to obtain these differences.

For example, the following is a case-specific profile for the two surface
full-Biot fixtures, with their existing mesh, material and hydraulic conditions:

```ini
cpml_reflection_target = 1e-4
cpml_kappa_max = 1
cpml_poly_order = 3
cpml_alpha_max_hz = 0.5
cpml_alpha_floor_fraction = 0.1
cpml_thickness_m = 128
```

The best tested BrickUP and buried full-Biot profiles instead use
`cpml_kappa_max = 2`. These are measured candidates for these fixtures, not new
defaults or universally optimal profiles. The search is a finite grid with
selected interactions, not a complete multivariable optimization. Source
bandwidth, mesh spacing, material, hydraulic surface and layer thickness matter.

The interactions are substantial. For BrickUP surface saturation, combining
the individually best $R=10^{-8}$, kappa 2 and order 4 gives a score of
0.23415%, compared with 0.08114% for the moderate combination. For inviscid
full Biot, $R=10^{-8}$, kappa 1 and order 4 gives 0.78657%, compared with
0.45728% for $R=10^{-4}$, kappa 1 and order 3. The viscous full-Biot surface
case prefers order 2 when changing order alone, but its best tested combined
profile uses order 3. Consequently, the gains from individual parameter changes
must not be multiplied or assumed to persist when combined.

Two exact short repeats of each moderate buried-Biot profile check its earlier
run-to-run variation without changing solver tolerance. Across the original
and two repeats, the score is 0.108353–0.108404% for kappa 1 and
0.081485–0.081512% for kappa 2. The measured preference for kappa 2 persists
and exceeds this observed repeat spread. This tests repeatability of that
selection; it does not establish a PCG-tolerance convergence bound or continuum
accuracy.

## Ten-second controls

Each family's best tested short profile was extended to 10 s, changing only
the end time. The same force table reaches zero and stays zero after 0.5 s.
The corresponding default-profile 10 s runs from the alpha study are reused
with their provenance rechecked. All ten long records have 100000 steps and
finite eligible fields; no sustained receiver-envelope growth is observed
over this measured duration. The broad decay is not monotonic in every window.

The following values are **ratios, not percentages**: the maximum peak in
9.5–10 s divided by that receiver/field's peak in 0–0.5 s, over solid velocity
magnitude, absolute pressure and Darcy-flux magnitude where available.
Zero-pressure traces are excluded. The maximum is then taken over receivers
and fields.

| Family | Default profile | Selected profile |
|---|---:|---:|
| BrickUP, surface | 2.61985e-2 | 9.69527e-3 |
| BrickUP, buried | 9.72296e-3 | 3.40975e-3 |
| Full Biot, viscous surface | 1.90970e-11 | 1.96137e-11 |
| Full Biot, viscous buried | 4.60489e-4 | 6.12263e-5 |
| Full Biot, inviscid surface | 3.20374e-4 | 1.22000e-5 |

The selected profiles reduce these late envelopes by approximately 2.7–2.9
times for BrickUP, 7.5 times for buried full Biot, and 26 times for inviscid
full Biot. The viscous surface full-Biot tail remains at a similarly tiny
level with either profile. BrickUP still has appreciable residual oscillations;
better absorption does not imply exact silence. These tail measurements have
no long-duration reference and do not establish stability beyond 10 s, global
energy decay or stability at an acoustic interface.

The first 0.5 s of each long record is also compared with its short counterpart.
The largest relative L2 difference among finite station columns is 1.4956e-4
for the selected buried-Biot profile and 2.5605e-4 for its imported default
versus the fresh baseline. Those discrepancies remain flagged, consistent with
the observed run variation; they are not claimed as exact agreement. All other
prefix comparisons are below 2.41e-6. The separate repeated short comparisons
above retain the same kappa ranking. Per-solve PCG convergence does not guarantee
bitwise reproducibility or bound accumulated waveform differences.

## Meaning of the controls

For normalized distance $r=\xi/L$ into the layer and polynomial order $m$, the
current profile builder uses

```math
d(r)=-\frac{(m+1)c_{\rm ref}\ln R}{2L}\,r^m,
\qquad \kappa(r)=1+(\kappa_{\max}-1)r^m,
\qquad S(s)=\kappa+\frac{d}{s+\alpha}.
```

Decreasing the target $R$ increases the damping profile through $-\ln R$.
The target is a profile-design input, not the measured waveform error.
Increasing $m$ moves more damping toward the outer layer while preserving
the continuum integral of $d$ for fixed $R$, $c_{\rm ref}$ and $L$; it also
changes the spatial kappa profile. Kappa changes the real stretch. In the
high-frequency limit, a plane-wave spatial phase varies with this stretch,
so large kappa can increase spatial resolution demands within the layer.
These relations explain why the controls can interact on a finite grid;
they do not prove that one setting is universally optimal or stable.

## Reproducibility and limits

The binary is unchanged from v13.12, matching all 107 frozen source/build
hashes. Its SHA256 is
`1777c3925b8184fd9ecbfd5d6dd6a5cee116f1a734998bcebbdb3062553dc893`.
Jobs use `ml reset`, checked GSL/PROJ and other library dependencies, CUDA-aware
MPI, one node/four ranks per compact case and `cuda-615-SC`. The imported
references use four nodes/16 ranks and the identical executable. Existing
job reservations are not changed.

The one-parameter job batches are 97714–97716 (BrickUP surface), 97717–97719
(BrickUP buried), 97720–97722 (Biot viscous surface), 97723–97725 (Biot viscous
buried), and 97726–97728 (Biot inviscid surface). Within each triplet, batches
cover reflection target including the baseline, kappa, and order. The references
remain jobs 97598–97602 in that family order. All 75 initial runs complete with
finite eligible fields, correct sample counts and verified frozen inputs.
Every reported full-Biot true residual satisfies the fixed 1e-6 tolerance.

Combined profiles run in jobs 97731–97735 in the same family order, with
6, 6, 2, 6 and 5 cases respectively. All 25 complete and pass the same audits.
Buried-Biot paired repeats run in jobs 97746 and 97747 and pass those audits
with unchanged parameter/source hashes. Each job runs kappa 1 followed by
kappa 2 for the moderate combined profile.

Selected-profile long controls are jobs 97740–97744 in the same family order.
Imported default long controls are 97617, 97619, 97621, 97623 and 97625.
All **109 new executions in 27 allocations** finish with Slurm exit `0:0`:
75 one-parameter cases, 25 combinations, four repeats and five long controls.
Together with five imported short references and five imported long defaults,
119 records pass the binary/input/runtime, field/sample and reported PCG audits.
The separate prefix repeatability observations are retained as described above.

The [preceding repeatability check](CPML_ALPHA_ONLY_VALIDATION.md#repeatability-check)
found small run-to-run differences for buried full Biot. Displayed digits are
provenance, not a universal precision guarantee. The common discretization,
stabilization, source and hydraulic-boundary errors are not removed by these
numerical comparisons. This study does not qualify other bandwidths, layer
thicknesses, meshes, topography, acoustic interfaces, heterogeneous/nonlinear/
viscoelastic interiors or purely elastic/acoustic CPML. Investigation helpers,
fixtures, raw histories and generated figures remain local-only and are not
available in a clean checkout.

## Complete single-parameter sweep

All entries below are the same worst-field relative L2 score in percent.
Only the named input changes from the default profile.

### Reflection target

| Value | BrickUP surface | BrickUP buried | Biot viscous surface | Biot viscous buried | Biot inviscid surface |
|---:|---:|---:|---:|---:|---:|
| 1e-08 | 0.76897% | 1.14645% | 0.09899% | 0.27224% | 4.43688% |
| 1e-06 | 1.31523% | 1.95144% | 0.09240% | 0.28509% | 5.62313% |
| 0.0001 | 3.06949% | 4.11734% | 0.18095% | 0.51774% | 7.94715% |
| 0.001 | 5.30304% | 6.61750% | 0.37772% | 0.90326% | 9.90539% |
| 0.01 | 10.21738% | 11.68406% | 1.21231% | 1.90392% | 12.61246% |
| 0.1 | 22.96026% | 24.55055% | 6.57190% | 7.03885% | 17.26300% |

### Kappa maximum

| Value | BrickUP surface | BrickUP buried | Biot viscous surface | Biot viscous buried | Biot inviscid surface |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.44627% | 0.58424% | 0.84821% | 0.93082% | 0.77306% |
| 1.5 | 0.30738% | 0.47224% | 0.87234% | 0.84398% | 0.89514% |
| 2 | 0.30085% | 0.37144% | 0.88150% | 0.77872% | 1.60744% |
| 3 | 0.92923% | 1.31077% | 0.92444% | 0.83305% | 5.01251% |
| 5 | 10.21738% | 11.68406% | 1.21231% | 1.90392% | 12.61246% |
| 8 | 30.84904% | 34.53003% | 5.56504% | 8.98977% | 20.85705% |
| 10 | 41.10304% | 51.75435% | 10.09700% | 17.14724% | 30.03359% |

### Polynomial order

| Value | BrickUP surface | BrickUP buried | Biot viscous surface | Biot viscous buried | Biot inviscid surface |
|---:|---:|---:|---:|---:|---:|
| 1 | 27.49331% | 31.87990% | 1.78309% | 5.02512% | 20.45590% |
| 2 | 10.21738% | 11.68406% | 1.21231% | 1.90392% | 12.61246% |
| 3 | 7.09380% | 7.68910% | 1.27361% | 1.46412% | 8.88648% |
| 4 | 5.57503% | 5.88762% | 1.41694% | 1.55632% | 6.39670% |


## Complete combination results

Alpha is 0.5 Hz in every row. Scores retain the same definition as above.

| Family | Reflection target | Kappa | Order | Worst-field score |
|---|---:|---:|---:|---:|
| brick_surface | 1e-08 | 2 | 2 | 0.32796% |
| brick_surface | 0.01 | 2 | 4 | 0.30029% |
| brick_surface | 1e-08 | 5 | 4 | 0.19651% |
| brick_surface | 1e-08 | 2 | 4 | 0.23415% |
| brick_surface | 0.0001 | 1 | 3 | 0.08527% |
| brick_surface | 0.0001 | 2 | 3 | 0.08114% |
| brick_buried | 1e-08 | 2 | 2 | 0.27257% |
| brick_buried | 0.01 | 2 | 4 | 0.42498% |
| brick_buried | 1e-08 | 5 | 4 | 0.23054% |
| brick_buried | 1e-08 | 2 | 4 | 0.23416% |
| brick_buried | 0.0001 | 1 | 3 | 0.09523% |
| brick_buried | 0.0001 | 2 | 3 | 0.08778% |
| biot_surface | 1e-06 | 1 | 2 | 0.09520% |
| biot_surface | 0.0001 | 1 | 3 | 0.05097% |
| biot_buried | 1e-08 | 2 | 2 | 0.21859% |
| biot_buried | 0.01 | 2 | 3 | 0.84062% |
| biot_buried | 1e-08 | 5 | 3 | 0.16240% |
| biot_buried | 1e-08 | 2 | 3 | 0.19377% |
| biot_buried | 0.0001 | 1 | 3 | 0.10838% |
| biot_buried | 0.0001 | 2 | 3 | 0.08150% |
| biot_inviscid | 1e-08 | 1 | 2 | 0.70206% |
| biot_inviscid | 0.01 | 1 | 4 | 0.76918% |
| biot_inviscid | 1e-08 | 5 | 4 | 2.39576% |
| biot_inviscid | 1e-08 | 1 | 4 | 0.78657% |
| biot_inviscid | 0.0001 | 1 | 3 | 0.45728% |

## Separate Izmir acoustic–elastic stability check

The 32-node test used the 1 Hz Izmir input with
`cpml_reflection_target=1e-4`, `cpml_kappa_max=2`, and
`cpml_poly_order=2`. This is a separate follow-up to the homogeneous porous
study above; it is not included in its run counts or error tables.

Horizon job **97896** completed the nominal **100 s** record on **32 nodes /
128 GPU solver ranks**, using `cuda-615-SC`, with Slurm `COMPLETED`, exit
`0:0`, and **11 min 47 s** allocation elapsed time. No divergence was observed
in the saved fields over this duration. This is an observed stability result,
not an absorption-error or mesh-convergence measurement.

The original input already specified reflection target `1e-4` and omitted
polynomial order, whose default is `2`. The isolated input changes only kappa
from `1` to `2` and adds the explicit order `2` line. Source, material and
geometry files were copied and hashed; the original run directory was not
modified. The supported unified launcher regenerates the mesh and hands it
from 4,608 CPU mesh ranks to 128 GPU solver ranks.

Relevant retained settings and effective values:

| Control | Tested value |
|---|---|
| Physics | Acoustic–elastic, BKT2 attenuation, FEM topography; poroelasticity off |
| Acoustic interface | `field-specific-constraints`, default `point` projection, all face orientations |
| Water free surface | Linearized gravity enabled |
| Mesh frequency / points per wavelength / minimum Vs | 1 Hz / 10 / 120 m/s |
| Source | Unchanged discrete SRFH source; 1 Hz is the mesh frequency, not a Ricker central frequency |
| Maximum octant level | 13 |
| Surface refinement | Omitted in input; current default `maximum_level` |
| Timestep | 0.0015 s; 66,666 steps |
| Reported critical timestep | 0.002706137 s, including CPML CFL factor 0.9 |
| CPML layers / derived thickness | 15 / 5,273.4375 m, from 351.5625 m boundary cells |
| CPML alpha maximum / floor fraction | Automatic 0.5 Hz / 0.1 |
| CPML partition weight | Explicit 2 for all CPML cells |
| Acoustic bulk sponge | 0 |
| CPML material extension | Enabled, default `elevation` reference |
| Outer CPML clamp | Enabled |

The mesh has 518,778,463 elements, including air, and 590,920,401 nodes.
The solver reports 102,813,727 elastic and 14,708,467 acoustic CPML elements;
automatic CPML reference Vp is 8,080 m/s. The solver's wall timer is 614.14 s,
excluding part of the allocation's MPI startup time. Solver stderr is empty.

Every saved value was checked: all 19 stations have 3,334 samples through
99.99 s, and all seven planes have 667 complete snapshots through 99.9 s.
There are no non-finite station displacement, velocity or acceleration
values, and no non-finite plane displacements. Maximum absolute components
over all stations are 0.2121863 m, 0.1446047 m/s and 0.6895292 m/s².
The corresponding velocity and acceleration maxima during 90–100 s are
0.01471973 m/s and 0.07039726 m/s², below their earlier maxima. The whole
surface plane, which includes reconstructed acoustic displacement and CPML,
peaks at 3.52330 m at 52.5 s and ends at a spatial maximum of 1.85119 m.
These amplitude histories do not show runaway late growth through 100 s.

The test uses the same validated v13.12 executable as the profile study,
SHA-256 `1777c3925b8184fd9ecbfd5d6dd6a5cee116f1a734998bcebbdb3062553dc893`,
with `ml reset` and CUDA-aware MPI. It does **not** reproduce the older binary
referenced by the original Izmir launcher. Older Izmir logs also show
`wavelength` surface refinement, whereas this unchanged input selects the
current `maximum_level` default. Therefore this result does not rule out a
failure with that older executable/mesh or a different alpha, timestep or
CPML thickness. There is no matched kappa-1 control, remote-boundary reference,
pressure output, all-DOF audit or test beyond 100 s in this follow-up.

Inputs, manifests, logs, numerical audits and plots are retained in the private validation archive. Production solver code and
defaults are unchanged.

## 4 Hz Izmir counterexample: job 98047

### Evidence and provenance

The maintainer submitted jobs 97979 and 98047. Both used 32 Horizon nodes,
128 GPU solver ranks, `cuda-615-SC`, the updated unified MPI core mapping
and binding, and `parameters_4Hz_notaper.in`. This investigation submitted
no additional solver jobs and changed no solver code or production input.

Both Slurm jobs report `COMPLETED`, exit `0:0`; this is not a numerical
stability test. Job 98047 ran for 1 h 45 min 41 s. Its 19 station histories
and seven planes have complete, finite records, yet develop catastrophic
late growth. Each station has 12,500 samples through 99.992 s; each plane
has 2,500 snapshots through 99.96 s. All saved displacement components,
station velocities and accelerations were scanned, not just final values.

The input/log differences include 15 layers in 97979 versus 20 in 98047:

| Quantity | Job 97979 | Job 98047 |
|---|---:|---:|
| Mesh frequency | 4 Hz | 4 Hz |
| Timestep | 0.4 ms | 0.4 ms |
| Reflection target / kappa / polynomial order | `1e-4` / `2` / `2` | `1e-4` / `2` / `2` |
| Alpha maximum / floor fraction | 0.5 Hz / 0.1 | 0.5 Hz / 0.1 |
| CPML layers | 15 | 20 |
| Derived physical thickness | 1,318.359375 m | 1,757.8125 m |
| Reported critical timestep, including CPML CFL factor 0.9 | 1.561161 ms | 1.561161 ms |
| Acoustic–elastic fragments with centres inside CPML | 394,932 | 532,065 |
| Gravity free-surface acoustic CPML elements | 368,416 | 492,726 |

The original 97979 wavefields were not found in a separately preserved
directory; both jobs wrote the same output paths. The inspected station
and plane modification times match the end of 98047. Only 98047's growth
is quantified here; the logs alone do not establish 97979's stability or
the measured change in growth rate from increasing the layer count.

Job 98047 retains FEM topography, `wavelength` surface refinement, maximum
octant level 15, BKT2 interior attenuation, gravity water free surface,
`field-specific-constraints`, `point` interface projection, and coupling
on all face orientations. Poroelasticity is off. CPML uses full Hex8
operators, no additional material damping or acoustic bulk sponge, no
M-PML cross mixing, a clamped outer boundary, automatic reference Vp
8,080 m/s, and per-element-centre profile sampling. Material extension
uses the `elevation` reference. The mesh has 4,391,499,448 elements and
4,422,967,867 nodes, including 338,421,234 elastic and 8,584,033 acoustic
CPML elements.

### Measured growth and location

Maximum absolute station displacement components are about 0.206 m during
10–20 s and 0.0764 m during 50–60 s. They subsequently reach 26.0 m during
70–80 s, 12,023.6 m during 80–90 s, and 2,823,235 m by the last sample.
The late log-envelope slope over 80–100 s is approximately 0.578/s,
equivalent to a 1.20 s doubling time. This is clear numerical divergence
despite the absence of NaNs, infinities or a nonzero process exit.

The saved whole-surface plane places the earliest clearly unphysical growth
in the eastern CPML, near solver X = 41–42.25 km, Y = 89.5–89.75 km. X is
north, Y east; the domain spans 0–90 km in both directions and the eastern
CPML begins at Y = 88.2421875 km. The surface samples are 250 m apart.

| First saved surface amplitude exceeding threshold | Land inside CPML | Anywhere outside CPML | Over water |
|---|---:|---:|---:|
| 10 m | 43.12 s | 69.00 s | 76.48 s |
| 100 m | 49.36 s | 73.60 s | 80.56 s |
| 1,000 m | 53.36 s | 77.64 s | 84.52 s |
| 1,000,000 m | 65.88 s | 90.16 s | 96.48 s |

At 60 s the land-CPML maximum is approximately 24,863 m, while the maximum
over water is 6.66 m and that outside CPML is 3.05 m. The overall surface
maximum reaches approximately 1.52e14 m at 99.76 s.

Water/land classification was obtained by bilinearly sampling the original
30 m topography and seabed grids at the plane points, using the same array
ordering and interpolation convention as the solver. The dominant location
at X = 41 km, Y = 89.75 km has equal surface and seabed heights, about
29.85 m; the nearest sampled water point is approximately 24.8 km away.
This is a physical geometry classification, not a retained element-owner
dump. The surface histories localize the dominant observed growth; they do
not identify the first unstable internal DOF or prove its generating operator.

### Interpretation and next discriminating checks

The evidence supports a CPML-related instability. It does **not** establish
acoustic–elastic coupling or the shoreline as the primary cause: the dominant
observed growth is on land, followed later by large amplitudes in the interior
and water-covered surface. Elastic CPML, its free-surface discretization,
material continuation and adaptive-mesh profile sampling warrant direct
investigation. Interface coupling could still contribute indirectly or at
unsampled internal locations.

The interface stretch is not simply missing. The run explicitly reports
the `J/S_j` cofactor on both sides of all 532,065 CPML interface fragments,
with common fragment-centre coefficients, six CFS histories per quadrature
point per side, and one quadrature point per fragment. Gravity surface
corrections are also active. Their presence does not prove stability of the
combined discrete operator. Published PML analysis likewise requires careful
treatment of auxiliary equations, interfaces and boundaries; it is not
enough to establish bulk absorption in a homogeneous test. See
[Duru, Gabriel and Kreiss](https://arxiv.org/abs/1802.06388).
Instability in some heterogeneous 3-D elastic configurations also predates
acoustic coupling; see the stability section of
[Marcinkovich and Olsen](https://doi.org/10.1029/2002JB002235).
These papers concern other discretizations and do not diagnose this kernel.

For the implemented polynomial profile,

```math
d_{\max}=-\frac{(m+1)c_{\rm ref}\ln R}{2L}.
```

At fixed thickness and reference velocity, changing `R=0.01` to `1e-4`
doubles `d_max`. Changing kappa from `5` to `2` at the same time makes the
outer-edge `d_max/kappa_max` coefficient five times as large. For 98047's
thickness, the values are approximately 31.75/s for the tested combination
versus 6.35/s for the default profile. This is an inverse-stretch coefficient,
not a measured attenuation rate or proof of the instability mechanism.
Stronger nominal absorption is not a stability guarantee.

The 0.4 ms timestep is below the reported 1.561161 ms estimate and below
the separate printed damping-related minima (about 1.001 and 0.986 ms).
There is no obvious violation of those reported limits, but they are not a
spectral stability proof for the complete CPML/interface/free-surface update.
A smaller-timestep control remains necessary to distinguish mechanisms.

The baseline at the time of run 98047 was `R=0.01`, `kappa_max=5`, order `2`
and automatic alpha 0.5 Hz. The porous profile candidates remain specific to
their validated homogeneous cases. Useful next checks, **not run here**, are:

1. Repeat the exact 98047 configuration with the default R/kappa profile,
   retaining its 20-layer physical thickness, timestep, geometry and binary;
   then vary R and kappa separately if the default control is stable.
2. Halve the timestep with the failing profile and compare growth versus
   physical time, including the eastern-CPML surface points.
3. Reproduce the eastern land segment in a smaller purely elastic model,
   retaining its material and topographic profiles. This tests whether
   fluid coupling is necessary before changing the coupled interface.

The successful 1 Hz run 97896 is not a matched frequency-only control:
it had a 5.273 km CPML, different geometry/material/source inputs, a different
surface-refinement policy, and a separately pinned executable. No conclusion
about a universal stable frequency or layer count follows from these runs.

Logs, copied station records, fixture-specific read-only analysis helpers,
sampled geometry/wavefield maps and figures are retained locally. Compact remote summaries
are in the matching Horizon validation analysis directory; the original
production outputs are unchanged. Helpers and raw fixtures are local-only
and unavailable in a clean checkout.

## Homogeneous elastic BKT2 control

On 2026-11-18, the default CPML profile gave a **6.586%** worst-receiver
velocity error in the homogeneous BKT2 solid. The two order-3 profiles
selected in the porous study gave **1.601%** and **1.608%**. All eight
receivers improved with either profile. The order-2 profile from Izmir job
98047 gave **1.563%** in this model.

### Model and numerical reference

The 512 m cube has a flat, traction-free top, uniform 8 m elements and a
128 m CPML on the four sides and bottom (16 elements). Elastic material
properties are uniform: Vp = 3296.27 m/s, Vs = 1961.16 m/s and density =
1664 kg/m³. The run uses the former porous benchmark's drained-frame
velocities, with poroelasticity disabled. The 25 Hz central-frequency
Ricker force is applied at 256 m depth. The timestep is 100 µs. Eight
receivers are 80 m north of the source at depths of 0, 0.5, 8, 32, 59,
67.5, 75 and 99 m.

The physical solid uses `type_of_damping=bkt2` and synthetic
`Qp=Qs=50`, with `use_parametric_q_factor=no` and `use_infinite_qk=no`.
The solver skips BKT material histories in CPML elements. The comparison
therefore measures the combined effect of the finite layer, the outer
boundary, and the change in the attenuation treatment at the CPML entrance.
Every compact run used the same layer thickness, alpha floor fraction of
0.1 and automatic reference Vp of 3296.27 m/s. The four primary profiles use
the automatic alpha maximum of 0.5 Hz. The default case omitted the three
profile inputs and the log confirmed their resolved values. All compact
cases used a zero multiaxial cross fraction and a zero-displacement clamp
on the outer CPML faces. The mesh-frequency input was 24.2718 Hz and the
force time history was copied directly from the porous benchmark.

The numerical reference is a 2048 m cube filled with the same BKT2 solid,
with the same free top, mesh spacing, timestep, source history and relative
receiver positions. Its lateral and bottom boundaries use dashpots. The
shortest lateral return estimated from the input Vp is at least 0.597 s.
Waveform errors use the first 0.5 s, sampled every 0.2 ms. For each receiver
the error is the L2 norm of the three-component difference divided by the
reference norm, using the original amplitudes and times. The table gives
the maximum over receivers for each field independently. The late column
uses only 0.25–0.5 s in the numerator and the full record in the denominator.

| Profile | Reflection target | Kappa max | Order | Velocity error | Displacement error | Late velocity error |
|---|---:|---:|---:|---:|---:|---:|
| v14.0 default | 0.01 | 5 | 2 | 6.586% | 4.196% | 6.460% |
| Porous candidate, kappa 1 | 1e-4 | 1 | 3 | 1.601% | 1.664% | 0.911% |
| Porous candidate, kappa 2 | 1e-4 | 2 | 3 | 1.608% | 1.670% | 0.916% |
| Izmir job 98047 profile | 1e-4 | 2 | 2 | 1.563% | 1.595% | 0.897% |

The default's largest velocity error occurred at 75 m depth. All three
alternative profiles had their largest velocity error at 8 m. Both porous
profiles reduced the worst-receiver error by about 76%. Most of the
default's discrepancy was visible after 0.25 s.

### Changing one parameter at a time

Each pair of columns below is a separate sweep around the default
`R=0.01`, `kappa_max=5`, order `2`, alpha `0.5` Hz. Each error is the
worst-receiver velocity error over 0–0.5 s. The rows across different
sweeps do not represent combined inputs.

| R | Error | Kappa max | Error | Order | Error | Alpha max (Hz) | Error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 17.288% | 1 | 1.549% | 1 | 16.605% | 0.0078125 | 6.568% |
| 0.01 | 6.586% | 1.5 | 1.710% | 2 | 6.586% | 0.03125 | 6.569% |
| 0.001 | 3.412% | 2 | 1.607% | 3 | 4.656% | 0.125 | 6.572% |
| 1e-4 | 2.186% | 3 | 1.553% | 4 | 3.814% | 0.5 | 6.586% |
| 1e-6 | 1.561% | 5 | 6.586% | | | 1 | 6.605% |
| 1e-8 | 1.529% | 8 | 24.031% | | | 2 | 6.646% |
| | | 10 | 37.724% | | | 6.25 | 6.877% |
| | | | | | | 12.5 | 7.386% |
| | | | | | | 25 | 9.071% |

Changing only the reflection target to `1e-8` gave the lowest measured
velocity error, 1.529%. With the default target, changing only kappa to
`1` gave 1.549%. Lowering alpha to 0.0078125 Hz improved velocity error
by only 0.017 percentage points. Raising alpha to the 25 Hz source
frequency increased the error to 9.071%.

### Long records and run provenance

All four primary profiles produced finite station displacement and velocity
through 10 s. For each receiver, the peak speed in a half-second window was
divided by that receiver's first-half-second peak speed. The largest final
window ratio over the eight receivers was 3.02e-7 for the default, 7.53e-10
for the kappa-1 porous candidate, 3.39e-9 for the kappa-2 porous candidate
and 2.32e-9 for the Izmir profile.

The default and Izmir profiles were each extended to **100 s** (one million
steps). Only the end time and station output interval changed. The 100 s
runs recorded every 1 ms. Over the common 10 s interval their velocity
traces agreed with the earlier runs sampled at those times to within
2.8e-11 relative L2 at every receiver. All displacement and velocity
samples remained finite.

| Profile | Largest peak-speed ratio anywhere in 40–100 s | Peak-speed ratio in 99.5–100 s |
|---|---:|---:|
| Default | 1.09e-9 | 2.49e-13 |
| Izmir job 98047 profile | 2.14e-13 | 4.20e-16 |

Neither profile showed sustained late growth at the receivers. In particular,
the profile that diverged in the 4 Hz Izmir model decayed through 100 s
in this homogeneous control. The complex model also changes topography,
material heterogeneity, mesh, source bandwidth, physical thickness, timestep
and acoustic coupling. The comparison cannot isolate which of those differences
produced the instability in job 98047.

All runs used the same v13.12 solver binary (SHA256
`1777c3925b8184fd9ecbfd5d6dd6a5cee116f1a734998bcebbdb3062553dc893`)
on Horizon's `cuda-615-SC` reservation with CUDA-aware MPI. The reference
used four nodes and 16 GPU ranks. Each compact model used one node and
four GPU ranks.

| Horizon job | Purpose |
|---|---|
| 98402 | 2048 m numerical reference |
| 98403, 98404, 98405, 98406 | 10 s default, porous kappa 1, porous kappa 2 and Izmir profiles, respectively |
| 98407, 98408, 98409 | Reflection-target, kappa and polynomial-order sweeps, respectively |
| 98412, 98413 | 100 s default and Izmir profiles, respectively |
| 98415, 98416, 98417, 98419 | Alpha-only sweep |

Inputs, native station records, analysis scripts, receiver CSV and figures
are retained privately.
