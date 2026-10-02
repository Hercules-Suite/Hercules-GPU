# Izmir gravity-surface investigation — 1 October 2026

This investigation separates a confirmed surface-CPML discretization error
from surface-gradient recovery error and sensitivity to boundary settings.
It does not establish that the remaining long-period water motion is physical.
The source is unchanged and unfiltered in every case.

## Experiment and provenance

Horizon job 40658 ran four cases sequentially on 128 nodes; job 40685 ran a
separate nodal-memory candidate, which failed validation. Both used 18,432
meshing ranks and 512 GPU solver ranks, MPI-IO plane output, and `IO_PES=0`. The original user's repository,
source directory, and completed output directories were preserved.

| Setting | Value |
| --- | --- |
| `simulation_wave_max_freq_hz` | 1 Hz |
| `maximum_octant_level` | 12 |
| `simulation_delta_time_sec` | 0.0005 s |
| Duration | 100 s |
| Source | Original SRFH, 160 point sources, `source_is_filtered=0` |
| Acoustic formulation | `velocity_potential` |
| Acoustic numerical regularizer | `acoustic_quality_factor=100` |
| Plane sampling | 361 × 361 surface grid, 250 m spacing, 1,000 frames at 0.1 s |
| Station sampling | 200,000 rows at every solver step |

The four baseline-binary cases change one setting at a time:
`clamp_yes` uses gravity, a 439.453125 m CPML, and the outer clamp;
`clamp_no` disables that clamp; `thick_cpml` doubles thickness to 878.90625 m;
`gravity_no` selects the pressure-release surface instead of gravity.
`nodal_cpml` uses the experimental solver with byte-identical `clamp_yes` parameters.

Baseline executable source: `8b048243e644662e767ec3714a155c9903d22d24`.
Baseline binary SHA-256:
`36aa645de7fc9092deac6bfc7889880bf4d950426c23d4debbc2cb14bf54a298`.
Rejected candidate binary SHA-256:
`5930c7d0bac015e7f45fc57f145b8b75ecd8da8e68d0e56f47c25ca8bab1f907`.
The candidate build uses that baseline plus the recorded `nodal_cpml.patch`;
local intervening solver edits only change comments. Diagnostic tools are from
`83b0c26`; the candidate-only CFS regression is retained in the evidence
bundle. Builds used NVHPC 26.5, CUDA 13.2, Open MPI 5.0.11 and `sm_100` on Horizon.

Remote evidence is under
`/scratch/08458/wz_tacc/hercules/Izmir/gravity_diagnostic_1hz_20261001`.
`source_sha256.json` verifies all original source `.in` files against the run
copy/references. `provenance.json` records job IDs, input hashes, executable
identities and timestamps. Compact station histories, peak reports, logs,
inputs, patch and build records are retained there and in the local report
bundle. Raw output remains on Horizon.

## What coarsening changed

The reduced mesh has 195,589,222 elements versus 1,112,763,417 previously.
Retained acoustic probes use 21.97265625 m elements, previously 10.986328125 m.
Frequency-dependent settings also change: automatic CPML alpha becomes
pi rad/s, and the acoustic stiffness regularizer changes with the frequency.
This is a reduced diagnostic configuration, not a clean mesh-convergence pair.
Keeping the original source does not remove frequencies above 1 Hz.

On the same 130,321-point surface observation grid, 2,741 formerly acoustic
points become solid, 1,985 inactive points become solid, and 91 solid points
become inactive. Stations 24–27 at `(100,180)` and `(83,270)` change from water
to elastic elements and are excluded from acoustic comparisons. Mixed
water/solid elements at the refinement limit are kept elastic by the current
material rule. The remaining ten acoustic probes retain Vp=1500 m/s and
rho=1025 kg/m³. In the thicker-CPML case stations 18–21 are inside CPML;
case-specific membership is recorded in the probe manifest.

## Surface operator and reconstruction

The existing surface mass and damping are lumped at four face nodes, but the
current CFS restoring/memory kernel uses their face average. It therefore
applies damping to a zero-mean face mode without its paired restoring/memory
force. An experimental correction kept four histories per tangential axis
and face,
using node-major GPU storage, without adding MPI exchanges or host/device
transfers. It was rejected after the full Izmir test and removed from the
production branch. The isolated Horizon build is marked `REJECTED_CANDIDATE.txt`.

The candidate-kernel host regression checks all four bilinear face modes,
two tangential axes, nontrivial face mapping, and zero damping against
closed-form CFS step responses. The shifted-surface GPU fixture passes with
one and four ranks and both MPI backends. Maximum plane differences are
2.04e-17 m between rank counts and 2.42e-17 m between backends; the old/new
pressure-release control differs by at most 8.39e-19 m. These are local
consistency and implementation checks, not a CPML reflection benchmark.

A separate analytical tool extracts production Hex8 stiffness and lumped
volume/surface masses. For flat, lossless water of depth 175.78125 m and
250 m wavelength, h=21.97265625 m gives −1.253% frequency error but 19.65%
error in surface-gradient motion amplitude. Halving h gives −0.316% and
11.58%. A two-element water column can have very accurate frequency while
its recovered surface gradient remains inaccurate. A static bottom-motion
check also finds a spurious recovered vertical offset that converges with h.
See [the analytical method](Extended%20README.md#check-gravity-wave-dispersion-and-reconstruction).

The production station and plane time integrators pass an independent
potential-gradient identity check, including a 100 s slow ramp. With zero
initial state, they reproduce the time integral of the discrete potential
gradient. This check found no additional drift caused by output time
integration; it does not validate spatial gradient recovery or the evolved
potential field.

The candidate's first checked nonfinite station sample occurs at 4.8325 s;
station 18 follows at 4.833 s. The first nonfinite plane frame is at 4.9 s;
the last fully finite frame at 4.8 s peaks at 0.0015775 m. These surface
observations do not identify the first failing internal state. This is why
the successful small tests do not justify adopting the candidate or claiming
that it fixes the original slow motion. Its patch and finite-prefix histories
are retained only as rejected experimental evidence.

## Reduced-run observations

| Case | Surface station 18 uX (m) | uZ (m) | Pressure-derived surface uZ (m) | Station 19 uX at 5 m depth (m) | Station 19 late mean vX (m/s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Clamp on | 6.875537 | 1.312784 | 0.117471 | 7.292068 | 0.12176873 |
| Clamp off | 5.743611 | 1.013036 | 0.044942 | 6.096241 | 0.069688202 |
| Gravity off | 0.000000 | 0.080242 | Not applicable | -0.003636 | 2.7043466e-06 |

Station 18 is at `(2,93)`: X=500 m north, Y=23,250 m east, on the water
surface. Station 19 is 5 m below it. All station values above are raw,
unfiltered values at 99.9995 s; late velocity is averaged over 80–100 s.
The pressure-based surface value uses positive-down `-p/(rho*g)` outside
CPML. Plane peaks use the separate 0–99.9 s saved time grid.

The same-mesh metadata audit found no medium, material, element-size or CPML
membership changes among clamp on/off, gravity off, and the nodal candidate.
The doubled layer changes 2,194 surface points' material/element properties
and 2,840 CPML memberships (including 16 acoustic-to-solid and 59
solid-to-acoustic transitions); common interior materials are unchanged.

Turning the clamp off changes solid plane histories by up to
0.287213 m outside CPML and 0.100998 m at points at least 20 km from a lateral boundary. These are maximum pointwise
component differences over all saved times, not differences between global
peaks. Thus the boundary setting affects the interior solid response as well
as artificial CPML samples. The gravity-off difference at the same distant
solid set is at most 0.001114 m.
`case_comparisons.json` retains per-axis values, locations, sample times and
common-material selection counts.

The doubled-CPML case is invalid: station 18 first becomes nonfinite at
7.861 s and plane 0 first contains nonfinite values at 7.9 s. The last fully
finite plane frame at 7.8 s has a 0.32275 m peak. The saved surface histories
do not reveal the first failing internal field or element, so these data
cannot identify its cause. The normal solver exit does not make its output
valid. No peak or reflection conclusion is taken from its nonfinite output.
The `thick_cpml_failure` report retains finite prefixes and startup metadata.

## Interpretation and limits

The outer clamp affects the solution, but turning it off leaves metre-scale
horizontal motion and a large surface pressure/motion discrepancy. Disabling
gravity greatly reduces horizontal displacement at the submerged station 19,
which avoids the trivial zero horizontal surface gradient imposed by the
pressure-release condition. This associates the slow response with the
gravity-enabled coupled system; it does not establish a unique defect.

Spatial recovery error is demonstrated independently, whereas extra output
integration drift was not found in the synthetic identity tests. The proposed
surface-CPML correction failed the field case and was rejected. Thus this
investigation identifies separate numerical limitations but does not supply
a validated fix for the remaining slow Izmir motion.

The next solver experiment should locate the first nonfinite internal field
and element in short reproductions of the rejected cases, before changing
more CPML terms. A separate manufactured gravity-wave test of surface flux
and element-gradient recovery is needed before replacing the reported surface
motion. The two investigations should preserve raw histories and vary one
operator or resolution control at a time.

The 1 Hz/level-12 set is useful for controlled debugging, but changed shoreline
classification and source content beyond its target mesh frequency limit its
use as a production result. A gravity-off comparison changes the physical
free-surface condition, so smaller motions in that case are not by themselves
a reason to disable gravity. Similarly, changing the outer clamp or CPML
thickness does not establish improved absorption merely by reducing a peak.

Local validation of the retained branch: 41 Python tests, documentation
lint/type checking, and 40 rendered-documentation/parser-audit tests passed. The experimental solver
compiled and GPU/MPI fixtures ran on Horizon, but it failed the field-case
check. Documentation tests are not solver validation.
