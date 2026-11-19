# Poroelastic CPML: alpha-only comparison, 2026-11-17

Changing only `cpml_alpha_max_hz` does not recover the large improvement seen
with the earlier modified CPML profile. In these five matched controls,
increasing alpha above the default 0.5 Hz increases the measured boundary
differences. Lowering it approaches a plateau: the best tested value,
0.0078125 Hz, reduces the worst-field error by only 0.04–0.25% **relative to
the original error**, not by that many percentage points. The production
default remains 0.5 Hz. These results do not establish a universal optimal alpha.

The [subsequent reflection/kappa/order study](CPML_PROFILE_VALIDATION.md) keeps
alpha fixed and measures those controls separately, then tests combinations.

## Controlled inputs

The preceding [consistent-CPML comparisons](WOLF_2022_VALIDATION.md#consistent-full-biot-cpml-extension-2026-11-17)
changed five profile controls together. This follow-up holds those other inputs
fixed and changes only `cpml_alpha_max_hz` within each family:

`0.0078125`, `0.03125`, `0.125`, `0.5`, `1`, `2`, `6.25`, `12.5`, `25` Hz.

The original seven-value sweep started at 0.125 Hz; the two lower values were
added because its minimum was at the lower endpoint. A nonpositive alpha input
selects the 0.5 Hz fallback and therefore is not a zero-alpha experiment.
Internally, alpha is converted to inverse seconds by multiplying by $2\pi$.
The alpha-floor **fraction** remains 0.1; its absolute value consequently scales
with alpha. No other input parameter or auxiliary source/material file changes
between the alpha variants of a family.

| Control | Fixed value |
|---|---|
| Total compact model | 512 m cube, including side/bottom CPML |
| CPML thickness | 128 m (16 cells); no top CPML |
| Physical interior | 256 m by 256 m horizontally, 384 m deep |
| Mesh and timestep | Uniform 8 m; 100 microseconds |
| Short comparison | 0–0.5 s; 5000 steps; output every two steps |
| Source | Discrete 25 Hz central-frequency Ricker force at (256, 256, 256) m |
| Force amplitudes | (1e9, 4e8, -2e8) N in solver X/Y/Z |
| Reflection target / polynomial order | 0.01 / 2 |
| Kappa maximum / alpha-floor fraction | 5 / 0.1 |
| Multiaxial cross fraction | 0 |
| Outer termination | Displacement clamp; regular Cartesian CPML |
| Partition weight | Explicit 4.38, preserving the earlier effective value |
| Full-Biot mass | Consistent in both physical and CPML elements |
| Full-Biot stabilization / PCG tolerance | 0.0025 / 1e-6, with production optimized solver controls |
| Surface drainage | BrickUP natural/weak, penalty 1; full Biot strong |

Pinning weight 4.38 prevents the recent material-dependent default of 1 for
porous CPML from changing this comparison. This is an accuracy experiment,
not a new partition-weight calibration. The mesh input
`simulation_wave_max_freq_hz=24.271798903240597` and 10 nodes per wavelength
produce the chosen uniform 8 m grid; this meshing control is distinct from
the source's 25 Hz central frequency.

The five families are BrickUP viscous surface saturation, BrickUP viscous
buried groundwater, consistent full-Biot viscous surface saturation, consistent
full-Biot viscous buried groundwater, and consistent full-Biot inviscid surface
saturation. Buried groundwater is 67 m below the flat physical ground, retaining
the existing whole-cell saturation and each formulation's hydraulic boundary
discretization. Geometry, saturation, drainage and all material inputs are
fixed within a family. This does not compare the continuum equivalence of the
two formulations or their different hydraulic closures.

The homogeneous material has total density 1664 kg/m3, drained Vp
3296.26828631782 m/s and Vs 1961.16135138184 m/s. Porosity is 0.4, Biot alpha
0.5226666666666667, fluid bulk modulus 2.5e9 Pa, fluid density 1040 kg/m3,
permeability 6e-9 m2 in each direction, and full-Biot tortuosity 2.
Viscosity is 1e-3 Pa s in the viscous families and zero in the inviscid family.
Skeleton attenuation and nonlinear constitutive response are disabled; Darcy
drag remains active in the viscous families.

Receivers are 80 m horizontally from the source at depths 0, 0.5, 8 and 32 m
for surface saturation, and 59, 67.5, 75 and 99 m for buried groundwater.
Each family has a newly rerun 2048 m cube reference with distant dashpots,
the same uniform grid, timestep, executable and relative source/receiver
positions. The continuum fast-wave travel-time estimate puts the first lateral
boundary return beyond 0.5 s (at least 1968/3590.329 = 0.548 s).
The reference retains the same free/drained physical top. Distant physical
dashpot damping still selects two inertia stages for the inviscid reference;
the inviscid CPML cases use one, as in the preceding validation.

## Measured errors

Each entry is the **maximum over receivers** of full-record relative L2 error
against the matching numerical reference, in percent. Solid velocity and Darcy
flux use all three vector components. There is no filtering, fitted time shift,
gain correction or resampling. Undefined dry pressure is excluded; exactly-zero
strong surface pressure is checked in absolute terms and excluded from relative
errors. Values at different columns may correspond to different receivers.

| Family | Velocity, 0.5 Hz | Velocity, 0.0078125 Hz | Pressure, 0.5 Hz | Pressure, 0.0078125 Hz | Flux, 0.5 Hz | Flux, 0.0078125 Hz |
|---|---:|---:|---:|---:|---:|---:|
| BrickUP, surface | 10.21738% | 10.19844% | 4.15158% | 4.14330% | — | — |
| BrickUP, buried | 11.68406% | 11.66400% | 2.51392% | 2.50584% | — | — |
| Full Biot, viscous surface | 1.04168% | 1.04050% | 1.21231% | 1.21178% | 1.02772% | 1.02653% |
| Full Biot, viscous buried | 1.90384% | 1.89918% | 1.37422% | 1.37215% | 1.83131% | 1.82789% |
| Full Biot, inviscid surface | 2.12192% | 2.11927% | 12.61246% | 12.58330% | 4.84308% | 4.83560% |

The ranking minimizes the largest of velocity-vector, pressure and flux-vector
errors over all receivers. All five minima occur at the smallest tested alpha;
there is no resolved interior optimum in this range. Relative reductions in
this score are 0.185%, 0.172%, 0.044%, 0.245% and 0.231%, respectively.
Further reduction from 0.03125 to 0.0078125 Hz changes the largest score by
less than 0.02% relative in every family. The small gains do not justify
changing the default from this experiment alone.

The complete sweep of this worst-field score is:

| Alpha (Hz) | BrickUP surface | BrickUP buried | Biot viscous surface | Biot viscous buried | Biot inviscid surface |
|---:|---:|---:|---:|---:|---:|
| 0.0078125 | 10.19844% | 11.66400% | 1.21178% | 1.89918% | 12.58330% |
| 0.03125 | 10.19931% | 11.66492% | 1.21180% | 1.89952% | 12.58464% |
| 0.125 | 10.20283% | 11.66865% | 1.21189% | 1.90029% | 12.59008% |
| 0.5 | 10.21738% | 11.68406% | 1.21231% | 1.90384% | 12.61246% |
| 1 | 10.23799% | 11.70585% | 1.21302% | 1.90881% | 12.64384% |
| 2 | 10.28337% | 11.75370% | 1.21494% | 1.91953% | 12.71187% |
| 6.25 | 10.62585% | 12.02113% | 1.23047% | 1.97718% | 13.07697% |
| 12.5 | 11.42445% | 12.60466% | 1.27479% | 2.09970% | 13.81359% |
| 25 | 13.73486% | 14.45202% | 1.43703% | 2.49996% | 15.74414% |

For comparison, raising only alpha to 6.25 Hz gives velocity/pressure errors
of 10.62585%/4.28187% (BrickUP surface), 12.02113%/2.62678% (BrickUP buried),
1.06411%/1.23047% (full-Biot viscous surface), 1.97718%/1.41314% (full-Biot
viscous buried), and 2.16871%/13.07697% (full-Biot inviscid surface).
At 25 Hz the inviscid pressure error rises to 15.74414%.

Restricting the error numerator to 0.25–0.5 s, while retaining the full-record
reference norm, does not reveal a large improvement hidden by the early pulse.
For example, inviscid pressure changes from 12.60706% at 0.5 Hz to 12.57789%
at 0.0078125 Hz. Individual peak-error metrics can select a different alpha:
the inviscid late flux maximum component difference, normalized by the largest
reference component amplitude,
is smallest at 2 Hz (3.58968%, versus 3.59718% at 0.5 Hz). Such small,
metric-dependent differences do not establish a robust new default.

The earlier modified profile used reflection target 1e-4, order 3, kappa 1,
alpha 6.25 Hz and floor fraction 0. Its much smaller errors cannot be attributed
to the alpha change alone. This experiment does not isolate which of those
other controls, or their interaction, produced that improvement.

Alpha is a frequency-shift parameter, not simply an absorption-strength
multiplier. For the implemented stretch $S(s)=\kappa+d/(s+\alpha)$, at positive
angular frequency $\omega$, the imaginary-part magnitude is
$d\omega/(\alpha^2+\omega^2)$. Raising alpha while holding the other profile
controls fixed reduces this quantity. This is consistent with the measured
trend but does not by itself predict discrete reflection or stability in all
materials, frequency bands and geometries.

## Ten-second controls

Each family was extended to 10 s at both 0.5 and 0.0078125 Hz. Within a pair,
only alpha differs; relative to its short counterpart, only the end time changes.
The same force table ends at zero and remains zero after 0.5 s. All ten controls
complete 100000 steps with finite eligible station fields and no sustained
growth of the receiver envelopes over the measured duration.

The table gives the largest final-half-second peak divided by that receiver
and field's first-half-second peak, over solid velocity magnitude, absolute
pressure and Darcy-flux magnitude where available. Values are **ratios**, not
percentages; zero-pressure traces are excluded.

| Family | Default alpha, 0.5 Hz | Low alpha, 0.0078125 Hz |
|---|---:|---:|
| BrickUP, surface | 2.61985e-2 | 2.62006e-2 |
| BrickUP, buried | 9.72296e-3 | 9.71884e-3 |
| Full Biot, viscous surface | 1.90970e-11 | 6.53631e-11 |
| Full Biot, viscous buried | 4.60489e-4 | 4.59978e-4 |
| Full Biot, inviscid surface | 3.20374e-4 | 3.20551e-4 |

There is no material improvement in these late envelopes from lowering alpha.
BrickUP still has appreciable residual oscillations, despite broad decay; the
inviscid and buried full-Biot envelopes also fluctuate rather than decreasing
monotonically. These are receiver-tail checks, not long-time accuracy against
a reference or a guarantee beyond 10 s. The larger-domain reference comparisons
stop at 0.5 s. No acoustic-interface stability restriction is changed.

## Repeatability check

The first 0.5 s of each long run was compared with its short counterpart.
For four families, the largest relative L2 difference among all finite station
columns is below 4.1e-9. Buried full Biot differs more: the maximum scalar-column
difference is 1.8701e-4 at default alpha and 2.5299e-4 at low alpha. These
diagnostic discrepancies were retained, not treated as exact agreement.

Two additional paired short repeats, with byte-identical solver inputs,
confirm that the buried case also varies between runs of the same duration.
Across the original, long-run prefix and two short repeats, its worst-field
reference-error score ranges from 1.903838% to 1.903994% at default alpha and
1.899168% to 1.899308% at low alpha. The small low-alpha advantage persists
in both repeated pairs. However, same-alpha waveform differences reach
0.02010% for velocity vectors, 0.04583% for pressure and 0.01422% for flux
vectors, relative to the original waveform norms. The cause of this
repeat variation has not been isolated; a per-solve residual tolerance is not
a bound on accumulated waveform differences or a bitwise reproducibility
guarantee. The last digits in the tables should not be treated as precise
estimates of a universal optimum.

## Reproducibility and scope

The executable matches all 107 frozen production/build source hashes at
v13.12 (`6113cb10ab2e83d4ba4d96de4bfa9f6de0c76069`). It reuses the preceding
clean Horizon build, verified before execution, with binary SHA256
`1777c3925b8184fd9ecbfd5d6dd6a5cee116f1a734998bcebbdb3062553dc893`.
All jobs use `ml reset`, checked library dependencies, CUDA-aware MPI and
reservation `cuda-615-SC`. Compact cases use one node/four ranks; numerical
references use four nodes/16 ranks. No existing job reservation was changed.

Initial sweeps are jobs 97593–97597, references 97598–97602, and lower-alpha
extensions 97608–97612, in the family order above. The 45 short CPML runs and
five references complete with finite eligible fields, correct step/sample
counts, matching frozen inputs and the recorded executable. Every full-Biot
solve satisfies its reported 1e-6 true-relative-residual tolerance. Exact
input-byte checks allow only the alpha line to differ within each sweep;
auxiliary files are byte-identical. The shared discrete force history SHA256 is
`4a4b97038b4cb89bc1a913c9f4e89bfe3e3e6f17e6391376953bc5020bf30a98`.

Long controls are jobs 97617/97618 (BrickUP surface), 97619/97620 (BrickUP
buried), 97621/97622 (Biot viscous surface), 97623/97624 (Biot viscous buried),
and 97625/97626 (Biot inviscid surface), with default/low alpha in each pair.
They use the same one-node/four-rank execution settings as the short controls.
The two buried-Biot repeat pairs are jobs 97665 and 97666. All 64 runs complete;
frozen-input/binary/runtime checks, field finiteness, sample counts and reported
PCG tolerance checks pass. The separate repeatability observations above remain
part of the result. No solver correction was made during this investigation.

These are finite-mesh numerical boundary comparisons, not new Wolf analytical
validations or measured plane-wave reflection coefficients. Common dispersion,
pressure stabilization and hydraulic-boundary approximation errors remain.
Results do not qualify terrain, acoustic open-pore interfaces, sealed media,
arbitrary heterogeneity, nonlinear/viscoelastic interiors or elastic/acoustic
CPML tuning. No solver, parser, interface or runtime default changes are made.
Investigation scripts, fixtures, raw histories and generated figures remain
local-only and are not available in a clean checkout.
