# OpenSees Companion Model

This directory contains a first OpenSees comparison model for the Hercules
homogeneous BrickUP column.

Run it from this directory:

```bash
OpenSees model.tcl
```

The script reads the same `../sourcefiles/baseforcefunction.in` velocity pulse
used by Hercules, converts it to total compliant-base force using

```text
F(t) = rho_base * Vs_base * A_base * v(t)
```

and distributes that force over the four bottom nodes. The model keeps the
same 1 m x 1 m x 32 m 3D block geometry as Hercules. It does not impose
OpenSees `equalDOF` column constraints; instead it approximates Hercules'
dashpot boundaries with nodal dashpot forces on the lateral and bottom faces.
Those forces are updated once per OpenSees time step as `F = -C v`.

Outputs are written to `output/station.0` through `output/station.4` at the
same station depths and output interval as the Hercules example.

Important comparison notes:

- The top pore-pressure DOF is fixed to zero to mimic
  `poroelastic_surface_pressure_bc = drained`.
- Permeability is zero, matching the current Hercules smoke test.
- OpenSees `brickUP` takes a single storage/combined bulk modulus, not a
  separate porosity. For the Hercules settings `K_f = 2.2e9 Pa`,
  `porosity = 0.4`, and `alpha = 1`, the script passes `K_f / n = 5.5e9 Pa`.
- This is a waveform and sign/magnitude comparison target, not a bitwise match:
  OpenSees and Hercules use different boundary implementations, element
  assembly paths, and time integrator details.
- To compare an undrained top, remove or comment the top-surface `fix ... 1`
  pressure constraints in `model.tcl` and change the Hercules case accordingly.
