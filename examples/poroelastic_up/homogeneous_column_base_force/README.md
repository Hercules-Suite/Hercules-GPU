# Homogeneous BrickUP Column With Base Force

This is a minimal OpenSees-style validation column for the poroelastic U-P
implementation and the distributed compliant-base source.

Use `parameters.in` as the single input file for the current solver interface:

```bash
mpirun -np <N> <path-to-psolve> parameters.in
```

The split `physics.in` and `numerical.in` files are kept only as reference for
older scripts.
This example assumes a PROJ-enabled build, so it does not include legacy
`domain_surface_corners` blocks.

Geometry and material:

- Domain: 1 m north-south x 1 m east-west x 32 m depth
- Homogeneous profile: `Vs = 200 m/s`, `Vp = 400 m/s`,
  `rho = 2000 kg/m^3`
- Poroelastic depth: 32 m, so the full column is eligible when
  `include_poroelastic = yes`
- Groundwater: `surface_groundwatertable = no`, meaning the groundwater level
  is the local ground surface and the full shallow column is saturated

The timestep is intentionally conservative for the U-P case. With
`poroelastic_fluid_bulk_modulus_pa = 2.2e9` and `poroelastic_porosity = 0.4`,
the undrained storage stiffness is much larger than the dry skeleton bulk
modulus, so the fast compressional speed is roughly 1.7 km/s rather than the
dry `Vp = 400 m/s`. Do not choose `simulation_delta_time_sec` from the dry
elastic CFL alone.

The source is `type_of_source = base_force`. It reads
`sourcefiles/baseforcefunction.in` as a base velocity history in m/s and applies
a nodal force distribution on the deepest mesh-node plane. The included table
is a sampled Ricker velocity pulse with peak frequency 10 Hz, peak time 0.12 s,
peak velocity 0.01 m/s, and a 0.30 s table duration so both tails are near zero:

```text
sum(F_base) = scale * rho_base * Vs_base * A_base * v(t)
```

For this 1 m x 1 m column, `A_base = 1 m^2`. The example sets
`base_force_rho_kgm3 = 2000` and `base_force_vs_mps = 200`, matching the
homogeneous material. Keep these explicit for OpenSees comparisons: they are
the dashpot impedance used to convert velocity to force, and in layered
problems they may intentionally represent the bedrock below the modeled soil.

Use `base_force_scale = 1.0` for an incident/within velocity history and
`2.0` for an outcrop velocity history. The default direction here is horizontal
`x` (north), matching a simple shear-wave site-response test. To do an
undrained pore-pressure sign check, switch the direction to vertical
`z` and use a compressional velocity pulse.

Suggested first checks on Vista:

1. Run as provided and compare horizontal station velocity/displacement against
   an OpenSees elastic/BrickUP column with the same base dashpot force.
2. Change `include_poroelastic = no` and `print_station_pressures = no`; for
   this shear-only, homogeneous, linear smoke test, the elastic motion should be
   very close unless the pore pressure field is intentionally excited.
3. Change `base_force_direction_x = 0` and `base_force_direction_z = 1` for a
   vertical compression check; compression should produce positive excess pore
   pressure under the current sign convention.

No local build is expected from this workspace; compile and run on Vista.
