#!/usr/bin/env python3
"""Generate a Hercules topography.in: a gentle planar tilt that reaches
ALL FOUR lateral boundaries -- purpose-built to EXERCISE the C-PML
flatten-PML re-route (cpml_pml_flattened_element).

Why a tilt (not a hill): the earlier Gaussian-hill example decays to
~0 over the PML band, so NO PML element is topography-intersected and
the flatten-PML path never fires (flattened-pml-skipped ~ 0). A planar
tilt is non-zero and sloping EVERYWHERE, including the entire lateral
PML perimeter, so essentially every near-surface lateral-PML element is
topography-deformed -> for computation_method=vt the flatten-PML
predicate returns 1 for all of them (flattened-pml-skipped is large and
uniform). That makes this the actual functional test of the re-route.

Kept deliberately GENTLE (default ~1 deg slope) so we test the re-route
MECHANISM, not steep-boundary pathology (steep relief at the boundary
is the separately-documented risky regime).

Hercules topography.in format
-----------------------------
  line 1            : grid spacing d (m), SAME in x and y
  lines 2 .. N*M+1  : one elevation (m) per grid point
Ordering (x = south->north, y = west->east): first point = SW corner
(x=0,y=0); finish ALL points sharing a y before advancing y, i.e.
  for j in 0..My-1:        # y : west -> east   (OUTER)
      y = j*d
      for i in 0..Mx-1:    # x : south -> north (INNER)
          x = i*d
          write elevation(x,y)

Surface (default): pure plane
  elev(x,y) = z0 + sx*x + sy*y          (sx,sy >= 0, elevation >= 0)
defaults sx = sy = 0.0125 m/m, z0 = 0  on the 40000 x 40000 m domain
=> SW corner 0 m, NE corner 1000 m, |grad| ~ 0.0177 -> ~1.01 deg.
Elevation along every lateral boundary is non-zero & sloping (e.g. the
east edge y=40000: elev = sx*x + 500, i.e. 500..1000 m), so all four
lateral PML bands are clearly topography-deformed.

Optional interior Gaussian bump (default OFF, --bump-h 0) if you also
want an interior feature; for a clean re-route mechanism test leave it
off (pure tilt = simplest to interpret).

Optional --edge-taper (default OFF) smoothly brings the surface to a
flat reference inside an interior buffer, so it is already flat -- with
margin -- BEFORE the PML inner edge (the recommended real-model
workflow; see ../taper_topography.py). Default OFF keeps this file's
deliberate worst-case (boundary-reaching) test surface exactly
reproducible.
"""
import argparse, math


def _edge_taper_factor(x, y, lx, ly, buffer_w, flat_margin):
    """2-D raised-cosine edge window: 0 in the boundary collar, smoothly
    up to 1 at the taper-start, exactly 1 in the interior. See
    ../taper_topography.py for the standalone (real-DEM) version."""
    r = min(x, lx - x, y, ly - y)
    if r >= buffer_w:
        return 1.0
    if r <= flat_margin:
        return 0.0
    s = (r - flat_margin) / (buffer_w - flat_margin)
    return 0.5 - 0.5 * math.cos(math.pi * s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=float, default=40000.0, help="domain S-N (x), m")
    ap.add_argument("--ly", type=float, default=40000.0, help="domain W-E (y), m")
    ap.add_argument("--spacing", type=float, default=100.0,
                    help="grid spacing d (m); must divide lx & ly")
    ap.add_argument("--sx", type=float, default=0.0125,
                    help="south->north slope (m/m)")
    ap.add_argument("--sy", type=float, default=0.0125,
                    help="west->east slope (m/m)")
    ap.add_argument("--z0", type=float, default=0.0,
                    help="elevation at SW corner (m)")
    ap.add_argument("--bump-h", type=float, default=0.0,
                    help="optional interior Gaussian bump height (m); 0=off")
    ap.add_argument("--bump-sigma", type=float, default=3000.0)
    ap.add_argument("--bump-x0", type=float, default=20000.0)
    ap.add_argument("--bump-y0", type=float, default=20000.0)
    ap.add_argument("--edge-taper", action="store_true",
                    help="taper the surface to a flat reference before "
                         "the PML (recommended real-model workflow; "
                         "OFF by default so the worst-case test surface "
                         "stays reproducible)")
    ap.add_argument("--taper-buffer", type=float, default=6000.0,
                    help="total taper band width from each lateral "
                         "boundary inward (m): flat collar + cosine ramp")
    ap.add_argument("--taper-flat-margin", type=float, default=1500.0,
                    help="guaranteed-flat collar width (m); must be >= "
                         "the PML thickness, with margin")
    ap.add_argument("--taper-zref", type=float, default=0.0,
                    help="flat reference elevation (m) the taper goes to")
    ap.add_argument("--out", default="topography.in")
    a = ap.parse_args()

    nx = int(round(a.lx / a.spacing)) + 1
    ny = int(round(a.ly / a.spacing)) + 1
    if abs(nx - 1 - a.lx / a.spacing) > 1e-9 or abs(ny - 1 - a.ly / a.spacing) > 1e-9:
        raise SystemExit("spacing must divide lx and ly exactly")

    binv = (1.0 / (2.0 * a.bump_sigma ** 2)) if a.bump_h != 0.0 else 0.0
    slope_deg = math.degrees(math.atan(math.hypot(a.sx, a.sy)))
    with open(a.out, "w") as f:
        f.write(f"{a.spacing:g}\n")                       # line 1: spacing
        for j in range(ny):                               # y: west -> east
            y = j * a.spacing
            for i in range(nx):                           # x: south -> north
                x = i * a.spacing
                z = a.z0 + a.sx * x + a.sy * y
                if a.bump_h != 0.0:
                    z += a.bump_h * math.exp(
                        -(((x - a.bump_x0) ** 2) + ((y - a.bump_y0) ** 2)) * binv)
                if a.edge_taper:
                    fac = _edge_taper_factor(x, y, a.lx, a.ly,
                                             a.taper_buffer,
                                             a.taper_flat_margin)
                    z = a.taper_zref + fac * (z - a.taper_zref)
                f.write(f"{z:.4f}\n")

    ne = a.z0 + a.sx * a.lx + a.sy * a.ly
    print(f"wrote {a.out}: spacing={a.spacing:g} m, grid {nx} x {ny} "
          f"= {nx*ny} points (+1 spacing line)")
    print(f"plane: z0={a.z0:g}  sx={a.sx:g}  sy={a.sy:g} m/m  "
          f"=> SW={a.z0:g} m, NE={ne:g} m, slope ~ {slope_deg:.2f} deg"
          + (f"  + interior bump H={a.bump_h:g} m" if a.bump_h else ""))
    if a.edge_taper:
        print(f"EDGE-TAPER ON: flat to z={a.taper_zref:g} m within "
              f"buffer={a.taper_buffer:g} m (collar={a.taper_flat_margin:g} m "
              f"+ ramp={a.taper_buffer - a.taper_flat_margin:g} m) of every "
              f"lateral boundary -> surface is flat BEFORE the PML "
              f"(recommended real-model setup; the boundary is no longer a "
              f"flatten-PML stress test).")
    else:
        print("Non-zero & sloping along ALL four lateral boundaries -> the "
              "flatten-PML re-route (computation_method=vt) is fully "
              "exercised (expect a LARGE flattened-pml-skipped count). "
              "This is a DELIBERATE worst-case test, NOT a model-setup "
              "template -- for real models add --edge-taper (or run "
              "../taper_topography.py).")
    print("ordering: line1=spacing; OUTER y(W->E), INNER x(S->N); "
          "first point = SW corner (x=0,y=0)")

if __name__ == "__main__":
    main()
