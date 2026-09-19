#!/usr/bin/env python3
"""Generate a Hercules topography.in for a smooth Gaussian hill.

Hercules topography.in format
-----------------------------
  line 1            : grid spacing d (metres), SAME in x and y
  lines 2 .. N*M+1  : one elevation (metres) per grid point

Point ordering (per Hercules convention; x = south->north, y = west->east):
  the FIRST point is the SW (bottom-left) corner (x=0, y=0); finish ALL
  points sharing a y before advancing y.  i.e.

      for j in 0..My-1:          # y : west -> east   (OUTER)
          y = j*d
          for i in 0..Mx-1:      # x : south -> north (INNER)
              x = i*d
              write elevation(x, y)

Model domain (this example): 40000 m (E-W, y) x 40000 m (S-N, x), matching
region_length_east_m = region_length_north_m = 40000 in the elastic-
halfspace / cps benchmark, so topo-vs-flat runs are directly comparable.

Feature: a Gaussian hill, centred and INTERIOR by design so it decays to
~0 well before the lateral boundaries.  That (a) keeps relief gentle at
the absorbing boundary -- the well-behaved regime -- and (b) minimises
the topography<->PML overlap, so the CPML flatten-PML path
(computation_method=vt) is exercised cleanly and the startup diagnostic
'flattened-pml-skipped=N' should report a small/zero N.

    elev(x,y) = H * exp( -((x-x0)^2 + (y-y0)^2) / (2*sigma^2) )

Defaults: H=800 m, sigma=3000 m, centre=(20000,20000), d=100 m
=> 401 x 401 = 160801 points; max slope ~ H/(sigma*sqrt(e)) ~ 9.2 deg
(gentle/moderate -- safe for a first VT+CPML example).  Increase H or
decrease sigma for a steeper stress test; keep the feature interior.

Optional --edge-taper (default OFF) smoothly forces the surface to a
flat reference inside an interior buffer so it is flat -- with margin --
BEFORE the PML inner edge (the recommended real-model workflow; see
../../script/taper_topography.py). For THIS default interior hill it is a near
no-op (the hill already decays to ~0 well before the boundary); it
matters once you widen/relocate the hill so it reaches the boundary.
Default OFF keeps this example's topography.in exactly reproducible.
"""
import argparse, math


def _edge_taper_factor(x, y, lx, ly, buffer_w, flat_margin):
    """2-D raised-cosine edge window: 0 in the boundary collar, smoothly
    up to 1 at the taper-start, exactly 1 in the interior. See
    ../../script/taper_topography.py for the standalone (real-DEM) version."""
    r = min(x, lx - x, y, ly - y)
    if r >= buffer_w:
        return 1.0
    if r <= flat_margin:
        return 0.0
    s = (r - flat_margin) / (buffer_w - flat_margin)
    return 0.5 - 0.5 * math.cos(math.pi * s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=float, default=40000.0,
                    help="domain length S-N (x), m")
    ap.add_argument("--ly", type=float, default=40000.0,
                    help="domain length W-E (y), m")
    ap.add_argument("--spacing", type=float, default=100.0,
                    help="grid spacing d (m), same in x and y; must divide lx & ly")
    ap.add_argument("--height", type=float, default=800.0,
                    help="Gaussian hill peak elevation H (m)")
    ap.add_argument("--sigma", type=float, default=3000.0,
                    help="Gaussian half-width sigma (m)")
    ap.add_argument("--x0", type=float, default=20000.0, help="hill centre x (m)")
    ap.add_argument("--y0", type=float, default=20000.0, help="hill centre y (m)")
    ap.add_argument("--edge-taper", action="store_true",
                    help="taper the surface to a flat reference before "
                         "the PML (recommended real-model workflow; OFF "
                         "by default so this example stays reproducible)")
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

    # inclusive of both endpoints (0 and L)
    nx = int(round(a.lx / a.spacing)) + 1
    ny = int(round(a.ly / a.spacing)) + 1
    if abs(nx - 1 - a.lx / a.spacing) > 1e-9 or abs(ny - 1 - a.ly / a.spacing) > 1e-9:
        raise SystemExit("spacing must divide lx and ly exactly")

    inv = 1.0 / (2.0 * a.sigma * a.sigma)
    max_slope = a.height / (a.sigma * math.sqrt(math.e))
    with open(a.out, "w") as f:
        f.write(f"{a.spacing:g}\n")                     # line 1: spacing
        for j in range(ny):                             # y: west -> east
            y = j * a.spacing
            dy2 = (y - a.y0) ** 2
            for i in range(nx):                         # x: south -> north
                x = i * a.spacing
                z = a.height * math.exp(-((x - a.x0) ** 2 + dy2) * inv)
                if a.edge_taper:
                    fac = _edge_taper_factor(x, y, a.lx, a.ly,
                                             a.taper_buffer,
                                             a.taper_flat_margin)
                    z = a.taper_zref + fac * (z - a.taper_zref)
                f.write(f"{z:.4f}\n")

    print(f"wrote {a.out}: spacing={a.spacing:g} m, grid {nx} x {ny} "
          f"= {nx*ny} points (+1 spacing line)")
    print(f"Gaussian hill: H={a.height:g} m, sigma={a.sigma:g} m, "
          f"centre=({a.x0:g},{a.y0:g}); max slope ~ {math.degrees(math.atan(max_slope)):.2f} deg")
    if a.edge_taper:
        print(f"EDGE-TAPER ON: flat to z={a.taper_zref:g} m within "
              f"buffer={a.taper_buffer:g} m (collar={a.taper_flat_margin:g} m "
              f"+ ramp={a.taper_buffer - a.taper_flat_margin:g} m) of every "
              f"lateral boundary -> surface is flat BEFORE the PML "
              f"(near no-op for the default interior hill; matters once "
              f"the hill reaches the boundary).")
    print("ordering: line1=spacing; then OUTER y(W->E), INNER x(S->N); "
          "first point = SW corner (x=0,y=0)")

if __name__ == "__main__":
    main()
