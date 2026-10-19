#!/usr/bin/env python3
"""Edge-taper any Hercules topography.in so it is EXACTLY FLAT in (and
well before) the absorbing-boundary / C-PML band.

WHY YOU WANT THIS
-----------------
Hercules' C-PML absorbing layer (and SPECFEM-style codes generally)
assumes a regular, axis-aligned, flat-topped buffer: real relief is NOT
carried into the PML.  For VT (`computation_method=vt`) Hercules has a
"flatten-PML" re-route that *forces* the lateral-PML elements flat so the
run stays stable -- but that is a STABILITY BACKSTOP, not free accuracy:
if your real topography is still several hundred metres high right at the
PML inner edge, flattening it there manufactures an artificial
~few-hundred-metre "cliff" (a strong spurious scatterer) at that edge.

The standard, robust fix (SCEC/SPECFEM common practice) is to TAPER your
topography smoothly down to a flat reference inside an interior buffer,
so that it is already flat -- with margin -- *before* the PML inner edge.
Then the PML sees exactly what it expects, the flatten-PML re-route is a
no-op-like collar, and there is no cliff.  This utility does that taper
for ANY topography.in, including real DEM-derived files.

  recommended:  taper to flat  >>  rely on flatten-PML alone

WHAT IT DOES
------------
Reads a topography.in, multiplies the *deviation from a flat reference
z_ref* by a 2-D raised-cosine (Tukey) edge window keyed on the distance
to the nearest lateral boundary r = min(x, lx-x, y, ly-y):

    r >= buffer                      -> factor 1     (interior: UNCHANGED)
    flat_margin < r < buffer         -> 0.5-0.5cos() (smooth C1 ramp)
    r <= flat_margin                 -> factor 0     (DEAD-FLAT collar)

    z_out = z_ref + factor * (z_in - z_ref)

so the outer `flat_margin`-wide collar (which must cover the whole PML
band, with slack) is bit-exactly the flat reference, the ramp is smooth
(C1, no slope discontinuity -> no ramp-induced scattering), and the
deep interior is untouched.

topography.in FORMAT (Hercules convention; unchanged on output)
---------------------------------------------------------------
  line 1            : grid spacing d (m), SAME in x and y
  lines 2 .. Nx*Ny+1: one elevation (m) per grid point
Ordering: first point = SW corner (x=0,y=0); OUTER loop x (south->north),
INNER loop y (west->east). The solver uses flat_index = x_index * ny + y_index.

DEFAULTS (40 km x 40 km benchmark domain; surface mesh ~78 m, PML ~10
layers ~ 780 m): buffer 6000 m, flat_margin 1500 m (>= the ~780 m PML
band, with margin), z_ref 0.0.  Tune to your domain / PML thickness.
"""

# USAGE (Python 3 standard library; run this example from the repository root)
#   mkdir -p inputfiles/topography_pml
#   python3 script/taper_topography.py inputfiles/topography/topography.in \
#       --out inputfiles/topography_pml/topography.in \
#       --lx 40000 --ly 40000 --buffer 6000 --flat-margin 1500 \
#       --pml-thickness 780 --z-ref 0
#   python3 script/taper_topography.py --help
#
# Required: input filename. Supply your actual domain dimensions and resolved
# CPML width; defaults describe a benchmark, not an arbitrary production case.
# Defaults (metres): --lx/--ly 40000, --buffer 6000, --flat-margin 1500,
# --pml-thickness 780, --z-ref 0. --out defaults to topography_tapered.in.
# --z-ref auto uses the mean elevation of the first interior ring. It is an
# ELEVATION reference, positive up, not a solver Z coordinate.
# --buffer includes the flat collar plus the cosine ramp. Require
# 0 < flat-margin < buffer < min(lx,ly)/2. Choose flat-margin > actual CPML
# thickness with margin; --pml-thickness only controls a warning, NOT the mesh.
# Input/output: spacing followed by X-outer/Y-inner nodal elevations. Spacing
# and ordering are preserved; all output elevations round to four decimals,
# including interior values. Whole grids are held in memory.
# Create the output parent first. Existing output/in-place writes are refused
# unless --force is supplied. Relative paths use the working directory.
# Reports: grid size, taper settings, affected point count, maximum change,
# and perimeter flatness. These are geometry checks, not a stability proof.
# Set topography_directory to the new folder and rebuild the mesh when adopting
# it. This script does NOT update seabed_geometry_file or material/source files.
# Keep paired topography/bed elevations consistent: bed above the outer surface
# is invalid. Tapering both to zero removes water in the flat collar; preserving
# water requires a consistent submerged bed target. See script/README.md.

import argparse
import math
import os
import sys


def edge_factor(x, y, lx, ly, buffer_w, flat_margin):
    """2-D raised-cosine edge window: 0 in the boundary collar, smoothly
    up to 1 at the taper-start, exactly 1 in the interior."""
    r = min(x, lx - x, y, ly - y)          # dist to nearest lateral edge
    if r >= buffer_w:
        return 1.0                          # interior: untouched
    if r <= flat_margin:
        return 0.0                          # dead-flat collar (incl. edge)
    s = (r - flat_margin) / (buffer_w - flat_margin)   # 0..1
    return 0.5 - 0.5 * math.cos(math.pi * s)           # C1 raised cosine


def main():
    ap = argparse.ArgumentParser(
        description="Edge-taper a Hercules topography.in to flat before "
                    "the PML (recommended over relying on flatten-PML).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("infile", help="input topography.in")
    ap.add_argument("--out", default="topography_tapered.in",
                    help="output file (must differ from infile unless "
                         "--force)")
    ap.add_argument("--force", action="store_true",
                    help="allow overwriting infile / an existing --out")
    ap.add_argument("--lx", type=float, default=40000.0,
                    help="domain length S-N (x), m")
    ap.add_argument("--ly", type=float, default=40000.0,
                    help="domain length W-E (y), m")
    ap.add_argument("--buffer", type=float, default=6000.0,
                    help="total taper band width from each lateral "
                         "boundary inward (m): collar + ramp")
    ap.add_argument("--flat-margin", type=float, default=1500.0,
                    help="width of the GUARANTEED-FLAT collar adjacent "
                         "to each boundary (m); must be >= the PML "
                         "thickness, with margin")
    ap.add_argument("--pml-thickness", type=float, default=780.0,
                    help="PML band thickness (m), for the safety check "
                         "only (default ~10 layers x ~78 m)")
    ap.add_argument("--z-ref", default="0.0",
                    help="flat reference: a number (m), or 'auto' = mean "
                         "elevation of the first interior ring just "
                         "inside the taper band")
    a = ap.parse_args()

    if not all(math.isfinite(value) for value in
               (a.lx, a.ly, a.buffer, a.flat_margin, a.pml_thickness)):
        sys.exit("error: domain and taper dimensions must be finite")
    if min(a.lx, a.ly) <= 0 or a.pml_thickness < 0:
        sys.exit("error: domain lengths must be positive and PML thickness nonnegative")

    if not os.path.isfile(a.infile):
        sys.exit("error: input file not found: %s" % a.infile)
    if os.path.abspath(a.out) == os.path.abspath(a.infile) and not a.force:
        sys.exit("error: --out equals infile; pass --force to overwrite "
                 "in place (or choose a different --out)")
    if os.path.exists(a.out) and os.path.abspath(a.out) != os.path.abspath(
            a.infile) and not a.force:
        sys.exit("error: --out already exists: %s (pass --force)" % a.out)

    if a.flat_margin <= 0 or a.buffer <= a.flat_margin:
        sys.exit("error: need 0 < flat-margin < buffer (got flat-margin="
                 "%g, buffer=%g)" % (a.flat_margin, a.buffer))
    if 2.0 * a.buffer >= min(a.lx, a.ly):
        sys.exit("error: 2*buffer (%g) >= min domain extent (%g): the "
                 "opposite tapers would overlap. Reduce --buffer or use "
                 "a larger domain." % (2.0 * a.buffer, min(a.lx, a.ly)))
    if a.flat_margin < a.pml_thickness:
        print("WARNING: flat-margin %g m < PML thickness %g m -- the "
              "guaranteed-flat collar does not fully cover the PML band. "
              "Increase --flat-margin to >= the PML thickness (with "
              "margin) so the PML sees only flat topography."
              % (a.flat_margin, a.pml_thickness), file=sys.stderr)

    with open(a.infile) as f:
        values = [float(value) for value in f.read().split()]
    if len(values) < 2 or not all(math.isfinite(value) for value in values):
        sys.exit("error: need spacing and finite elevation values")
    spacing, z = values[0], values[1:]
    if spacing <= 0:
        sys.exit("error: spacing must be positive")

    nx = int(round(a.lx / spacing)) + 1
    ny = int(round(a.ly / spacing)) + 1
    if abs(nx - 1 - a.lx / spacing) > 1e-9 or abs(
            ny - 1 - a.ly / spacing) > 1e-9:
        sys.exit("error: spacing %g must divide --lx %g and --ly %g "
                 "exactly" % (spacing, a.lx, a.ly))
    if len(z) != nx * ny:
        sys.exit("error: read %d elevations but --lx/--ly/spacing imply "
                 "%d x %d = %d. Pass the correct --lx/--ly (this file's "
                 "spacing is %g)." % (len(z), nx, ny, nx * ny, spacing))

    # resolve z_ref (number, or 'auto' = mean of the first interior ring
    # just inside the taper band -> the ramp connects to a representative
    # interior level, minimising the tapered amplitude for offset DEMs)
    if str(a.z_ref).lower() == "auto":
        ring, lo, hi = [], a.buffer, a.buffer + spacing
        for k, zk in enumerate(z):
            i, j = divmod(k, ny)
            r = min(i * spacing, a.lx - i * spacing,
                    j * spacing, a.ly - j * spacing)
            if lo <= r < hi:
                ring.append(zk)
        z_ref = (sum(ring) / len(ring)) if ring else (sum(z) / len(z))
        zref_src = "auto (mean of first interior ring = %.4f m)" % z_ref
    else:
        z_ref = float(a.z_ref)
        if not math.isfinite(z_ref):
            sys.exit("error: z-ref must be finite")
        zref_src = "%.4f m (fixed)" % z_ref

    max_dz, n_mod, max_perim = 0.0, 0, 0.0
    out = [None] * len(z)
    for k, zk in enumerate(z):
        i, j = divmod(k, ny)
        x, y = i * spacing, j * spacing
        fac = edge_factor(x, y, a.lx, a.ly, a.buffer, a.flat_margin)
        zo = z_ref + fac * (zk - z_ref)
        out[k] = zo
        if fac < 1.0:
            n_mod += 1
            max_dz = max(max_dz, abs(zo - zk))
        if i == 0 or i == nx - 1 or j == 0 or j == ny - 1:
            max_perim = max(max_perim, abs(zo - z_ref))

    with open(a.out, "w") as f:
        f.write("%g\n" % spacing)
        for zo in out:
            f.write("%.4f\n" % zo)

    ramp = a.buffer - a.flat_margin
    print("wrote %s: spacing=%g m, grid %d x %d = %d points (+1 spacing "
          "line)" % (a.out, spacing, nx, ny, nx * ny))
    print("taper: buffer=%g m  (flat collar=%g m + cosine ramp=%g m)  "
          "z_ref=%s" % (a.buffer, a.flat_margin, ramp, zref_src))
    print("modified %d / %d points (interior %d untouched); max |dz| = "
          "%.4f m" % (n_mod, len(z), len(z) - n_mod, max_dz))
    print("perimeter flatness check: max |z_boundary - z_ref| = %.6g m "
          "(must be 0)  ->  %s"
          % (max_perim, "OK" if max_perim == 0.0 else "FAIL"))
    print("The outer %g m (>= PML band) is now EXACTLY flat; the PML / "
          "flatten-PML re-route sees no relief -> no artificial cliff."
          % a.flat_margin)


if __name__ == "__main__":
    main()
