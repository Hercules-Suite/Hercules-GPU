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

SEABED PAIRS
------------
If seabed.in exists beside the input surface (or --seabed supplies a file),
check the tapered pair with the original wet/dry nodes before publishing.
Keep an already valid bed byte-for-byte; otherwise repair originally dry beds
to the new surface, retaining wet beds. Invalid original geometry, changed
water levels and invalid repaired wet cells stop the command without publishing.
Both outputs and repair-summary.json go beside --out. Keep original inputs in
a separate directory. --no-seabed explicitly requests surface-only operation.
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
# unless --force is supplied. Paired outputs must never replace either input,
# even with --force. Output symlinks are refused. Relative paths use the cwd.
# Reports: grid size, taper settings, affected point count, maximum change,
# and perimeter flatness. These are geometry checks, not a stability proof.
# Auto-detect infile's sibling seabed.in, or use --seabed PATH. --no-seabed skips
# detection/check/repair. Paired mode requires repair_tapered_geometry.py beside
# this script (both ship in script/). It uses that helper's fixed 1e-8 m checks,
# independent of parameters.in. --topography-base optionally supplies the run's
# nonnegative datum for the absolute top-plane check; omitted means skip that
# check. Wet/dry, crossing and wet-cell flatness checks still run.
# Paired output: --out, sibling seabed.in, repair-summary.json (valid/repaired,
# counts, maximum adjustment and hashes). Originally wet beds retain precision;
# dry beds follow the tapered surface when repair is needed. Changed water
# levels, malformed originals or invalid repaired cells stop before publishing.
# Set BOTH topography_directory and seabed_geometry_file to the output pair;
# rebuild meshes/sources and review placements. No parameters.in or material/
# source inputs are edited. Validation is geometric, not a propagation test.

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile


def same_path(first, second):
    first, second = Path(first), Path(second)
    return first.resolve() == second.resolve() or (
        first.exists() and second.exists() and os.path.samefile(first, second))


def check_outputs(paths, force, protected=()):
    for i, path in enumerate(paths):
        if not path.parent.is_dir():
            raise ValueError('Output parent directory does not exist: ' + str(path.parent))
        if any(same_path(path, other) for other in paths[:i]):
            raise ValueError('Output filenames must be distinct: ' + str(path))
        if any(same_path(path, source) for source in protected):
            raise ValueError('Paired output would replace an input: %s; '
                             'choose --out in a separate directory, even with --force' % path)
        if path.is_symlink() or path.is_dir():
            raise ValueError('Output cannot be a symlink or directory: ' + str(path))
        if path.exists() and not force:
            raise ValueError('Output already exists: %s (pass --force to replace outputs)' % path)


def publish_outputs(staged, destinations, force, protected=()):
    """Publish only validated files; roll back ordinary I/O failures."""
    check_outputs(destinations, force, protected)
    with tempfile.TemporaryDirectory(prefix='.taper-backup-', dir=destinations[0].parent) as tmp:
        backups, published = {}, []
        try:
            for i, (source, destination) in enumerate(zip(staged, destinations)):
                if force:
                    if destination.exists():
                        backup = Path(tmp)/str(i)
                        shutil.copy2(destination, backup)
                        backups[destination] = backup
                    os.replace(source, destination)
                else:
                    # Exclusive creation: do not overwrite a path created since preflight.
                    os.link(source, destination)
                published.append(destination)
        except OSError:
            for destination in reversed(published):
                if destination in backups:
                    os.replace(backups[destination], destination)
                else:
                    destination.unlink()
            raise


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
                    help="replace existing outputs; in-place input replacement "
                         "is allowed only without seabed handling")
    bed_options = ap.add_mutually_exclusive_group()
    bed_options.add_argument("--seabed", type=Path,
                             help="original bed file; otherwise auto-detect seabed.in beside infile")
    bed_options.add_argument("--no-seabed", action="store_true",
                             help="explicitly disable seabed detection, checking and repair")
    ap.add_argument("--topography-base", type=float,
                    help="optional topographybase_zcoord (m, nonnegative); "
                         "with seabed, check wet surfaces against this datum when positive")
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
    if a.topography_base is not None and (
            not math.isfinite(a.topography_base) or a.topography_base < 0):
        ap.error('--topography-base must be finite and nonnegative')
    bed = None
    if not a.no_seabed:
        candidate = a.seabed if a.seabed is not None else Path(a.infile).parent/'seabed.in'
        if a.seabed is not None or candidate.exists() or candidate.is_symlink():
            if not candidate.is_file():
                ap.error('Seabed input is not a readable file: ' + str(candidate))
            bed = candidate
    destination = Path(a.out)
    destinations = [destination]
    protected = []
    if bed is not None:
        destinations += [destination.parent/'seabed.in', destination.parent/'repair-summary.json']
        protected = [Path(a.infile), bed]
    check_outputs(destinations, a.force, protected)

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

    with open(a.infile, 'rb') as f:
        raw = f.read()
    surface_hash = hashlib.sha256(raw).hexdigest()
    values = [float(value) for value in raw.split()]
    del raw
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

    summary = None
    with tempfile.TemporaryDirectory(prefix='.taper-geometry-', dir=destination.parent) as tmp:
        tmp = Path(tmp)
        staged_surface = tmp/'topography.in'
        with staged_surface.open('w') as f:
            f.write("%g\n" % spacing)
            for zo in out:
                f.write("%.4f\n" % zo)
        staged = [staged_surface]
        if bed is not None:
            # Import only for paired geometry: surface-only use remains standalone.
            try:
                from repair_tapered_geometry import repair_geometry
            except ImportError as error:
                raise ValueError('Paired mode requires repair_tapered_geometry.py beside '
                                 'taper_topography.py; copy both utilities from script/') from error
            print('Seabed detected: %s; checking the tapered pair before publishing.' % bed,
                  flush=True)
            pair = tmp/'pair'
            summary = repair_geometry(Path(a.infile), bed, staged_surface, pair,
                                      a.lx, a.ly, a.topography_base, only_if_needed=True)
            if summary['input_sha256'][str(Path(a.infile))] != surface_hash:
                raise ValueError('Original surface changed during tapering')
            summary['tapered_surface_sha256'] = summary['input_sha256'].pop(str(staged_surface))
            summary['output_files'] = {name: str(path.resolve()) for name, path in zip(
                ('topography.in', 'seabed.in', 'repair-summary.json'), destinations)}
            (pair/'repair-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            staged = [pair/name for name in ('topography.in', 'seabed.in', 'repair-summary.json')]
        publish_outputs(staged, destinations, a.force, protected)

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
    if summary is not None:
        if summary['status'] == 'valid':
            print('Seabed check: valid; original bed copied byte-for-byte.')
        else:
            print('Seabed check: repaired %d dry bed nodes; maximum adjustment %.9g m; '
                  'all %d original wet bed values retained.' % (
                      summary['changed_dry_bed_nodes'], summary['maximum_dry_bed_adjustment_m'],
                      summary['preserved_wet_bed_nodes']))
        print('Paired outputs: %s and %s; summary: %s' % tuple(destinations))
        if a.topography_base is None:
            print('Absolute top-plane check skipped; supply --topography-base for this check.')
    elif a.no_seabed:
        print('Seabed handling explicitly disabled (--no-seabed); surface output only.')
    else:
        print('No sibling seabed.in found; surface output only.')


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        sys.exit('error: ' + str(error))
