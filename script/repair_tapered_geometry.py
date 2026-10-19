#!/usr/bin/env python3
"""Preserve original wet/dry nodes when pairing a tapered surface with a bed.

Python 3 standard library only. Elevations are metres positive upward;
ordering is X/north outer, Y/east inner, matching Hercules. See the USAGE
comments below and script/README.md for a complete command.

All arguments are required. Refuses an existing output directory. Original
files are read-only. The new directory contains topography.in (byte copy),
seabed.in, and repair-summary.json. Dry bed tokens are copied from tapered
surface tokens, and wet bed tokens from the original bed; no re-rounding.
Validates the original pair, every output cell, and unchanged water levels
at wet nodes using this helper's fixed 1e-8 m tolerance. It does not read the
solver's seabed_height_tolerance_m setting. This does not support intentional
changes to water levels or sloping acoustic free surfaces. Input is streamed;
memory scales with grid-row width and the longest input line, not grid area.
No case settings or jobs are changed.
Rebuild the mesh and review source placement before using the new geometry.
"""

# USAGE (Python 3 standard library; run from the repository root)
#   python3 script/repair_tapered_geometry.py \
#       --original-surface inputfiles/geometry/topography.in \
#       --original-bed inputfiles/geometry/seabed.in \
#       --tapered-surface inputfiles/geometry_tapered/topography.in \
#       --out-dir inputfiles/geometry_tapered_consistent \
#       --lx 40000 --ly 40000 --topography-base 1250
#   python3 script/repair_tapered_geometry.py --help
#
# All seven options are required; dimensional values above are examples.
# --lx/--ly: positive north/east domain lengths in metres. --topography-base:
# the run's nonnegative topographybase_zcoord in metres, NOT an elevation.
# Each grid: one positive common spacing, then exactly (lx/spacing+1) *
# (ly/spacing+1) finite ASCII elevations, positive UP, X-outer/Y-inner.
# Whitespace separates values; no count header, coordinates or comments.
# Both original grids must use the same origin, datum and grid as the taper.
# Original wet nodes: gap > 1e-8 m. Their bed tokens are preserved; the new
# water level must agree within 1e-8 m and the node must remain wet.
# Original dry nodes: copy the tapered-surface token into the repaired bed.
# This removes rounding-induced water gaps without re-rounding elevations.
# The fixed 1e-8 m check is independent of solver settings; there are no
# optional repair controls. Existing output directories/symlinks are refused,
# with no --force. Create the output parent first; original inputs are read-only.
# Output: topography.in (byte copy of taper), seabed.in, repair-summary.json.
# Summary: wet/dry counts, changed dry points, maximum adjustment, wet-cell
# counts and input/output SHA-256 hashes. Validation failures publish no pair.
# Limits: preserves node classifications, not an identical interpolated
# coastline. Does not flatten changed water levels, check acoustic mesh-face
# alignment, or establish simulation accuracy/stability. See script/README.md.

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile

TOL = 1e-8


class Grid:
    def __init__(self, path):
        self.path = Path(path)
        self.digest = hashlib.sha256()
        self.stream = self.tokens()
        self.spacing = self.take()[0]
        if self.spacing <= 0:
            raise ValueError('Grid spacing must be positive: ' + str(self.path))

    def tokens(self):
        with self.path.open('rb') as source:
            for line in source:
                self.digest.update(line)
                for token in line.split():
                    yield token.decode('ascii')

    def take(self):
        try:
            token = next(self.stream)
        except StopIteration:
            raise ValueError('Too few grid values: ' + str(self.path)) from None
        try:
            # Python permits digit separators; Hercules' numeric reader does not.
            if '_' in token:
                raise ValueError
            value = float(token)
        except ValueError:
            raise ValueError('Invalid numeric token %r in %s' % (token, self.path)) from None
        if not math.isfinite(value):
            raise ValueError('Nonfinite grid value: ' + str(self.path))
        return value, token

    def finish(self):
        if next(self.stream, None) is not None:
            raise ValueError('Extra grid values: ' + str(self.path))
        return self.digest.hexdigest()


def validate_cells(previous, current, row, label, base):
    wet_count = 0
    if previous is None:
        return wet_count
    for j in range(len(current)-1):
        nodes = (previous[j], previous[j+1], current[j], current[j+1])
        if any(surface-bed > TOL for surface, bed in nodes):
            wet_count += 1
            low = min(surface for surface, bed in nodes)
            high = max(surface for surface, bed in nodes)
            if high-low > TOL or (base > 0 and high >= base):
                raise ValueError('%s: invalid wet cell (%d,%d), surface range %g..%g m'
                                 % (label, row-1, j, low, high))
    return wet_count


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    for name, help_text in (
        ('original-surface', 'untapered topography.in; elevations in metres, positive up'),
        ('original-bed', 'matching original seabed.in; elevations in metres, positive up'),
        ('tapered-surface', 'tapered topography.in on the same grid and datum'),
        ('out-dir', 'new output directory; parent must exist; never overwritten'),
    ):
        parser.add_argument('--'+name, type=Path, required=True, help=help_text)
    for name, help_text in (
        ('lx', 'north/X domain length in metres, positive'),
        ('ly', 'east/Y domain length in metres, positive'),
        ('topography-base', 'topographybase_zcoord in metres, nonnegative'),
    ):
        parser.add_argument('--'+name, type=float, required=True, help=help_text)
    args = parser.parse_args()
    if not all(math.isfinite(v) for v in (args.lx, args.ly, args.topography_base)):
        parser.error('Domain dimensions and base must be finite')
    if min(args.lx, args.ly) <= 0 or args.topography_base < 0:
        parser.error('Domain dimensions must be positive; base must be nonnegative')
    out = args.out_dir
    if out.exists() or out.is_symlink():
        parser.error('Output already exists; choose a new directory: '+str(out))
    if not out.parent.is_dir():
        parser.error('Output parent directory does not exist: '+str(out.parent))
    old_surface, old_bed = Grid(args.original_surface), Grid(args.original_bed)
    new_surface = Grid(args.tapered_surface)
    spacing = old_surface.spacing
    if any(abs(g.spacing-spacing) > 1e-10*max(1, spacing)
           for g in (old_bed, new_surface)):
        raise ValueError('Grid spacing mismatch')
    counts = [round(args.lx/spacing), round(args.ly/spacing)]
    if min(counts) < 1 or any(abs(count-value/spacing) > 1e-9
                            for count, value in zip(counts, (args.lx, args.ly))):
        raise ValueError('Grid spacing must divide both domain dimensions')
    nx, ny = [v+1 for v in counts]
    summary = {'grid': [nx, ny], 'spacing_m': spacing, 'original_wet_nodes': 0,
               'original_dry_nodes': 0, 'changed_dry_bed_nodes': 0,
               'preserved_wet_bed_nodes': 0, 'maximum_dry_bed_adjustment_m': 0.0,
               'original_wet_cells': 0, 'output_wet_cells': 0}
    with tempfile.TemporaryDirectory(prefix='.repair-geometry-', dir=out.parent) as tmp:
        tmp = Path(tmp)
        previous_old = previous_new = None
        with (tmp/'seabed.in').open('x') as sink:
            sink.write(format(spacing, '.17g')+'\n')
            for i in range(nx):
                row_old, row_new = [], []
                for j in range(ny):
                    s0, _ = old_surface.take()
                    b0, bed_token = old_bed.take()
                    s, surface_token = new_surface.take()
                    if b0-s0 > TOL:
                        raise ValueError('Original bed is above surface at (%d,%d)' % (i,j))
                    wet = s0-b0 > TOL
                    if wet:
                        if abs(s-s0) > TOL or s-b0 <= TOL:
                            raise ValueError('Taper changes original water at (%d,%d)' % (i,j))
                        bed = b0
                        token = bed_token
                        summary['original_wet_nodes'] += 1
                        summary['preserved_wet_bed_nodes'] += 1
                    else:
                        bed = s
                        token = surface_token
                        summary['original_dry_nodes'] += 1
                        summary['changed_dry_bed_nodes'] += bed != b0
                        summary['maximum_dry_bed_adjustment_m'] = max(
                            summary['maximum_dry_bed_adjustment_m'], abs(bed-b0))
                    sink.write(token+'\n')
                    row_old.append((s0, b0))
                    row_new.append((s, bed))
                summary['original_wet_cells'] += validate_cells(
                    previous_old, row_old, i, 'Original', args.topography_base)
                summary['output_wet_cells'] += validate_cells(
                    previous_new, row_new, i, 'Corrected', args.topography_base)
                previous_old, previous_new = row_old, row_new
        summary['input_sha256'] = {str(g.path): g.finish()
                                   for g in (old_surface, old_bed, new_surface)}
        # Copy after validation and check that no concurrent edit changed the
        # surface used to select exact dry tokens.
        copied_hash = hashlib.sha256()
        with args.tapered_surface.open('rb') as source, (tmp/'topography.in').open('xb') as sink:
            for block in iter(lambda: source.read(1024*1024), b''):
                sink.write(block); copied_hash.update(block)
        if copied_hash.hexdigest() != summary['input_sha256'][str(args.tapered_surface)]:
            raise ValueError('Tapered surface changed during repair')
        summary['output_sha256'] = {'topography.in': copied_hash.hexdigest()}
        bed_hash = hashlib.sha256()
        with (tmp/'seabed.in').open('rb') as source:
            for block in iter(lambda: source.read(1024*1024), b''):
                bed_hash.update(block)
        summary['output_sha256']['seabed.in'] = bed_hash.hexdigest()
        summary['invalid_wet_cells'] = 0
        summary['wet_bed_values_changed'] = 0
        (tmp/'repair-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        # mkdir is exclusive, so even a concurrent output creator is protected.
        out.mkdir()
        for name in ('topography.in', 'seabed.in', 'repair-summary.json'):
            shutil.move(str(tmp/name), str(out/name))
    print(json.dumps(summary, indent=2))
    print('Wrote validated pair to '+str(out.resolve()))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        sys.exit('error: ' + str(error))
