#!/usr/bin/env python3
"""Extract production station routines for the host/MPI regression harness.

This compiles the actual interpolation and output dispatch without the CUDA
solver. Tensor output and GPU transfers are outside this harness's scope.
"""
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
psolve = (root / 'quake/forward/psolve.cc').read_text()
topography = (root / 'quake/forward/topography.cc').read_text()


def function(source, signature):
    start = source.index(signature)
    # These routines have balanced braces in comments/strings as well as code.
    brace = source.index('{', start)
    level = 1
    end = brace + 1
    while level:
        level += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[start:end] + '\n'


Path(sys.argv[1]).write_text(
    function(topography, 'int compute_tetra_displ (') +
    function(psolve, 'static int interpolate_station_displacements(int32_t step) {') +
    function(psolve, 'static void solver_output_stations(int step) {') +
    function(psolve, 'static void solver_close_station_outputs(void) {')
)
