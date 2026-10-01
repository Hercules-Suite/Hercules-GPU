#!/usr/bin/env python3
"""Check surface enrollment, finite/nonzero fields and station/plane parity.

This is a CUDA/MPI smoke test, not an analytical wavefield acceptance test.
Python 3 and NumPy are required. Run from the copied case directory.
"""
from pathlib import Path
import re
import numpy as np

log = Path('output/monitor.txt').read_text()
surface = re.search(r'Acoustic free-surface geometry: solver Z=([\d.e+-]+) m, (\d+) exposed', log)
assert surface and float(surface[1]) == 250 and int(surface[2]) > 0, 'shifted free surface missing'
if 'gravity (linearized' in log:
    pml = re.search(r'CPML gravity free-surface correction: (\d+) top-face', log)
    assert pml and int(pml[1]) > 0, 'gravity/CPML surface correction missing'
output = re.search(r'Acoustic output CPML gradient correction: (\d+) stations, (\d+) plane points', log)
assert output and int(output[1]) > 0 and int(output[2]) > 0, 'output CPML correction missing'
plane = np.fromfile('output/planes/planedisplacements.0', dtype='=f8').reshape(800,11,11,3)
probe = np.fromfile('output/planes/planedisplacements.1', dtype='=f8').reshape(800,3)
station = np.loadtxt('output/stations/station.0')
assert np.isfinite(plane).all() and np.isfinite(station).all(), 'nonfinite output'
assert np.max(np.abs(plane)) > 1e-12, 'no excited plane motion'
assert np.max(np.abs(plane)) < 1.0, 'gross wavefield growth in this weak-source fixture'
assert np.max(np.abs(station[:,1:4])) < 1.0, 'gross station displacement growth'
assert station.shape[0] == 800, 'incomplete station output'
np.testing.assert_allclose(station[:,0], np.arange(800)*.01, atol=1e-8, rtol=0)
# Both observations use the same latitude/longitude at x=0, y~500 m, depth=25 m.
# Unlike a doubly clamped corner, the normal gradient here can be nonzero.
assert np.max(np.abs(probe)) > 1e-12, 'no excited CPML probe motion'
np.testing.assert_allclose(station[:,1:4], probe, atol=1e-10, rtol=1e-5)
print('PASS: shifted surface and CPML enrolled; finite excited fields; station/plane parity')
print('Plane peak XYZ (m):', np.max(np.abs(plane), axis=(0,1,2)))
