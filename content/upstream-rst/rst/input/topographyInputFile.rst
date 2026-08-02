=====================
Topography Input File
=====================
The only file needed in ``topography_directory`` is ``topography.in``. This file contains the elevations across the whole simulation domain. The first line is the spacing between the points in both x and y directions. For the rest of lines, each line contains the elevation of the corresponding point in the domain. 

The first point is the bottom left corner of the domain and the points are ordered from south to north (x direction in Hercules) and from west to east (y direction in Hercules). In other words, you should finish writing elevations from all points with the same y coordinate before moving to the next y coordinate. The elevation should be in meters.