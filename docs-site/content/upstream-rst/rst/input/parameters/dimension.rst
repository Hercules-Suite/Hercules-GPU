=========
Dimension
=========

Defining the Domain
===================

Domain Origin
-------------
To define the domain origin, the following parameter should be set.

``region_origin_latitude_deg``
    * Type: Float
    * Description: The latitude (in degree) of the origin. Note that the origin is the left bottom corner of the domain. 

``region_origin_longitude_deg``
    * Type: Float
    * Description: The longitude (in degree) of the origin. Note that the origin is the left bottom corner of the domain. 

``region_depth_shallow_m``
    * Type: Float
    * Description: The starting depth (in meter) of the region.

Note that for DRM :ref:`Part 2` analysis, the origin here is the left bottom corner of the box that is created by :ref:`DRM Offsets`.

Domain Dimension
----------------
To define the domain dimension, the following parameter should be set.

``region_length_east_m``
    * Type: Float
    * Description: The length (in meter) from the origin to the local east boundary.

``region_length_north_m``
    * Type: Float
    * Description: The length (in meter) from the origin to the local north boundary.

``region_depth_deep_m``
    * Type: Float
    * Description: The length (in meter) from the origin to the depth limit.

Note that all these dimensions have to be :ref:`Octree Dimension` and integer.

.. TODO: Have to investigate why they have to be integers (if they are not, the program will use a ton of memory and not running at all).

Domain Surface Corners
======================
``domain_surface_corners``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: ``longitude latitude``
    * Description: The coordinates (latitudes and longitudes) of four corners of the domain. This parameter is required when Hercules is not compiled with :ref:`PROJ Support`. The data should be written as many as four lines. The first line contains the same information as parameters ``region_origin_longitude_deg`` and ``region_origin_latitude_deg`` define. And the rest of three corners should be defined clockwise.


DRM Dimensions
==============
For DRM (Domain Reduction Method) dimensions, check :ref:`DRM Offsets` and :ref:`DRM Boundary`.