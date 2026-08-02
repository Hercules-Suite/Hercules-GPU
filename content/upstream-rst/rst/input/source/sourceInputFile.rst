=================
Source Input File
=================

(This section is under construction.)

``source.in``
-------------
The following is a working example of a point source input file:

.. code-block:: 

    source_is_filtered   = 0
    threshold_frequency  = 1.2
    number_of_poles      = 22
    type_of_source       = point

    source_function_type = quadratic
    average_risetime_sec = .5
    moment_amplitude     = 1e15

    # Source Location
    lonlat_or_cartesian  = 0
    hypocenter_lat_deg   = 38.419
    hypocenter_long_deg  = 27.129
    hypocenter_depth_m   = 500
    source_strike_deg    = 0
    source_dip_deg       = 90
    source_rake_deg      = 0


As you can see, the source input file shares the same :ref:`Basic Rules` as the parameter input file.