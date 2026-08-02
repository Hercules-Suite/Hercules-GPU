===============
Station Outputs
===============

When ``number_output_stations`` is set to a positive integer, Hercules will generate responses at the stations (observation points) specified by the parameter ``output_stations``. For the definition of the related parameters, see the section :ref:`Output for Stations`.

In the folder defined by the parameter ``output_stations_directory``, Hercules will create files with the name ``station.X`` where ``X`` is the station number. For example, if ``number_output_stations`` is set to 3 and ``output_stations`` are properly defined, Hercules will create files ``station.0``, ``station.1`` and ``station.2``. Where ``station.0`` contains the responses at the first station defined in ``output_stations`` (the first line after ``output_stations``), ``station.1`` contains the responses at the second station defined in ``output_stations`` and so on. All the ``station.X`` files are ASCII files. Below is an example of a ``station.X`` file with both ``print_station_velocities`` and ``print_station_accelerations`` set to ``yes``:

.. code-block::
    
    #  Time(s)         X|(m)         Y-(m)         Z.(m)       X|(m/s)       Y-(m/s)       Z.(m/s)      X|(m/s2)      Y-(m/s2)      Z.(m/s2)
      0.000000  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
      0.020000  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
      0.040000  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
      0.060000 -2.857821e-22  4.389021e-23 -1.320382e-21 -1.428910e-20  2.194511e-21 -6.601910e-20 -7.144552e-19  1.097255e-19 -3.300955e-18
      0.080000  5.397648e-22 -6.574493e-23  2.056479e-21  4.127734e-20 -5.481757e-21  1.688431e-19  2.778322e-18 -3.838134e-19  1.174311e-17
      0.100000 -1.863683e-21  3.730192e-22 -6.838559e-21 -1.201724e-19  2.193821e-20 -4.447519e-19 -8.072487e-18  1.370998e-18 -3.067975e-17

Note that the time step in the ``station.X`` files is defined by the parameters ``output_stations_print_rate`` and ``simulation_delta_time_sec``.