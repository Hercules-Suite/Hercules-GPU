=============
Plane Outputs
=============

Once the parameter ``number_output_planes`` is set to a positive integer, Hercules will generate displacement history at grid points defined by ``output_planes``. For the definition of the related parameters, see the section :ref:`Output for Planes`.

In the folder defined by ``output_planes_directory``, Hercules will create files with the name ``planedisplacements.X`` where ``X`` is the number of the plane. ``planedisplacements.X`` is a binary file.

planedisplacements.X
--------------------
Hercules writes the displacement history at the grid points defined by ``output_planes``. The file is a binary file written in double precision. 

The first three numbers are the displacements in the x, y, and z directions at the first grid point (the origin of the defined plane) at the first time step. The next three numbers are the displacements in the x, y, and z directions at the second grid point, which has the same x coordinate and different y coordinate with one ``dy`` away from the first grid point along the y direction, at the first time step. 

It continues to write the displacements at the grid points along the y direction until the last grid point along the y direction. After writing displacements along this first row, Hercules moves on to write displacement along the second row (one ``dx`` away from the first row along the x direction) and so on. 

After writing the displacements at the last grid point, Hercules moves on to the next time step and writes the displacements for all the grid points again.

Similar to the output for mesh coordinates, you can also use NumPy to read the binary file. The following code reads the displacement history at the grid points on the first plane.

.. code-block:: python

    import numpy as np
    with open('planedisplacements.0', 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
        data = data.reshape((number_of_time_steps, number_of_grid_points, 3))

The variable ``number_of_time_steps`` is defined as 

.. math:: 
    \frac{\text{simulation_end_time_sec}}{\text{simulation_delta_time_sec} \times \text{output_planes_print_rate}}

where ``simulation_end_time_sec``, ``simulation_delta_time_sec``, and ``output_planes_print_rate`` are the parameters defined in the parameter input file. The variable ``number_of_grid_points`` is defined as ``nx`` times ``ny``, which are the number of grid points along the x and y directions as defined in ``output_planes``, respectively.