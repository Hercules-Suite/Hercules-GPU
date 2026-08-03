=======================================
Mesh Coordinates Output for A Subdomain
=======================================

Once the parameter ``mesh_coordinates_for_matlab`` is set to ``yes`` in the parameter input file, Hercules will output the mesh coordinates and the velocity profile (material properties) of each element in a subdomain defined by the parameter ``mesh_corners_matlab``. For the definition of the related parameters, see the section :ref:`Output Mesh Coordinates for A Subdomain`.

In the folder defined by the parameter ``mesh_coordinates_directory_for_matlab``, Hercules will create files with the names ``mesh_coordinates.X`` and ``mesh_data.X`` where ``X`` is the rank number when running Hercules with MPI. For example, if the number of processors when running Hercules is 4, Hercules will create files named ``mesh_coordinates.0``, ``mesh_coordinates.1``, ``mesh_coordinates.2``, and ``mesh_coordinates.3``. The same rule applies to the file ``mesh_data.X``. Both ``mesh_coordinates.X`` and ``mesh_data.X`` are binary files. 

mesh_coordinates.X
------------------
Hercules writes the coordinates of each node in the subdomain sequentially. The first four numbers are the global element ID, the x, y, and z coordinates of the first node of the first element. The next four numbers are the global element ID, the x, y, and z coordinates of the second node of the first element, and so on. Hercules uses 8-node hexahedral elements. The global element ID is written in 64-bit integer, and the coordinates of the nodes are written in double precision.

To read the file ``mesh_coordinates.X`` (``mesh_coordinates.0``, for example), you can use the following Python code:

.. code-block:: python

    import numpy as np
    with open('mesh_coordinates.0', 'rb') as f:
        data = np.fromfile(f, dtype=[('geid', np.int64), ('x', np.float64), ('y', np.float64), ('z', np.float64)])

Then for each entity in the array ``data``, the first number is the global element ID, and the second, third, and fourth numbers are the x, y, and z coordinates of the node, respectively. The first eight entities are the eight nodes of the first element with the same global element ID, the next eight entities are the eight nodes of the second element with another global element ID, and so on.


mesh_data.X
-----------
Hercules writes the velocity profile of each element in the subdomain sequentially. The first four numbers are the global element ID, the :math:`V_s`, :math:`V_p`, and :math:`\rho` of the first element. The next four numbers are the global element ID, the :math:`V_s`, :math:`V_p`, and :math:`\rho` of the second element, and so on. The global element ID is written in 64-bit integer, and the velocity profile of each element is written in single precision.

To read the file ``mesh_data.X`` (``mesh_data.0``, for example), you can use the following Python code:

.. code-block:: python

    import numpy as np
    with open('mesh_data.0', 'rb') as f:
        data = np.fromfile(f, dtype=[('geid', np.int64), ('Vs', np.float32), ('Vp', np.float32), ('rho', np.float32)])

Then for each entity in the array ``data``, the first number is the global element ID, and the second, third, and fourth numbers are the :math:`V_s`, :math:`V_p`, and :math:`\rho` of the element, respectively. The first entity is the velocity profile of the first element, the next entity is the velocity profile of the second element, and so on.


plotmesh.m
----------
In addition to the Python code above, you can also plot the mesh and the velocity profile of each element in the subdomain with the MATLAB script ``plotmesh.m``. The script is located in the Hercules folder under the ``matlab-utils/scripts/`` directory. 

To run the script, you also need another input file. An example input file for this script is provided in the Hercules folder under the ``matlab-utils/examples/`` directory with the name ``plotmeshinput.in``.

Once you set up the input file, the easiest way to use them is to copy the script and the input file to your working directory and run the script as a function in MATLAB with the following command:

.. code-block:: matlab

    plotmesh('plotmeshinput.in')

Note that it is not recommended to plot a large number of elements as it may take a long time to plot them all.