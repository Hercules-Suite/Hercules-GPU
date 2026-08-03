========
Material
========

.. TEMPLATE
.. ``parameter_name``
..     * Optional:
..     * Type:
..     * Description:

One-Dimensional Layered Profile
===============================

``cvmdb_input_file``
    * Type: Path/profile
    * Description: The path to the material information (crustal structure) database. A special way to assign the crustal structure is setting this parameter to ``profile``, i.e.

    .. code-block::

        cvmdb_input_file = profile

    * By doing so, Hercules will use the crustal structure provided by parameters ``number_profile_layers`` and ``profile_layers`` in the same input file.


``number_profile_layers``
    * Type: Integer
    * Description: How many layers are the parameter ``profile_layers`` defines.


``profile_layers``
    * Type: Multiple numbers
    * Format: ``z_top Vp Vs rho Qp Qs``
    * Description: The crustal structure information defined for a layered half-space. The data should be written as many lines as parameter ``number_profile_layers`` defines. ``z_top`` is the depth (in meter), ``Vp`` and ``Vs`` are velocities for compressional and shear waves (in meter/second), ``rho`` is the density (in kg/meter :superscript:`3`), and ``Qp`` and ``Qs`` are Q factor damping.


Istanbul Velocity Model
=======================

``Istanbul_velocity_model``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to use the velocity model developed for the Istanbul or not. If this parameter is set to ``yes``, Hercules will replace the crustal structure information in the region that is covered by the Istanbul velocity model. For the rest of the region, the material information is still provided by the database or the profile that is assigned with the parameter ``cvmdb_input_file``. This is an optional parameter and the default value is ``no``. 


Basin Velocity Model
====================

``basin_velocity_model``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to use the velocity model developed for a basin. If this parameter is set to ``yes``, Hercules will replace the crustal structure information within the basin. For the rest of the region, the material information is still provided by the database or the profile that is assigned with the parameter ``cvmdb_input_file``. This is an optional parameter and the default value is ``no``. Note that only one of the ``Istanbul_velocity_model``, ``basin_velocity_model``, and ``3D_velocity_model`` can be set to ``yes``. If multiple options are set to ``yes``, the system prioritizes the models in the following order: Istanbul model, Basin model, and then the 3D velocity model. The latter models won't be used if the former models are used.


``basin_input_file``
    * Optional: Conditional
    * Type: Path
    * Description: The path to the basin input file. This parameter is required when ``basin_velocity_model`` is set to ``yes``. The basin input file is a text file that describes the boundary of the basin. The format of the basin input file is described in the section :ref:`Basin Input File`.


``basin_origin_latitude_deg``
    * Optional: Conditional
    * Type: Float
    * Description: The latitude of the basin origin (in degree). Note that the origin is the left bottom corner of the domain of the basin. This parameter is required when ``basin_velocity_model`` is set to ``yes``.


``basin_origin_longitude_deg``
    * Optional: Conditional
    * Type: Float
    * Description: The longitude of the basin origin (in degree). Note that the origin is the left bottom corner of the domain of the basin. This parameter is required when ``basin_velocity_model`` is set to ``yes``.


Basin Input File
----------------

The basin input file should be a text file (the extension is not important) that describes the boundary of the basin. The format of the basin input file is described below.

    * The first line: number_of_x_points number_of_y_points number_of_z_points

Note the the number of z points should always be 1.

    * The second and the following lines: x y z
    
Repeat the line with x points written in ascending order, then y points in ascending order, and z in ascending. The following is a working basin input file example:

.. code-block::

    6 2 1
    0.0 0.0 0.0
    50.0 0.0 0.0
    100.0 0.0 3.0
    150.0 0.0 4.5
    200.0 0.0 5.0
    250.0 0.0 5.0
    0.0 50.0 0.0
    50.0 50.0 0.0
    100.0 50.0 2.8
    150.0 50.0 4.2
    200.0 50.0 4.8
    250.0 50.0 5.0

Note that the increments of x and y should be the same. The increment of z can be different. The unit of x, y, and z is meter.

The total number of lines in the basin input file should be ``number_of_x_points * number_of_y_points + 1``. In the example above, the total number of lines is ``13``.


basin.c
-------

Unfortunately, we don't have a universal way to read a text file to define the equations of the material properties used in the basin. Therefore, users have to modify the file ``basin.c`` to define these equations. The file ``basin.c`` is located in ``./quake/forward`` and contains a function ``getBasinMaterialProperties`` that defines these equations. The following is an example of ``basin.c``:

.. code-block:: c

    #include <math.h>
    #include "basin.h"

    int getBasinMaterialProperties(cvmpayload_t *g_props, double z_m) {
        g_props->rho = 2140.0 + 0.125*z_m;
        g_props->Vp = 1000.0 + 1.2*z_m;        
        g_props->Vs = 320.0 + 19.0*sqrt(z_m);
        g_props->Qs = 0.4*g_props->Vs;
        g_props->Qp = 0.5*g_props->Qs;
        return 0;
    }

Note that you should re-compile Hercules after modifying ``basin.c``.


3D Velocity Model
=================
Contrast to the Istanbul Velocity Model, the 3D Velocity Model is a general velocity model that covers any interested regions. Users need to provide the velocity model in the designated format. The format of the 3D Velocity Model is described in the section :ref:`3D Velocity Model Input Files`.

``3D_velocity_model``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to use the general 3D velocity model or not. If this parameter is set to ``yes``, Hercules will replace the crustal structure information in the region that is covered by the general 3D velocity model. For the rest of the region, the material information is still provided by the database or the profile that is assigned with the parameter ``cvmdb_input_file``. This is an optional parameter and the default value is ``no``. Note that only one of the ``Istanbul_velocity_model``, ``basin_velocity_model``, and ``3D_velocity_model`` can be set to ``yes``. If multiple options are set to ``yes``, the system prioritizes the models in the following order: Istanbul model, Basin model, and then the 3D velocity model. The latter models won't be used if the former models are used.

``3D_velocity_model_origin_latitude_deg``
    * Optional: Conditional
    * Type: Float
    * Description: The latitude of the origin of the 3D velocity model (in degree). Note that the origin is the left bottom corner of the domain of the 3D velocity model. This parameter is required when ``3D_velocity_model`` is set to ``yes``.

``3D_velocity_model_origin_longitude_deg``
    * Optional: Conditional
    * Type: Float
    * Description: The longitude of the origin of the 3D velocity model (in degree). Note that the origin is the left bottom corner of the domain of the 3D velocity model. This parameter is required when ``3D_velocity_model`` is set to ``yes``.


Misc
====
``softening_factor``
    * Optional: Yes
    * Type: Float
    * Description: The softening factor used to adjust material properties according to the element size. A factor of 1 means perfect compliance between the mesh and the elements' material properties resulting in strong changes to the results. A factor of 4 tends to double the simulation :math:`\Delta t` without affecting too much the results. Testing is needed, and it's recommended to use a factor larger than 4. The possible value for this parameter is 0 and any floating number larger than 1. This is an optional parameter and the default value is ``0``, which asks Hercules not to change material properties with softening factor.

``the_threshold_Vp_over_Vs``
    * Type: Float
    * Description: The threshold to limit the :math:`V_p/V_s` ratio. If the computed :math:`V_p/V_s` is larger than this threshold, :math:`V_p` will be set to this threshold multiplied by :math:`V_s` when correcting material properties for mesh elements and in the calculation of the first Lamé parameter :math:`\lambda`.
