===============
Output Settings
===============

``output_mesh``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should print the output information for the mesh. The printed result will be store in the file specified by ``mesh_etree_output_file``. This is an optional parameter and the default value is False (``output_mesh = 0``).

``output_parallel``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should print the output summary stats. The printed result will be store in the file specified by ``output_stats_file``. This parameter also controls whether Hercules should print the outputs for the displacement and velocity. When ``output_parallel`` is set to False (``output_parallel = 0``), ``output_displacement`` and ``output_velocity`` will be ignored. This is an optional parameter and the default value is False.

``output_displacement``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should print the output information for the displacement. The printed result will be store in the file specified by ``output_displacement_file``. This is an optional parameter and the default value is False (``output_displacement = 0``).

``output_velocity``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should print the output information for the velocity. The printed result will be store in the file specified by ``output_velocity_file``. This is an optional parameter and the default value is False (``output_velocity = 0``).

``simulation_output_rate``
    * Optional: Yes
    * Type: Integer
    * Description: The rate Hercules saves output data for the mesh, the displacement, and the velocity. For example, if you set ``simulation_output_rate = 2``, it means Hercules would store output data every 2 time steps. This is an optional parameter and the default value is ``1000000`` to avoid unexpected and gigantic output files.

``mesh_etree_output_file``
    * Optional: Conditional
    * Type: Path
    * Description: The path to mesh file in etree format. This parameter is required if ``output_mesh`` is set to True (``output_mesh = 1``).

``output_displacement_file``
    * Optional: Conditional
    * Type: Path
    * Description: The path to the displacement output file. This parameter is required if ``output_displacement`` is set to True (``output_displacement = 1``).

``output_velocity_file``
    * Optional: Conditional
    * Type: Path
    * Description: The path to the velocity output file. This parameter is required if ``output_velocity`` is set to True (``output_velocity = 1``).

``output_stats_file``
    * Optional: Yes
    * Type: Path
    * Description: The path to the stats output file. This parameter is only effective when ``output_parallel`` is set to True (``output_parallel = 1``). If ``output_stats_file`` is not set when ``output_parallel = 1``, Hercules will print the output summary stats to the file with the default file name ``output-stats.txt``.

``monitor_file``
    * Optional: Yes
    * Type: Path
    * Description: The path to a text file in which Hercules would write information as it runs for users to monitor the progress of the running job. This is an optional parameter and the default value is ``monitor.txt``.

``stat_mesh_filename``
    * Optional: Yes
    * Type: Path
    * Description: The path to the file where Hercules should store the statistics of the mesh. This is an optional parameter and the default value is ``stat_mesh.txt``.

``stat_schedule_filename``
    * Optional: Yes
    * Type: Path
    * Description: The path to the file where Hercules should store the solver communication schedule summary. This is an optional parameter and the default value is ``stat_sched.txt``.

``stat_topo_filename``
    * Optional: Yes
    * Type: Path
    * Description: The path to the file where Hercules should store the statistics of the topography. This is an optional parameter and the default value is ``stat_topo.txt``.

``print_matrix_k``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether Hercules should print the stiffness matrix to the standard output (usually the command line window running Hercules). This is an optional parameter and the default value is ``no``.


Output Mesh Coordinates for A Subdomain
=======================================
You can ask Hercules to store another mesh coordinates database which includes the information only within a specified subdomain. Here are the parameters you will need to set up if you want.

``mesh_coordinates_for_matlab``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to store an extra mesh coordinates database for the subdomain defined by the parameter ``mesh_corners_matlab``. This is an optional parameter and the default value is ``no``.

``mesh_coordinates_directory_for_matlab``
    * Optional: Yes
    * Type: Path
    * Description: Path to the folder where Hercules should store the mesh coordinates database for the subdomain. This is an optional parameter and the default value is ``outputfiles/For_Matlab``.

``mesh_corners_matlab``
    * Optional: Conditional
    * Type: Multiple Numbers
    * Format: ``min_x min_y max_x max_y min_z max_z``
    * Description: Values used to define the boundaries of the subdomain. Note that X is in the vertical direction of the plan view. This parameter is required if ``mesh_coordinates_for_matlab`` is set to ``yes``.

To better understand the format of the output files for mesh coordinates, refer to the section :ref:`Mesh Coordinates Output for A Subdomain`.


Output for Stations
===================
The following parameters are used for the output for stations (observation points).

``number_output_stations``
    * Optional: Yes
    * Type: Integer
    * Description: The number of output stations. This is an optional parameter and the default value is ``0``.

``output_stations_print_rate``
    * Optional: Conditional
    * Type: Integer
    * Description: The rate Hercules saves output data for stations. For example, if you set ``output_stations_print_rate = 2``, it means Hercules would store output data every 2 time steps. This parameter is required if ``number_output_stations`` is bigger than ``0``.

``print_station_velocities``
    * Optional: Conditional
    * Type: Yes or No
    * Description: Whether to print the velocities at each station. This parameter is required if ``number_output_stations`` is bigger than ``0``.

``print_station_accelerations``
    * Optional: Conditional
    * Type: Yes or No
    * Description: Whether to print the accelerations at each station. This parameter is required if ``number_output_stations`` is bigger than ``0``. Also, if ``print_station_accelerations = yes``, not only the accelerations will be printed, but also the velocities will be printed regardless of the value of ``print_station_velocities``.

``output_stations_directory``
    * Optional: Yes
    * Type: Path
    * Description: Path to the folder where Hercules should store the related information of stations. This is an optional parameter and the default value is ``outputfiles/stations``.

``output_stations``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: ``latitude longitude depth``
    * Description: The location information of stations. The data should be written as many lines as parameter ``number_output_stations`` defines. This parameter is required if ``number_output_stations`` is bigger than ``0``. Note that if the number of listed ``output_stations`` is more than ``number_output_stations``, Hercules will only use the first ``number_output_stations`` stations.

    .. TODO: The difference between engineering depth and bedrock depth needs to be explained.

To better understand the format of the station output files, refer to the section :ref:`Station Outputs`.


Output for Planes
=================
A plane contains multiple grid points storing displacements on them.

``number_output_planes``
    * Optional: Yes
    * Type: Integer
    * Description: The number of output planes. This is an optional parameter and the default value is ``0``.

``planes_input_file``
    * Optional: Yes
    * Type: Path
    * Description: The path to plane input file. The plane information is usually written in the same parameter input file for running Hercules. So, chances are this path is the same as the path to the parameter input file. Nevertheless, it is possible to write the plane information in a separate file. This parameter is an optional parameter and the default value is the same as the path to the parameter input file.

``output_planes_print_rate``
    * Optional: Conditional
    * Type: Integer
    * Description: The rate Hercules saves output data for planes. For example, if you set ``output_planes_print_rate = 2``, it means Hercules would store output data every 2 time steps. This parameter is required if ``number_output_planes`` is bigger than ``0``.

``output_planes_directory``
    * Optional: Conditional
    * Type: Path
    * Description: Path to the folder where Hercules should store the planes database. This is an optional parameter and the default value is ``outputfiles/planes``.

``output_planes``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: ``x_lat y_lon z_depth dx nx dy ny strk dp``
    * Description: Values used to define the planes. ``x_lat``, ``y_lon``, and ``z_depth`` define the origin of a plane. ``dx`` and ``dy`` are the spacing along X and Y directions. ``nx`` and ``ny`` are the number of grid points along X and Y directions. ``strk`` is the strike angle (in degree) measured from the true north clockwise, and ``dp`` is the dip angle (in degree). This parameter is required if ``number_output_planes`` is bigger than ``0``. Note that if the number of listed ``output_planes`` is more than ``number_output_planes``, Hercules will only use the first ``number_output_planes`` planes.

To better understand the format of the plane output files, refer to the section :ref:`Plane Outputs`.


Checkpoint
==========
Hercules can generate checkpoint files during the simulation. Checkpoint files can be used to restart the simulation from the time step when the checkpoint was generated.

``checkpointing_rate``
    * Optional: Yes
    * Type: Integer
    * Description: The rate Hercules saves checkpoints. For example, if you set ``checkpointing_rate = 1000``, it means Hercules would store checkpoint files every 1000 time steps. Hercules will keep saving and overwriting the checkpoint files ``checkpoint.out0`` and ``checkpoint.out1``, which are the last two checkpoint files generated during the simulation, in the folder specified by ``checkpoint_path``. This is an optional parameter and the default value is ``0``, which asks Hercules not to generate checkpoints.

``use_checkpoint``
    * Optional: Yes
    * Type: Boolean
    * Description: Whether Hercules should use the checkpoint file. Note that users should have the checkpoint file ready before setting this parameter to True (``use_checkpoint = 1``). Checkpoint files can be generated by setting ``checkpointing_rate`` to a integer bigger than ``0``, and users have to rename the checkpoint file they want to use to ``checkpoint.in`` and save it in the same folder specified by ``checkpoint_path``. This is an optional parameter and the default value is False (``use_checkpoint = 0``).

``checkpoint_path``
    * Optional: Yes
    * Type: Path
    * Description: Path to the folder where Hercules should store or read the checkpoint files. This is an optional parameter and the default value is ``outputfiles/checkpoints``.
