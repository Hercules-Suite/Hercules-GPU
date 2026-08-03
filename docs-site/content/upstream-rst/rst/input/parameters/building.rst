========
Building
========

.. TODO: some details should be added for some parameters.

``include_buildings``
    * Optional: Yes
    * Type: Yes or No
    * Description: Whether to consider buildings in the simulation. This is an optional parameter and the default value is ``no``. 

``number_of_buildings``
    * Optional: Conditional
    * Type: Integer
    * Description: How many buildings are the parameter ``building_properties`` defines. This parameter is required if ``include_buildings`` is set to ``yes``.

``min_octant_size_m``
    * Optional: Conditional
    * Type: Float
    * Description: The minimum octant size for buildings. This parameter is required if ``include_buildings`` is set to ``yes``.

.. TODO: is this surface shift upward or downward?

``surface_shift_m``
    * Optional: Conditional
    * Type: Float
    * Description: The surface shift for buildings. It can be any floating number not less than 0. This parameter is required if ``include_buildings`` is set to ``yes``.

.. TODO: The explanations for all the fixed-base related parameters.

``consider_fixed_base``
    * Optional: Conditional
    * Type: Yes or No
    * Description: This parameter is required if ``include_buildings`` is set to ``yes``.

``fixedbase_input_dt``
    * Optional: Conditional
    * Type: Float
    * Description: This parameter is required if ``consider_fixed_base`` is set to ``yes``.

``fixedbase_input_dir``
    * Optional: Conditional
    * Type: Path
    * Description: The directory where the fixed-base input files are located. This parameter is required if ``consider_fixed_base`` is set to ``yes``.

``fixedbase_input_sufix``
    * Optional: Conditional
    * Type: String
    * Description: The file name of the fixed-base inputs. This parameter is required if ``consider_fixed_base`` is set to ``yes``.

``fixedbase_input_startindex``
    * Optional: Conditional
    * Type: Integer
    * Description: The starting index of the fixed-base input files. For example, if ``fixedbase_input_sufix = station`` and ``fixedbase_input_startindex = 448``, then Hercules will read the fixed-base input files as ``station.448``, ``station.449``, ``station.450``, etc. from the directory ``fixedbase_input_dir``. This parameter is required if ``consider_fixed_base`` is set to ``yes``.

``building_properties``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: ``min_x max_x min_y max_y depth height Vp_building Vs_building Rho_building Vp_foundation Vs_foundation Rho_foundation``
    * Description: Material properties for buildings. The data should be written as many lines as parameter ``number_of_buildings`` defines. This parameter is required if ``number_of_buildings`` is bigger than ``0``. Note that if the number of listed ``building_properties`` is more than ``number_of_buildings``, Hercules will only use the first ``number_of_buildings`` stations.
