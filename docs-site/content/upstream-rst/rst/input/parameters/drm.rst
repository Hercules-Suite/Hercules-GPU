============================================
Domain Reduction Method (DRM) Implementation
============================================

Hercules is capable of implementing the Domain Reduction Method (DRM). You can implement it by setting 

.. code-block::

    implement_drm = yes

And set parameter ``drm_directory`` to the directory you want to store the DRM information. If you do not set ``drm_directory``, the default setting is ``drm_directory = outputfiles/DRM``. Folders named ``part0``, ``part1``, and ``part2`` will be created automatically in this directory when you run the corresponding part of DRM simulation.

There are three steps should be executed sequentially for implementing DRM in Hercules, which are part 0, part 1, and part 2. You can choose within these three parts by setting the value of ``which_drm_part`` in the input file. For example, setting the following line in the input file for running part 0 in Hercules:

.. code-block::
    
    which_drm_part = part0

Other parameters that are needed to be assigned for DRM will be explained in this documentation. Note that if ``implement_drm = no`` or it is not set (it is an optional parameter), then all the parameters related to DRM in the input file will be neglected, including 

.. code-block::

    drm_directory
    which_drm_part
    drm_edgesize
    part1_delta_t
    part1_simulation_time
    drm_print_rate
    drm_offset_x
    drm_offset_y
    drm_boundary
    DRMBox_Noelem_Halfwidth_EW
    DRMBox_Noelem_Halfwidth_NS
    DRMBox_Noelem_depth
    DRMBox_element_size
    DRM_xo
    DRM_yo

If ``implement_drm = yes``, all the DRM parameters, except ``which_drm_part``, should be the same for all three parts.

.. TODO: it seems other parameters can be changed for different parts. For the 57 scenarios Doriam run (Wenyang said those are run by Doriam, not him), the domain actually shrink to a smaller region in part2 (the DRM box is always the same). Need to check and test this.

Part 0, Part 1, and Part 2
==========================

Part 0
------
In part 0, the information of the coordinates of the DRM box would be generated.

To execute part 0 analysis, simply setting the following line in the input file.

.. code-block::
    
    which_drm_part = part0

Since the information of the coordinates of the DRM box is the only data that are generated, the end time of the simulation can be set as a small value (e.g., 1 sec). You can do so by setting the following line in the input file.

.. code-block::
    
    simulation_end_time_sec = 1

Part 1
------
In part 1, the displacements of DRM nodes would be saved.

To execute part 1 analysis, setting the following line in the input file.

.. code-block::
    
    which_drm_part = part1

Part 2
------
If you are not going to run nonlinear analysis, you don't need to run part 2 unless you want to store effective forces information for local use (e.g., for Abaqus).

In part 2, Hercules no longer get the information from the source. Instead, it gets the information from the results that are calculated in part 0 and part 1.

To execute part 2 analysis, setting the following line in the input file.

.. code-block::
    
    which_drm_part = part2


Other Parameters
================

Defining the Element Size
-------------------------
``drm_edgesize``
    * Optional: Conditional
    * Type: Float
    * Description: The element edge size (in meter) you want to define for elements in the DRM box.

Note that all the dimensions in Hercules have to be "octree dimension." For more information, please refer to :ref:`Octree Dimension`

:math:`\Delta t` and Simulation Time for Part 1
-----------------------------------------------
These parameters only used for part 2 analysis.

You should set the same values for these two parameters, ``part1_delta_t`` and ``part1_simulation_time``, to ``simulation_delta_time_sec`` and ``simulation_end_time_sec``, respectively. You can do so by setting the following lines in the input file.

.. code-block::
    
    part1_delta_t = delta_t_used_in_part1_in_second
    part1_simulation_time = simulation_time_used_in_part1_in_second

``delta_t_used_in_part1_in_second`` and ``simulation_time_used_in_part1_in_second`` mentioned above should be replaced with the :math:`\Delta t` and the simulation end time that are used in part 1.

Print Rate
----------
``drm_print_rate``
    * Optional: Conditional
    * Type: Integer
    * Description: The rate Hercules saves output data for DRM. For example, if you set :code:`drm_print_rate = 2`, it means Hercules would store output data every 2 time step (:math:`\Delta t`).

DRM Offsets
-----------
DRM offsets are the distances (in meter) from the origin that is defined by parameters ``region_origin_latitude_deg`` and ``region_origin_longitude_deg``. You can define these parameters by setting the following parameters.

``drm_offset_x``
    * Optional: Conditional
    * Type: Float
    * Description: The DRM offset (in meter) from the origin in x direction.

``drm_offset_y``
    * Optional: Conditional
    * Type: Float
    * Description: The DRM offset (in meter) from the origin in y direction.

Note that x direction is north-south direction (north is positive), while y direction is east-west direction (east is positive). And both offsets are still :ref:`Octree Dimension`.

DRM offsets create a box bigger than the DRM boundary. The space between these two boundaries is used for the absorbing boundary conditions (ABCs).


DRM Boundary
------------
``drm_boundary``
    * Optional: Conditional
    * Type: Multiple numbers
    * Format: ``min_x min_y max_x max_y depth``
    * Description: Dimensions in DRM boundary section are distances (in meter) from the box that is created with DRM offsets. There are five dimensions needed to be defined., which are ``min_x``, ``min_y``, ``max_x``, ``max_y``, and ``depth``.