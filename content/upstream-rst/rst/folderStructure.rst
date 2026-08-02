================
Folder Structure
================

Although many input and output paths in Hercules are configurable, the suggested (default) folder structure is as follows:

.. code-block::

    working_directory/
    ├── inputfiles/
    │   ├── materialfiles/
    │   │   └── ...
    │   ├── sourcefiles/
    │   │   ├── source.in
    │   │   └── ...
    │   ├── topography/
    │   │   └── topography.in
    │   └── parameters.in
    ├── outputfiles/
    │   ├── checkpoints/
    │   │   ├── checkpoint.out0
    │   │   └── checkpoint.out1
    │   ├── DRM/
    │   │   ├── part0/
    │   │   ├── part1/
    │   │   └── part2/
    │   ├── For_Matlab/
    │   │   ├── mesh_coordinates.0
    │   │   ├── mesh_coordinates.1
    │   │   ├── ...
    │   │   ├── mesh_data.0
    │   │   ├── mesh_data.1
    │   │   └── ...
    │   ├── planes/
    │   │   ├── planedisplacements.0
    │   │   ├── planedisplacements.1
    │   │   └── ...
    │   ├── srctmp/
    │   └── stations/
    │       ├── station.0
    │       ├── station.1
    │       └── ...
    ├── monitor.txt
    ├── stat-mesh.txt
    ├── stat-sched.txt
    └── stat-topo.txt


Folder Structure for Input Files
================================

All files under the ``inputfiles/`` folder should be prepared by the user. The ``parameters.in`` and ``source.in`` files are mandatory, while the other files are optional depending on the simulation setup. This section focuses on the structure of the ``inputfiles/`` folder. For more information on the format of these files, see the :ref:`Input Files` section.

The ``materialfiles/`` folder contains the input files for the 3D Velocity Model. For 1D layered soil profile, you will define the profile in the parameter input file, and you don't need this folder.

The structure of this folder is shown below. For more information on the format of these files, see the :ref:`3D Velocity Model Input Files` section.

.. code-block::

    materialfiles/
    ├── index.in
    ├── rho.in
    ├── vp.in
    ├── vs.in
    ├── x.in
    ├── y.in
    └── z.in


The ``sourcefiles/`` folder contains the input files for the source. The structure of this folder is shown below. The ``source.in`` is a mandatory file, while other files are optional depends on the source type.

.. code-block::

    sourcefiles/
    ├── source.in
    ├── area.in
    ├── coords.in
    ├── dip.in
    ├── rake.in
    ├── slip.in
    ├── slipfunction.in
    └── strike.in
    

The only needed file in the ``topography/`` folder is ``topography.in`` file. This file contains the topography data. For more information, see the :ref:`Topography Input File` section.


Folder Structure for Output Files
=================================
All the folders and files used to store the output will be created by Hercules automatically if they don't exist, including ``monitor.txt``, ``stat-mesh.txt``, ``stat-sched.txt``, and ``stat-topo.txt``, and everything under the ``outputfiles/`` folder. For more information on the format of these files, see the :ref:`Output Files` section.