=======
Running
=======

Basic Usage
===========
To run Hercules, the basic usage is

.. code-block:: shell

    psolve inputFile.in

The ``inputFile.in`` mentioned above should be changed to your actual input file (default: ``parameters.in``). The :code:`psolve` executable file is located in ``quake/forward/`` under the Hercules folder. Assuming your working directory is already changed to the Hercules folder, then the proper way to call :code:`psolve` is

.. code-block:: shell

    quake/forward/psolve inputFile.in

Alternatively, you could also add the ``quake/forward/`` to your system ``$PATH``.

Running with Open-MPI
=====================
Now you should have Open-MPI installed on your machine, you can utilize the maximum power of your machine with the use of it. The basic usage with Open-MPI is

.. code-block:: shell

    mpirun -n numProcesses psolve inputFile.in

``numProcesses`` mentioned above should be changed to the number of processes (cores or threads, depending on your system). For example, assuming we want to run Hercules with 4 cores, and the current working directory is the Hercules folder, then the proper command is

.. code-block:: shell

    mpirun -n 4 quake/forward/psolve inputFile.in

