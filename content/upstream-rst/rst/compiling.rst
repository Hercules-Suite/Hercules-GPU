=========
Compiling
=========

Prerequisites
=============
To compile Hercules, you need to have Open-MPI and GNU Scientific Library (GSL) installed. If you haven't installed them, the easiest way to do so is installing them via `Homebrew`_. Assuming Homebrew is installed on your machine, all you need to do is typing the following commands in the terminal one by one.

.. code-block:: shell

    brew install open-mpi
    brew install gsl

.. _Homebrew: https://brew.sh

Compiling with ``user.mk``
==========================
There are some predefined system configurations in ``systemdef.mk``. On some well-defined systems, you only need to create a ``user.mk`` as simple as

.. code-block:: 

    SYSTEM = Frontera

Replace ``Frontera`` with the system name defined in ``systemdef.mk`` that you are using, place ``user.mk`` next to ``systemdef.mk`` in the same directory, and you are good to go. However, it's very likely you will find that it's insufficient, especially when you are building Hercules on your own machine. 

The following is an example of a more complex ``user.mk`` file used for Intel-based Macs. Your might need to modify it according to your system configuration, but it should give you a good starting point.

.. code-block::

    # System definition
    # -----------------

    SYSTEM = MAC

    # Main compilation setup
    # ----------------------

    CFLAGS += -g -ggdb -O3

    # Compilation flags
    # -----------------

    # CVM and IO related:

    IO_CPPFLAGS = -DUSECVMDB -DSCEC  -DPROCPERNODE=4

    ifeq ($(SYSTEM), MAC)
        MPI_DIR      = /usr/local
        MPI_INCLUDE  = $(MPI_DIR)/include/openmpi/ompi/mpi/cxx
        CC           = $(MPI_DIR)/bin/mpicc
        CXX          = $(MPI_DIR)/bin/mpicxx
        LD           = $(MPI_DIR)/bin/mpicxx
        CXXFLAGS    += -DMPICH_IGNORE_CXX_SEEK
        CFLAGS      += -Wall
        CPPFLAGS    += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
        LDFLAGS     += $(MPI_DIR)/lib/libgsl.a
    endif

    LDFLAGS += -lgslcblas -lgsl

Note that the indentation used in a makefile should be tabs, not spaces. If you copy the code block above, make sure to replace all spaces with tabs.

PROJ Support
============
Hercules supports `PROJ`_. The easiest way to install PROJ is via Homebrew by running the following command in the terminal.

.. code-block:: shell

    brew install proj

If you want to use PROJ in Hercules, you must add the following lines to your ``user.mk``.

.. code-block::

    CPPFLAGS += -DPROJ # Define the PROJ macro
    CFLAGS += -I/opt/homebrew/Cellar/proj/9.4.1/include
    LDFLAGS += -L/opt/homebrew/Cellar/proj/9.4.1/lib -lproj

The path to the PROJ library and include directory might be different on your machine. Make sure to replace them with the correct path.

Although PROJ is optional, it is highly recommended to compile Hercules with PROJ support to get a more accurate simulation result.

.. _PROJ: https://proj.org/

Compiling Procedure
===================
To compile Hercules, make sure to :code:`cd` into Hercules' folder and simply run the following lines in the terminal one by one.

.. code-block::

    make cleanall
    make

Now, you are all set.

Side Notes for Apple Silicon Macs
=================================
If you are installing Hercules on ARM-based Macs, which are equipped with M-series SoC (the Apple Silicon), here are some notes for you to install Hercules successfully.

* Install Open-MPI and GSL with native Apple Silicon support using natively supported Homebrew. The installed Open-MPI and GSL should be located under :code:`/opt/homebrew/`
* Change ``user.mk`` accordingly
* (Optional) The last step of compiling procedure (:code:`make`) might need to be replaced with :code:`arch -arm64e make`

``user.mk`` for Apple Silicon Macs can be written as follows.

.. code-block::

    # System definition
    # -----------------

    SYSTEM = ARMMAC

    # Main compilation setup
    # ----------------------

    CFLAGS += -g -ggdb -O3

    # Compilation flags
    # -----------------

    # CVM and IO related:

    IO_CPPFLAGS = -DUSECVMDB -DSCEC  -DPROCPERNODE=4

    ifeq ($(SYSTEM), ARMMAC)
        MPI_DIR      = /opt/homebrew
        MPI_INCLUDE  = $(MPI_DIR)/include/
        CC           = $(MPI_DIR)/bin/mpicc
        CXX          = $(MPI_DIR)/bin/mpicxx
        LD           = $(MPI_DIR)/bin/mpicxx
        CXXFLAGS    += -DMPICH_IGNORE_CXX_SEEK
        CFLAGS      += -Wall -I$(MPI_DIR)/include/
        CPPFLAGS    += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
        LDFLAGS     += -L$(MPI_DIR)/lib/
    endif

    LDFLAGS += -lgslcblas -lgsl

And yes, you can combine this section with the ``user.mk`` for Intel-based Macs. Just make sure you set the correct :code:`SYSTEM` variable. Alternatively, you can use :code:`uname -s` to determine the system type and :code:`uname -m` to determine the system architecture, and set the :code:`SYSTEM` and :code:`ARCH` variables accordingly. Here is an example:

.. code-block::

    # System definition
    # -----------------

    SYSTEM = $(shell uname -s | tr A-Z a-z)
    ARCH   = $(shell uname -m | tr A-Z a-z)

    # Main compilation setup
    # ----------------------

    CFLAGS += -g -ggdb -O3

    # Compilation flags
    # -----------------

    # CVM and IO related:

    IO_CPPFLAGS = -DUSECVMDB -DSCEC  -DPROCPERNODE=4

    ifeq ($(SYSTEM), darwin)
        ifeq ($(ARCH), arm64)
            MPI_DIR      = /opt/homebrew
            MPI_INCLUDE  = $(MPI_DIR)/include/
            CC           = $(MPI_DIR)/bin/mpicc
            CXX          = $(MPI_DIR)/bin/mpicxx
            LD           = $(MPI_DIR)/bin/mpicxx
            CXXFLAGS    += -DMPICH_IGNORE_CXX_SEEK
            CFLAGS      += -Wall -I$(MPI_DIR)/include/
            CPPFLAGS    += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
            LDFLAGS     += -L$(MPI_DIR)/lib/
        else
            MPI_DIR      = /usr/local
            MPI_INCLUDE  = $(MPI_DIR)/include/openmpi/ompi/mpi/cxx
            CC           = $(MPI_DIR)/bin/mpicc
            CXX          = $(MPI_DIR)/bin/mpicxx
            LD           = $(MPI_DIR)/bin/mpicxx
            CXXFLAGS    += -DMPICH_IGNORE_CXX_SEEK
            CFLAGS      += -Wall
            CPPFLAGS    += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
            LDFLAGS     += $(MPI_DIR)/lib/libgsl.a
        endif
    endif

    LDFLAGS += -lgslcblas -lgsl