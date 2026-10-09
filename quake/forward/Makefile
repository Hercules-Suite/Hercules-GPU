# -*- Makefile -*-
#
#
ifndef WORKDIR
  WORKDIR = $(PWD)/../..
endif

include $(WORKDIR)/systemdef.mk
include $(WORKDIR)/common.mk

# Ensure optimized builds for both CUDA translation units and host C++ code.
# User-provided CFLAGS from user.mk are propagated to C++ compilation and to
# NVCC host compilation. NVFLAGS can still be overridden from user.mk.
CXXFLAGS += $(CFLAGS)
NVFLAGS ?= -O3
NVFLAGS += -Xcompiler "$(CFLAGS)"

# Enable OpenMP by default; users can still override with USE_OPENMP=0.
USE_OPENMP ?= 1

ifeq ($(USE_OPENMP),1)
	CFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

ifeq ($(filter -O%,$(CXXFLAGS)),)
	CXXFLAGS += -O3
endif

# Earthquake related, set these in user.mk
#
# -DDAMPING
# -DDEBUG_OVERFLOW
# -DAUTO_DELTA_T
# -DVALIDATION
# -DHALFSPACE
# -DBOUNDARY
SOLVE_CPPFLAGS ?= -DHALFSPACE -DBOUNDARY


LOADLIBES = $(CVM_DIR)/cvm.o $(OCTOR_LIB) $(ETREE_LIB)

#
# Is vis enabled?
#
ifeq (${ENABLE_VIZ},1)
  VIS_CPPFLAGS ?= -DVISINTERVAL=10
  VIS_CPPFLAGS += -DVIS -I${VIS_DIR}
  LOADLIBES += $(VIS_LIB) $(JPEG_LIB)
endif

# I/O
# -DNO_OUTPUT
# -DUSERCVMDB
IO_CPPFLAGS ?= -DUSECVMDB


CPPFLAGS += -I$(ETREE_DIR) -I$(CVM_DIR) -I$(OCTOR_DIR) -I. -I$(TACC_MPI_DIR)/include
CPPFLAGS += $(SOLVE_CPPFLAGS) $(IO_CPPFLAGS) $(VIS_CPPFLAGS)

LDFLAGS += -lm -L$(TACC_CUDA_DIR)/lib64 -lcudart

ifeq ($(SYSTEM), BGW)
  LOADLIBES += /bgl/local/lib/libmpihpm_c.a /bgl/BlueLight/ppcfloor/bglsys/lib/libbgl_perfctr.rts.a
endif

#
# Object modules
#
NONLINEAR_GPU_OBJECTS = nonlinear/nonlinear_gpu_runtime.o \
	nonlinear/nonlinear_gpu_classic.o \
	nonlinear/nonlinear_gpu_manzari_dafalias.o \
	nonlinear/nonlinear_gpu_pdmy03.o \
	nonlinear/nonlinear_gpu_riva_sand.o \
	nonlinear/nonlinear_gpu_riva_sand_intermediate_bias_research.o \
	nonlinear/nonlinear_gpu_riva_sand_intermediate_bias_material_point.o \
	nonlinear/nonlinear_gpu_material_point.o

PSOLVE_OBJECTS = psolve.o psolve_gpu_launch.o nrutila.o quakesource.o geometrics.o nonlinear.o \
	$(NONLINEAR_GPU_OBJECTS) \
	commutil.o util.o output.o io_planes.o io_drm.o io_checkpoint.o stiffness.o \
	damping.o quake_util.o kernel.o timers.o buildings.o meshformatlab.o drm.o \
	topography.o drm_planewaves.o basin.o


OBJECTS = qmesh.o q4.o single_query.o ${PSOLVE_OBJECTS}

.PHONY: all tags clean cleanall move_cxx_repository

all: move_cxx_repository $(TARGET)

# These translation units share UTMZone_t and the horizontal frame helpers.
psolve.o geometrics.o quakesource.o io_planes.o io_drm.o topography.o \
drm_planewaves.o nonlinear.o: geometrics.h horizontal_frame.h

TARGET = qmesh psolve single_query q4showmeta


qmesh: qmesh.o

psolve: ${PSOLVE_OBJECTS} ${OCTOR_LIB}
	$(LD) $^ \
	$(LOADLIBES) \
	$(LDFLAGS) \
	-o $@

q4showmeta: q4showmeta.o util.o

psolve_output: psolve_output.o

psolve_output.o: psolve.o
	cp -f psolve.o psolve_output.o

psolve_debug: psolve_debug.o

psolve_debug.o: psolve.o nrutila.o quakesource.o geometrics.o io_planes.o io_checkpoint.o
	cp -f psolve.o psolve_debug.o

psolve_novis: psolve_novis.o

psolve_novis.o: psolve.o
	cp -f psolve.o psolve_novis.o

single_query: single_query.o q4.o

output.o: output.cc output.h psolve.h
io_drm.o: io_drm.cc io_drm.h drm_output_points.h drm_output_coordinates.h flat_output_mpi.h psolve.h psolve_gpu_launch.h
io_checkpoint.o: io_checkpoint.cc io_checkpoint.h psolve.h nonlinear_gpu.h
psolve.o: psolve.cc acoustic_interface_setup.h acoustic_interface_trace.h \
          acoustic_open_pore_quadrature.h \
          acoustic_receiver_setup.h acoustic_receiver_gradient.h psolve.h psolve_gpu_launch.h output.h quakesource.h nonlinear.h nonlinear_gpu.h io_planes.h \
          io_checkpoint.h io_drm.h flat_output_mpi.h stiffness.h damping.h quake_util.h timers.h \
		  buildings.h meshformatlab.h drm.h topography.h drm_planewaves.h \
		  basin.h kernel.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

psolve_gpu_launch.o: psolve_gpu_launch.cu psolve_gpu_launch.h psolve.h \
		kernel.h topography.h
	$(NVCC) $(CPPFLAGS) $(NVFLAGS) -diag-suppress=2464 -c $< -o $@

quakesource.o: quakesource.cc quakesource.h psolve.h

nonlinear.o: nonlinear.cc nonlinear.h nonlinear_gpu.h nonlinear_pdmy03.h \
	nonlinear_riva_sand.h nonlinear_riva_sand_intermediate_bias_research.h

stiffness.o: stiffness.cu stiffness.h psolve.h
	$(NVCC) $(CPPFLAGS) $(NVFLAGS) -c $< -o $@

NONLINEAR_GPU_DEPS = nonlinear/nonlinear_gpu_implementation.cuh \
	nonlinear_gpu.h nonlinear.h nonlinear_pdmy03.h psolve.h kernel.h \
	nonlinear_riva_sand.h \
	nonlinear_riva_sand_intermediate_bias_research.h util.h quake_util.h

nonlinear/nonlinear_gpu_%.o: nonlinear/nonlinear_gpu_%.cu $(NONLINEAR_GPU_DEPS)
	$(NVCC) $(CPPFLAGS) $(NVFLAGS) -c $< -o $@

kernel.o: kernel.cu kernel.h psolve.h acoustic_open_pore_quadrature.h
	$(NVCC) $(CPPFLAGS) $(NVFLAGS) -c $< -o $@

all: $(TARGET)

clean:
	rm -f $(OBJECTS)  core *.o *~

cleanall:
	rm -f $(OBJECTS)  $(TARGET) core *~

TAGS:
	etags *.c *.h

tags:	TAGS

# $Id: Makefile,v 1.20 2011/10/07 15:00:10 yisbilir Exp $
