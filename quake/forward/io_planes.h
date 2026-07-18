/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */

#include "psolve.h"
#include "geometrics.h"

extern int planes_print(int32_t myID, int IO_pool_pe_count, int theNumberOfPlanes, 
			mysolver_t* mySolver);
extern void planes_update_acoustic_kinematics(mysolver_t* mySolver,
    int theNumberOfPlanes, int32_t step, double deltaT,
    double deltaTSquared);
extern int32_t planes_displacement_node_count(void);
extern const int32_t* planes_displacement_nodes(void);
extern int32_t planes_acoustic_pressure_node_count(void);
extern const int32_t* planes_acoustic_pressure_nodes(void);
extern int32_t planes_gpu_mpiio_enabled(void);
extern int32_t planes_local_interpolation_point_count(void);
extern int32_t planes_host_only_point_count(void);
extern int planes_export_gpu_interpolation_metadata(
    int32_t pointCount, int32_t* nodes, double* coords, int32_t* modes);
extern int planes_patch_host_only_values(double* packedValues,
                                         int32_t pointCount);
extern int planes_print_gpu_mpiio(int32_t myID, int theNumberOfPlanes,
                                  const double* packedValues,
                                  int32_t pointCount);
extern void planes_setup(int32_t myID, int32_t *thePlanePrintRate, char *thePlaneDirOut, 
    int IO_pool_pe_count, int theNumberOfPlanes, const char *numericalin, 
    double surfaceShift, double *theSurfaceCornersLong, double *theSurfaceCornersLat,
    double theDomainX, double theDomainY, double theDomainZ, char* planes_input_file,
    UTMZone_t* utmZone, int reconstructAcousticPlaneDisplacements);
extern void planes_close(int32_t myID, int IO_pool_pe_count, int theNumberOfPlanes);
extern void planes_IO_PES_main(int32_t myID);
