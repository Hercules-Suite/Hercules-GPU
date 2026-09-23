#ifndef IO_DRM_H
#define IO_DRM_H
#include "psolve.h"
#include "geometrics.h"

void drm_output_setup(int32_t rank, const char* parameters, double surface_shift,
                      int topography, double* corner_lon, double* corner_lat,
                      double domain_x, double domain_y, double domain_z,
                      UTMZone_t* zone, mysolver_t* solver, double delta_t);
void drm_output_print(mysolver_t* solver, int32_t step);
void drm_output_close();
#endif
