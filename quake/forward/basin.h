#ifndef BASIN_H
#define BASIN_H

#include "cvm.h"

typedef enum basin_material_type_t {
    BASIN_MATERIAL_ELASTIC = 0,
    BASIN_MATERIAL_ACOUSTIC_WATER = 1
} basin_material_type_t;

int getBasinMaterialProperties(cvmpayload_t *g_props, double z_m,
                               basin_material_type_t material_type,
                               double water_vp, double water_rho);

#endif // BASIN_H
