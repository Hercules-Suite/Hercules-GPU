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

#ifndef DAMPING_H_
#define DAMPING_H_

/* physics related */

typedef enum
{

  UNSET_DAMP = -1, RAYLEIGH = 0, MASS, NONE, BKT, BKT2, BKT3, BKT3F

} damping_type_t;

/* Absorbing boundary condition selection.
 * ABC_DASHPOT: solid Lysmer dashpots and scalar acoustic radiation damping
 *              (folded into their respective mass coefficients at setup).
 * ABC_CPML:    unsplit Convolutional PML in ADE form (Komatitsch & Martin 2007).
 */
typedef enum
{

  ABC_DASHPOT = 0, ABC_CPML

} abc_type_t;

/* CPML damping-profile sampling along each absorbing normal
 * (cpml_profile_slab). NONE: legacy per-element-centre sampling. AUTO:
 * octree-aligned slabs of the coarsest non-air PML element edge; every
 * element size inside a slab receives identical d/alpha/kappa. */
typedef enum
{

  CPML_PROFILE_SLAB_NONE = 0, CPML_PROFILE_SLAB_AUTO = 1

} cpml_profile_slab_t;

void damping_addforce(mesh_t *myMesh, mysolver_t *mySolver, fmatrix_t (*theK1)[8], fmatrix_t (*theK2)[8]);
void calc_conv(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);
void constant_Q_addforce(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);
void damp_init(int32_t myID, mesh_t *myMesh);

void calc_conv_optimized(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);
void constant_Q_addforce_optimized(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);

void constant_Q_addforce_simplifiedConv_old(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);
void constant_Q_addforce_simplifiedConvolution(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);

void conv_and_bktForceCombined(mesh_t *myMesh, mysolver_t *mySolver, double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping);


void BKT_TU_transf( double theFreq, double theDeltaT, double theDeltaTSquared, damping_type_t typeOfDamping,
        fvector_t *convShear1, fvector_t *convKappa1,
        fvector_t *convShear2, fvector_t *convKappa2,
        fvector_t *convShear3, fvector_t *convKappa3,
        edata_t *edata,
        fvector_t  *tm1Disp, fvector_t *tm2Disp,
        fvector_t  *damping_vector_shear, fvector_t *damping_vector_kappa);


#endif /* DAMPING_H_ */
