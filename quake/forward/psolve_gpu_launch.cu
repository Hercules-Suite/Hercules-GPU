/* -*- CUDA -*- */

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
 * @copyright_notice_end */

#include "psolve_gpu_launch.h"
#include "kernel.h"

__global__ static void kernelSourceForceScatter(int32_t count,
                                                const int32_t* lnidList,
                                                const vector3D_t* srcForces,
                                                fvector_t* forceDevice,
                                                double dt2) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  int32_t lnid = lnidList[idx];
  const vector3D_t src = srcForces[idx];
  fvector_t out;
  out.f[0] = src.x[0] * dt2;
  out.f[1] = src.x[1] * dt2;
  out.f[2] = src.x[2] * dt2;
  forceDevice[lnid] = out;
}

__global__ static void kernelSparseForceScatter(int32_t count,
                                                const int32_t* lnidList,
                                                const fvector_t* srcForces,
                                                fvector_t* forceDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  int32_t lnid = lnidList[idx];
  forceDevice[lnid] = srcForces[idx];
}

__global__ static void kernelDispBoundaryGather(int32_t count,
                                                const int32_t* lnidList,
                                                const fvector_t* dispDevice,
                                                fvector_t* bufferDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  bufferDevice[idx] = dispDevice[lnidList[idx]];
}

__global__ static void kernelDispBoundaryScatter(int32_t count,
                                                 const int32_t* lnidList,
                                                 const fvector_t* bufferDevice,
                                                 fvector_t* dispDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  dispDevice[lnidList[idx]] = bufferDevice[idx];
}

__global__ static void kernelScalarBoundaryGather(int32_t count,
                                                  const int32_t* lnidList,
                                                  const solver_float* srcDevice,
                                                  solver_float* bufferDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  bufferDevice[idx] = srcDevice[lnidList[idx]];
}

__global__ static void kernelScalarBoundaryScatter(int32_t count,
                                                   const int32_t* lnidList,
                                                   const solver_float* bufferDevice,
                                                   solver_float* dstDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  dstDevice[lnidList[idx]] = bufferDevice[idx];
}

__global__ static void kernelScalarBoundaryScatterSerial(
    int32_t count, const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice) {
  if ((blockIdx.x != 0) || (threadIdx.x != 0)) {
    return;
  }
  for (int32_t idx = 0; idx < count; idx++) {
    dstDevice[lnidList[idx]] = bufferDevice[idx];
  }
}

__global__ static void kernelScalarBoundaryAccumulate(
    int32_t count, const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  atomicAdd(dstDevice + lnidList[idx], bufferDevice[idx]);
}

/* The host compute_adjust cache is already grouped by parent for the
 * contribution path. One GPU thread therefore owns each parent update. */
__global__ static void kernelScalarAdjustDistribution(
    int32_t runCount, const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values) {
  int32_t run = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (run >= runCount) {
    return;
  }

  solver_float sum = 0.0;
  for (int32_t edge = runOffsets[run]; edge < runOffsets[run + 1]; edge++) {
    const int32_t dnindex = edgeDn[edge];
    sum += values[danglingBase3[dnindex] / 3] * invDeps[dnindex];
  }
  values[runParentBase3[run] / 3] += sum;
}

__global__ static void kernelScalarAdjustAssignment(
    int32_t dnCount, const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values) {
  int32_t orderIndex = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (orderIndex >= dnCount) {
    return;
  }

  const int32_t dnindex = dnOrder[orderIndex];
  solver_float sum = 0.0;
  for (int32_t parent = parentOffsets[dnindex];
       parent < parentOffsets[dnindex + 1]; parent++) {
    sum += values[parentBase3[parent] / 3];
  }
  values[danglingBase3[dnindex] / 3] = sum * invDeps[dnindex];
}

__global__ static void kernelVectorBoundaryGather(int32_t count,
                                                  const int32_t* lnidList,
                                                  const fvector_t* srcDevice,
                                                  fvector_t* bufferDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx < count) {
    bufferDevice[idx] = srcDevice[lnidList[idx]];
  }
}

__global__ static void kernelVectorBoundaryScatter(int32_t count,
                                                   const int32_t* lnidList,
                                                   const fvector_t* bufferDevice,
                                                   fvector_t* dstDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx < count) {
    dstDevice[lnidList[idx]] = bufferDevice[idx];
  }
}

__global__ static void kernelVectorBoundaryScatterSerial(
    int32_t count, const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice) {
  if ((blockIdx.x != 0) || (threadIdx.x != 0)) {
    return;
  }
  for (int32_t idx = 0; idx < count; idx++) {
    dstDevice[lnidList[idx]] = bufferDevice[idx];
  }
}

__global__ static void kernelVectorBoundaryAccumulate(
    int32_t count, const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= count) {
    return;
  }

  const int32_t node = lnidList[idx];
  atomicAdd(&dstDevice[node].f[0], bufferDevice[idx].f[0]);
  atomicAdd(&dstDevice[node].f[1], bufferDevice[idx].f[1]);
  atomicAdd(&dstDevice[node].f[2], bufferDevice[idx].f[2]);
}

__global__ static void kernelVectorAdjustDistribution(
    int32_t runCount, const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values) {
  int32_t run = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (run >= runCount) {
    return;
  }

  fvector_t sum = {{0.0, 0.0, 0.0}};
  for (int32_t edge = runOffsets[run]; edge < runOffsets[run + 1]; edge++) {
    const int32_t dnindex = edgeDn[edge];
    const fvector_t value = values[danglingBase3[dnindex] / 3];
    const solver_float weight = invDeps[dnindex];
    sum.f[0] += value.f[0] * weight;
    sum.f[1] += value.f[1] * weight;
    sum.f[2] += value.f[2] * weight;
  }
  const int32_t parent = runParentBase3[run] / 3;
  values[parent].f[0] += sum.f[0];
  values[parent].f[1] += sum.f[1];
  values[parent].f[2] += sum.f[2];
}

__global__ static void kernelVectorAdjustAssignment(
    int32_t dnCount, const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values) {
  int32_t orderIndex = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (orderIndex >= dnCount) {
    return;
  }

  const int32_t dnindex = dnOrder[orderIndex];
  fvector_t sum = {{0.0, 0.0, 0.0}};
  for (int32_t parent = parentOffsets[dnindex];
       parent < parentOffsets[dnindex + 1]; parent++) {
    const fvector_t value = values[parentBase3[parent] / 3];
    sum.f[0] += value.f[0];
    sum.f[1] += value.f[1];
    sum.f[2] += value.f[2];
  }
  const solver_float weight = invDeps[dnindex];
  fvector_t* dst = &values[danglingBase3[dnindex] / 3];
  dst->f[0] = sum.f[0] * weight;
  dst->f[1] = sum.f[1] * weight;
  dst->f[2] = sum.f[2] * weight;
}

__device__ __constant__ int kTopoTetNodes[2][5][4] = {
    {{0, 2, 1, 4}, {3, 1, 2, 7}, {6, 4, 7, 2},
     {5, 7, 4, 1}, {2, 4, 7, 1}},
    {{0, 3, 1, 5}, {0, 2, 3, 6}, {4, 5, 6, 0},
     {6, 5, 7, 3}, {0, 6, 3, 5}}};

__device__ __constant__ double kTopoCubeCoord[8][3] = {
    {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
    {1.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
    {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

__device__ static double topo_det3(const double a[3][3]) {
  return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
         a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
         a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

__device__ static void topo_inv3(const double a[3][3], double inv[3][3]) {
  double d = topo_det3(a);

  inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / d;
  inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / d;
  inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / d;
  inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / d;
  inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / d;
  inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / d;
  inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / d;
  inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / d;
  inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / d;
}

__device__ static void topoTetraForcesEffective(const fvector_t* un,
                                                fvector_t* resVec,
                                                const topoconstants_t* topo_ec) {
  const int part = (topo_ec->cube_part == 1) ? 0 : 1;
  const double h = topo_ec->h;

  for (int k = 0; k < 5; k++) {
    if (topo_ec->tetraVol[k] == 0.0) {
      continue;
    }

    int nodes[4];
    double x[4][3];
    double a[3][3];
    double inv[3][3];
    double grad[4][3];
    double gradU[3][3];
    double sigma[3][3];

    for (int i = 0; i < 4; i++) {
      nodes[i] = kTopoTetNodes[part][k][i];
    }

    for (int i = 0; i < 4; i++) {
      int n = nodes[i];
      x[i][0] = kTopoCubeCoord[n][0] * h;
      x[i][1] = kTopoCubeCoord[n][1] * h;
      x[i][2] = kTopoCubeCoord[n][2] * h;
    }

    for (int col = 0; col < 3; col++) {
      for (int row = 0; row < 3; row++) {
        a[row][col] = x[col + 1][row] - x[0][row];
      }
    }

    topo_inv3(a, inv);

    for (int i = 1; i < 4; i++) {
      grad[i][0] = inv[i - 1][0];
      grad[i][1] = inv[i - 1][1];
      grad[i][2] = inv[i - 1][2];
    }
    for (int d = 0; d < 3; d++) {
      grad[0][d] = -(grad[1][d] + grad[2][d] + grad[3][d]);
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        gradU[i][j] = 0.0;
      }
    }

    for (int aidx = 0; aidx < 4; aidx++) {
      int n = nodes[aidx];
      for (int comp = 0; comp < 3; comp++) {
        for (int d = 0; d < 3; d++) {
          gradU[comp][d] += un[n].f[comp] * grad[aidx][d];
        }
      }
    }

    double trace = gradU[0][0] + gradU[1][1] + gradU[2][2];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sigma[i][j] = topo_ec->mu * (gradU[i][j] + gradU[j][i]);
        if (i == j) {
          sigma[i][j] += topo_ec->lambda * trace;
        }
      }
    }

    double volume = fabs(topo_det3(a)) * topo_ec->tetraVol[k] / 6.0;
    for (int aidx = 0; aidx < 4; aidx++) {
      int n = nodes[aidx];
      for (int comp = 0; comp < 3; comp++) {
        double value = 0.0;
        for (int d = 0; d < 3; d++) {
          value += sigma[comp][d] * grad[aidx][d];
        }
        resVec[n].f[comp] += (solver_float)(volume * value);
      }
    }
  }
}

__global__ static void kernelTopoForceEffective(int32_t lenum,
                                                gpu_data_t* gpuData,
                                                int32_t topoElemCount,
                                                const int32_t* topoElemMapper,
                                                const topoconstants_t* topoConstants,
                                                double dt2) {
  int32_t topoIdx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (topoIdx >= topoElemCount) {
    return;
  }

  int32_t eindex = topoElemMapper[topoIdx];
  topoconstants_t topo_ec = topoConstants[topoIdx];
  fvector_t localForce[8];
  fvector_t curDisp[8];
  double c1 = gpuData->c1ArrayDevice[eindex];
  double b_over_dt;
  int do_contrib = 0;

  if (c1 == 0.0) {
    return;
  }

  b_over_dt = gpuData->c3ArrayDevice[eindex] / c1;

  for (int i = 0; i < 8; i++) {
    int32_t lnid = gpuData->lnidArrayDevice[i * lenum + eindex];
    const fvector_t* tm1Disp = gpuData->tm1Device + lnid;
    const fvector_t* tm2Disp = gpuData->tm2Device + lnid;

    localForce[i].f[0] = 0;
    localForce[i].f[1] = 0;
    localForce[i].f[2] = 0;

    curDisp[i].f[0] =
        tm1Disp->f[0] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[0];
    curDisp[i].f[1] =
        tm1Disp->f[1] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[1];
    curDisp[i].f[2] =
        tm1Disp->f[2] * (1.0 + b_over_dt) - b_over_dt * tm2Disp->f[2];

    do_contrib =
        (fabs(curDisp[i].f[0]) > 1e-20) ? 1 : do_contrib;
    do_contrib =
        (fabs(curDisp[i].f[1]) > 1e-20) ? 1 : do_contrib;
    do_contrib =
        (fabs(curDisp[i].f[2]) > 1e-20) ? 1 : do_contrib;
  }

  if (!do_contrib) {
    return;
  }

  topoTetraForcesEffective(curDisp, localForce, &topo_ec);

  for (int i = 0; i < 8; i++) {
    int32_t lnid = gpuData->lnidArrayDevice[i * lenum + eindex];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;

    atomicAdd((solver_float*)&(nodalForce->f[0]),
              (solver_float)(-localForce[i].f[0] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[1]),
              (solver_float)(-localForce[i].f[1] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[2]),
              (solver_float)(-localForce[i].f[2] * dt2));
  }
}

/* ----------------------------------------------------------------------- *
 *  VT + BKT GPU stabilization: faithful CUDA port of the CPU cut-cell
 *  visco-elastic operator TetraForcesBKT (topography.cc). The body below
 *  is copied VERBATIM from the CPU source (the per-tet topoC1/C2/C3 / prs
 *  / resVec algebra, both cube_part variants, all 5 tetrahedra) so the
 *  GPU result matches the validated CPU reference to ported-kernel float
 *  tolerance. un/un_kappa are the same damping_vector_shear/kappa the BKT
 *  conv kernel already produces; mu/kappa are the CPU rho*Vs^2 /
 *  rho*(Vp^2-4/3 Vs^2) moduli (computed host-side, bit-identical).
 * ----------------------------------------------------------------------- */
__device__ static void topoTetraForcesBKT(const fvector_t* un,
                                           const fvector_t* un_kappa,
                                           fvector_t* resVec,
                                           double edgesize,
                                           double mu, double kappa,
                                           const double* tetraVol,
                                           int cube_part) {
  int k;
  int32_t N0, N1, N2, N3;
  double prs;
    double VTetr = edgesize * edgesize * edgesize / 6.0; /* full tetrahedron volume */

    /*  distribution for the first and third quadrants */
    if ( cube_part == 1 )
    {
        for ( k = 0; k < 5; k++ ) { /* for each tetrahedron */

            if ( k == 4 )
                VTetr = 2.0 * VTetr;

            if ( tetraVol[k] != 0 ) {

                double topoC1 = tetraVol[k] * VTetr * kappa / ( edgesize * edgesize );
                double topoC2 = tetraVol[k] * VTetr * mu / ( edgesize * edgesize );
                double topoC3 = tetraVol[k] * VTetr * ( -2. * mu / 3. ) / ( edgesize * edgesize );

                switch ( k ) {

                case ( 0 ):
                N0 = 0;
                N1 = 2;
                N2 = 1;
                N3 = 4;

                prs = topoC3 * (       un[N0].f[0] +       un[N0].f[1] +       un[N0].f[2] -       un[N1].f[1] -       un[N2].f[0] -       un[N3].f[2] ) +
                      topoC1 * ( un_kappa[N0].f[0] + un_kappa[N0].f[1] + un_kappa[N0].f[2] - un_kappa[N1].f[1] - un_kappa[N2].f[0] - un_kappa[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4. * un[N0].f[0] + un[N0].f[1] + un[N0].f[2] - un[N1].f[0] - 2. * un[N2].f[0] - un[N2].f[1] - un[N2].f[2] - un[N3].f[0] );
                resVec[N0].f[1] +=  prs + topoC2 * ( un[N0].f[0] + 4.  * un[N0].f[1] + un[N0].f[2] - un[N1].f[0] - 2.  * un[N1].f[1] - un[N1].f[2] - un[N2].f[1] - un[N3].f[1] );
                resVec[N0].f[2] +=  prs + topoC2 * ( un[N0].f[0] + un[N0].f[1] + 4. * un[N0].f[2] - un[N1].f[2] - un[N2].f[2] - un[N3].f[0] - un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[0] + un[N2].f[1] );
                resVec[N1].f[1] += -prs + topoC2 * ( -2. * ( un[N0].f[1] - un[N1].f[1] ) );
                resVec[N1].f[2] +=        topoC2 * ( -un[N0].f[1] - un[N0].f[2] + un[N1].f[2] + un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( -2. * ( un[N0].f[0] - un[N2].f[0] ) );
                resVec[N2].f[1] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[0] + un[N2].f[1] );
                resVec[N2].f[2] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[2] + un[N2].f[2] + un[N3].f[0]  );
                /* ================ */
                resVec[N3].f[0] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[2] + un[N2].f[2] + un[N3].f[0]  );
                resVec[N3].f[1] +=        topoC2 * ( -un[N0].f[1] - un[N0].f[2] + un[N1].f[2] + un[N3].f[1] );
                resVec[N3].f[2] += -prs + topoC2 * ( -2. * ( un[N0].f[2] - un[N3].f[2] ) );
                break;

                case ( 1 ):
                N0 = 3;
                N1 = 1;
                N2 = 2;
                N3 = 7;

                prs              = topoC3 * (       un[N0].f[0]       + un[N0].f[1]       - un[N0].f[2]       - un[N1].f[1]       - un[N2].f[0]       + un[N3].f[2] ) +
                                   topoC1 * ( un_kappa[N0].f[0] + un_kappa[N0].f[1] - un_kappa[N0].f[2] - un_kappa[N1].f[1] - un_kappa[N2].f[0] + un_kappa[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4. * un[N0].f[0] + un[N0].f[1] - un[N0].f[2] - un[N1].f[0] - 2. * un[N2].f[0] - un[N2].f[1] + un[N2].f[2] - un[N3].f[0] );
                resVec[N0].f[1] +=  prs + topoC2 * ( un[N0].f[0] + 4. * un[N0].f[1] - un[N0].f[2] - un[N1].f[0] - 2. * un[N1].f[1] + un[N1].f[2] - un[N2].f[1] - un[N3].f[1] );
                resVec[N0].f[2] += -prs + topoC2 * ( -un[N0].f[0] - un[N0].f[1] + 4. * un[N0].f[2] - un[N1].f[2] - un[N2].f[2] + un[N3].f[0] + un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[0] + un[N2].f[1] );
                resVec[N1].f[1] += -prs + topoC2 * ( -2. * ( un[N0].f[1] - un[N1].f[1] ) );
                resVec[N1].f[2] +=        topoC2 * ( un[N0].f[1] - un[N0].f[2] + un[N1].f[2] - un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( -2. * ( un[N0].f[0] - un[N2].f[0] ) );
                resVec[N2].f[1] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[0] + un[N2].f[1] );
                resVec[N2].f[2] +=        topoC2 * ( un[N0].f[0] - un[N0].f[2] + un[N2].f[2] - un[N3].f[0]  );
                /* ================ */
                resVec[N3].f[0] +=       -topoC2 * ( un[N0].f[0] - un[N0].f[2] + un[N2].f[2] - un[N3].f[0]  );
                resVec[N3].f[1] +=       -topoC2 * ( un[N0].f[1] - un[N0].f[2] + un[N1].f[2] - un[N3].f[1] );
                resVec[N3].f[2] +=  prs + topoC2 * ( -2. * ( un[N0].f[2] - un[N3].f[2] ) );
                break;

                case ( 2 ):
                N0 = 6;
                N1 = 4;
                N2 = 7;
                N3 = 2;

                prs              = topoC3 * (       un[N0].f[0] -       un[N0].f[1] -       un[N0].f[2] +       un[N1].f[1] -       un[N2].f[0] +       un[N3].f[2] ) +
                                   topoC1 * ( un_kappa[N0].f[0] - un_kappa[N0].f[1] - un_kappa[N0].f[2] + un_kappa[N1].f[1] - un_kappa[N2].f[0] + un_kappa[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4. * un[N0].f[0] - un[N0].f[1] - un[N0].f[2] - un[N1].f[0] - 2. * un[N2].f[0] + un[N2].f[1] + un[N2].f[2] - un[N3].f[0] );

                resVec[N0].f[1] += -prs + topoC2 * ( -un[N0].f[0] + 4. * un[N0].f[1] + un[N0].f[2] + un[N1].f[0] - 2. * un[N1].f[1] - un[N1].f[2] - un[N2].f[1] - un[N3].f[1] );
                resVec[N0].f[2] += -prs + topoC2 * ( -un[N0].f[0] + un[N0].f[1] + 4. * un[N0].f[2] - un[N1].f[2] - un[N2].f[2] + un[N3].f[0] - un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * ( -un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N2].f[1] );
                resVec[N1].f[1] +=  prs + topoC2 * ( -2. * ( un[N0].f[1] - un[N1].f[1] ) );
                resVec[N1].f[2] +=        topoC2 * ( -un[N0].f[1] - un[N0].f[2] + un[N1].f[2] + un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( -2. * ( un[N0].f[0] - un[N2].f[0] ) );
                resVec[N2].f[1] +=       -topoC2 * ( -un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N2].f[1] );
                resVec[N2].f[2] +=        topoC2 * (  un[N0].f[0] - un[N0].f[2] + un[N2].f[2] - un[N3].f[0]  );
                /* ================ */
                resVec[N3].f[0] +=       -topoC2 * ( un[N0].f[0] - un[N0].f[2] + un[N2].f[2] - un[N3].f[0]  );
                resVec[N3].f[1] +=        topoC2 * ( -un[N0].f[1] - un[N0].f[2] + un[N1].f[2] + un[N3].f[1] );
                resVec[N3].f[2] +=  prs + topoC2 * ( -2. * ( un[N0].f[2] - un[N3].f[2] ) );
                break;

                case ( 3 ):
                N0 = 5;
                N1 = 7;
                N2 = 4;
                N3 = 1;

                prs              = topoC3 * (       un[N0].f[0] -       un[N0].f[1] +       un[N0].f[2] +       un[N1].f[1] -       un[N2].f[0] -       un[N3].f[2] ) +
                                   topoC1 * ( un_kappa[N0].f[0] - un_kappa[N0].f[1] + un_kappa[N0].f[2] + un_kappa[N1].f[1] - un_kappa[N2].f[0] - un_kappa[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4. * un[N0].f[0] - un[N0].f[1] + un[N0].f[2] - un[N1].f[0] - 2. * un[N2].f[0] + un[N2].f[1] - un[N2].f[2] - un[N3].f[0] );
                resVec[N0].f[1] += -prs + topoC2 * ( -un[N0].f[0] + 4. * un[N0].f[1] - un[N0].f[2] + un[N1].f[0] - 2. * un[N1].f[1] + un[N1].f[2] - un[N2].f[1] - un[N3].f[1] );
                resVec[N0].f[2] +=  prs + topoC2 * (  un[N0].f[0] - un[N0].f[1] + 4. * un[N0].f[2] - un[N1].f[2] - un[N2].f[2] - un[N3].f[0] + un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * ( -un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N2].f[1] );
                resVec[N1].f[1] +=  prs + topoC2 * ( -2. * ( un[N0].f[1] - un[N1].f[1] ) );
                resVec[N1].f[2] +=        topoC2 * ( un[N0].f[1] - un[N0].f[2] + un[N1].f[2] - un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs - topoC2 * ( un[N0].f[0] - un[N2].f[0] ) * 2.0;
                resVec[N2].f[1] +=       -topoC2 * ( -un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N2].f[1] );
                resVec[N2].f[2] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[2] + un[N2].f[2] + un[N3].f[0] );
                /* ================ */
                resVec[N3].f[0] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[2] + un[N2].f[2] + un[N3].f[0] );
                resVec[N3].f[1] +=       -topoC2 * ( un[N0].f[1] - un[N0].f[2] + un[N1].f[2] - un[N3].f[1] );
                resVec[N3].f[2] += -prs + topoC2 * ( -2. * ( un[N0].f[2] - un[N3].f[2] ) );
                break;

                case ( 4 ):
                N0 = 2;
                N1 = 4;
                N2 = 7;
                N3 = 1;

                topoC1 = topoC1 / 4.;
                topoC2 = topoC2 / 4.;
                topoC3 = topoC3 / 4.;

                prs              = topoC1 * ( un_kappa[N0].f[0] - un_kappa[N0].f[1] + un_kappa[N0].f[2] + un_kappa[N1].f[0] + un_kappa[N1].f[1] - un_kappa[N1].f[2] - un_kappa[N2].f[0] - un_kappa[N2].f[1] - un_kappa[N2].f[2] - un_kappa[N3].f[0] + un_kappa[N3].f[1] + un_kappa[N3].f[2] ) +
                                   topoC3 * (       un[N0].f[0] -       un[N0].f[1] +       un[N0].f[2] +       un[N1].f[0] +       un[N1].f[1] -       un[N1].f[2] -       un[N2].f[0] -       un[N2].f[1] -       un[N2].f[2] -       un[N3].f[0] +       un[N3].f[1] +       un[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4. * un[N0].f[0] - un[N0].f[1] + un[N0].f[2] - un[N1].f[1] + un[N1].f[2] - 2. * un[N2].f[0] + un[N2].f[1] - un[N2].f[2] - 2. * un[N3].f[0] + un[N3].f[1] - un[N3].f[2] );
                resVec[N0].f[1] += -prs + topoC2 * ( -un[N0].f[0] + 4. * un[N0].f[1] - un[N0].f[2] + un[N1].f[0] - 2. * un[N1].f[1] + un[N1].f[2] - un[N2].f[0] - un[N2].f[2] + un[N3].f[0] - 2. * un[N3].f[1] + un[N3].f[2] );
                resVec[N0].f[2] +=  prs + topoC2 * ( un[N0].f[0] - un[N0].f[1] + 4. * un[N0].f[2] - un[N1].f[0] + un[N1].f[1] - 2. * un[N1].f[2] - un[N2].f[0] + un[N2].f[1] - 2. * un[N2].f[2] + un[N3].f[0] - un[N3].f[1] );
                /* ================ */
                resVec[N1].f[0] +=  prs + topoC2 * (  un[N0].f[1] - un[N0].f[2] + 4. * un[N1].f[0] + un[N1].f[1] - un[N1].f[2] - 2. * un[N2].f[0] - un[N2].f[1] + un[N2].f[2] - 2. * un[N3].f[0] - un[N3].f[1] + un[N3].f[2]  );
                resVec[N1].f[1] +=  prs + topoC2 * ( -un[N0].f[0] - 2. * un[N0].f[1] + un[N0].f[2] + un[N1].f[0] + 4. * un[N1].f[1] - un[N1].f[2] - un[N2].f[0] - 2. * un[N2].f[1] + un[N2].f[2] + un[N3].f[0] - un[N3].f[2]  );
                resVec[N1].f[2] += -prs + topoC2 * (  un[N0].f[0] + un[N0].f[1] - 2. * un[N0].f[2] - un[N1].f[0] - un[N1].f[1] + 4. * un[N1].f[2] - un[N2].f[0]  - un[N2].f[1] + un[N3].f[0] + un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 *  ( - 2. * un[N0].f[0] - un[N0].f[1] - un[N0].f[2] - 2. * un[N1].f[0] - un[N1].f[1] - un[N1].f[2] + 4. * un[N2].f[0] + un[N2].f[1] + un[N2].f[2] + un[N3].f[1] + un[N3].f[2] );
                resVec[N2].f[1] += -prs + topoC2 *  (  un[N0].f[0] + un[N0].f[2] - un[N1].f[0] - 2. * un[N1].f[1] - un[N1].f[2] + un[N2].f[0] + 4. * un[N2].f[1] + un[N2].f[2] - un[N3].f[0] - 2. * un[N3].f[1] - un[N3].f[2] );
                resVec[N2].f[2] += -prs + topoC2 *  ( -un[N0].f[0] - un[N0].f[1] - 2. * un[N0].f[2] + un[N1].f[0] + un[N1].f[1] + un[N2].f[0] + un[N2].f[1] + 4. * un[N2].f[2] - un[N3].f[0] - un[N3].f[1] - 2. * un[N3].f[2] );
                /* ================ */
                resVec[N3].f[0] += -prs + topoC2 *  ( - 2. * un[N0].f[0] + un[N0].f[1] + un[N0].f[2] - 2. * un[N1].f[0] + un[N1].f[1] + un[N1].f[2] - un[N2].f[1] - un[N2].f[2] + 4. * un[N3].f[0] - un[N3].f[1] - un[N3].f[2] );
                resVec[N3].f[1] +=  prs + topoC2 *  (  un[N0].f[0] - 2. * un[N0].f[1] - un[N0].f[2] - un[N1].f[0] + un[N1].f[2] + un[N2].f[0] - 2. * un[N2].f[1] - un[N2].f[2] - un[N3].f[0] + 4. * un[N3].f[1] + un[N3].f[2] );
                resVec[N3].f[2] +=  prs + topoC2 *  ( -un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N1].f[1] - 2.0 * un[N1].f[2] + un[N2].f[0] - un[N2].f[1] - 2.0 * un[N2].f[2] - un[N3].f[0] + un[N3].f[1] + 4.0 * un[N3].f[2] );
                break;
                }
            }
        }
    }  else  {
        /*  distribution for the second and fourth quadrants */
        for ( k = 0; k < 5; k++ ) { /* for each tetrahedron */

            if ( k == 4 )
                VTetr = 2.0 * VTetr;

            if ( tetraVol[k] != 0 ) {

                double topoC1 = tetraVol[k] * VTetr * kappa / ( edgesize * edgesize );
                double topoC2 = tetraVol[k] * VTetr * mu / ( edgesize * edgesize );
                double topoC3 = tetraVol[k] * VTetr * ( -2. * mu / 3. ) / ( edgesize * edgesize );

                switch ( k ) {

                case ( 0 ):
                N0 = 0;
                N1 = 3;
                N2 = 1;
                N3 = 5;

                prs              = topoC1 * ( un_kappa[N0].f[0] - un_kappa[N1].f[1] - un_kappa[N3].f[2] - un_kappa[N2].f[0] + un_kappa[N2].f[1] + un_kappa[N2].f[2] ) +
                                   topoC3 * (       un[N0].f[0] -       un[N1].f[1] -       un[N3].f[2] -       un[N2].f[0] +       un[N2].f[1] +       un[N2].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * 2.0 * ( un[N0].f[0] - un[N2].f[0] );
                resVec[N0].f[1] +=        topoC2 * ( un[N0].f[1] - un[N1].f[0] + un[N2].f[0] - un[N2].f[1] );
                resVec[N0].f[2] +=        topoC2 * ( un[N0].f[2]  + un[N2].f[0] - un[N2].f[2] - un[N3].f[0] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * ( -un[N0].f[1] + un[N1].f[0] - un[N2].f[0] + un[N2].f[1] );
                resVec[N1].f[1] += -prs + topoC2 * 2.0 * ( un[N1].f[1] - un[N2].f[1] );
                resVec[N1].f[2] +=        topoC2 * ( un[N1].f[2] - un[N2].f[1] - un[N2].f[2] + un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( -2.0 * un[N0].f[0] + un[N0].f[1] + un[N0].f[2] - un[N1].f[0] + 4.0 * un[N2].f[0] - un[N2].f[1] - un[N2].f[2] - un[N3].f[0] );
                resVec[N2].f[1] +=  prs + topoC2 * ( -un[N0].f[1] + un[N1].f[0] - 2.0 * un[N1].f[1] - un[N1].f[2] - un[N2].f[0] + 4.0 * un[N2].f[1] + un[N2].f[2] - un[N3].f[1] );
                resVec[N2].f[2] +=  prs + topoC2 * ( -un[N0].f[2] - un[N1].f[2] - un[N2].f[0] + un[N2].f[1] + 4.0 * un[N2].f[2] + un[N3].f[0] - un[N3].f[1] - 2.0 * un[N3].f[2] );
                /* ================ */
                resVec[N3].f[0] +=        topoC2 * ( -un[N0].f[2] - un[N2].f[0] + un[N2].f[2] + un[N3].f[0] );
                resVec[N3].f[1] +=        topoC2 * (  un[N1].f[2] - un[N2].f[1] - un[N2].f[2] + un[N3].f[1] );
                resVec[N3].f[2] += -prs + topoC2 * 2.0 * ( un[N3].f[2] - un[N2].f[2] );

                /* ================ */
                break;

                case ( 1 ):
                N0 = 0;
                N1 = 2;
                N2 = 3;
                N3 = 6;

                prs              = topoC1 * ( un_kappa[N0].f[1] - un_kappa[N2].f[0] - un_kappa[N3].f[2] + un_kappa[N1].f[0] - un_kappa[N1].f[1] + un_kappa[N1].f[2] ) +
                                   topoC3 * (       un[N0].f[1] -       un[N2].f[0] -       un[N3].f[2] +       un[N1].f[0] -       un[N1].f[1] +       un[N1].f[2] );

                resVec[N0].f[0] +=        topoC2 * ( un[N0].f[0] - un[N1].f[0] + un[N1].f[1] - un[N2].f[1] );
                resVec[N0].f[1] +=  prs + topoC2 * 2.0 * ( un[N0].f[1] - un[N1].f[1] );
                resVec[N0].f[2] +=        topoC2 * ( un[N0].f[2] + un[N1].f[1] - un[N1].f[2] - un[N3].f[1] );
                /* ================ */
                resVec[N1].f[0] +=  prs + topoC2 * ( -un[N0].f[0] + 4.0 * un[N1].f[0] - un[N1].f[1] + un[N1].f[2] - 2.0 * un[N2].f[0] + un[N2].f[1] - un[N2].f[2] - un[N3].f[0] );
                resVec[N1].f[1] += -prs + topoC2 * ( un[N0].f[0] - 2.0 * un[N0].f[1] + un[N0].f[2] - un[N1].f[0] + 4.0 * un[N1].f[1] - un[N1].f[2] - un[N2].f[1] - un[N3].f[1] );
                resVec[N1].f[2] +=  prs + topoC2 * ( -un[N0].f[2] + un[N1].f[0] - un[N1].f[1] + 4.0 * un[N1].f[2] - un[N2].f[2] - un[N3].f[0] + un[N3].f[1] - 2.0 * un[N3].f[2] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * 2.0 * ( un[N2].f[0] - un[N1].f[0] );
                resVec[N2].f[1] +=        topoC2 * ( -un[N0].f[0] + un[N1].f[0] - un[N1].f[1] + un[N2].f[1] );
                resVec[N2].f[2] +=        topoC2 * ( -un[N1].f[0] - un[N1].f[2] + un[N2].f[2] + un[N3].f[0] );
                /* ================ */
                resVec[N3].f[0] +=        topoC2 * ( -un[N1].f[0] - un[N1].f[2] + un[N2].f[2] + un[N3].f[0] );
                resVec[N3].f[1] +=        topoC2 * ( -un[N0].f[2] - un[N1].f[1] + un[N1].f[2] + un[N3].f[1] );
                resVec[N3].f[2] += -prs + topoC2 * 2.0 * ( un[N3].f[2] - un[N1].f[2] );
                break;

                case ( 2 ):
                N0 = 4;
                N1 = 5;
                N2 = 6;
                N3 = 0;

                prs              = topoC1 * ( un_kappa[N3].f[2] - un_kappa[N2].f[1] - un_kappa[N1].f[0] + un_kappa[N0].f[0] + un_kappa[N0].f[1] - un_kappa[N0].f[2] ) +
                                   topoC3 * (       un[N3].f[2] -       un[N2].f[1] -       un[N1].f[0] +       un[N0].f[0] +       un[N0].f[1] -       un[N0].f[2] );


                resVec[N0].f[0] +=  prs + topoC2 * ( 4.0 * un[N0].f[0] + un[N0].f[1] - un[N0].f[2] - 2.0 * un[N1].f[0] - un[N1].f[1] + un[N1].f[2] - un[N2].f[0] - un[N3].f[0] );
                resVec[N0].f[1] +=  prs + topoC2 * ( 4.0 * un[N0].f[1] + un[N0].f[0] - un[N0].f[2] - un[N1].f[1] - 2.0 * un[N2].f[1] - un[N2].f[0] + un[N2].f[2] - un[N3].f[1] );
                resVec[N0].f[2] += -prs + topoC2 * ( 4.0 * un[N0].f[2] - un[N0].f[0] - un[N0].f[1] - un[N1].f[2] - un[N2].f[2] - 2.0 * un[N3].f[2] + un[N3].f[0] + un[N3].f[1] );
                /* ================ */
                resVec[N1].f[0] += -prs + topoC2 * 2.0 * ( un[N1].f[0] - un[N0].f[0] );
                resVec[N1].f[1] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[1] + un[N2].f[0] );
                resVec[N1].f[2] +=        topoC2 * (  un[N0].f[0] - un[N0].f[2] + un[N1].f[2] - un[N3].f[0] );
                /* ================ */
                resVec[N2].f[0] +=        topoC2 * ( -un[N0].f[0] - un[N0].f[1] + un[N1].f[1] + un[N2].f[0] );
                resVec[N2].f[1] += -prs + topoC2 * 2.0 * ( un[N2].f[1] - un[N0].f[1] );
                resVec[N2].f[2] +=        topoC2 * (  un[N0].f[1] - un[N0].f[2] + un[N2].f[2] - un[N3].f[1] );
                /* ================ */
                resVec[N3].f[0] +=        topoC2 * ( -un[N0].f[0] + un[N0].f[2] - un[N1].f[2] + un[N3].f[0] );
                resVec[N3].f[1] +=        topoC2 * ( -un[N0].f[1] + un[N0].f[2] - un[N2].f[2] + un[N3].f[1] );
                resVec[N3].f[2] +=  prs + topoC2 * 2.0 * ( un[N3].f[2] - un[N0].f[2] );
                break;

                case ( 3 ):
                N0 = 6;
                N1 = 5;
                N2 = 7;
                N3 = 3;

                prs              = topoC1 * ( un_kappa[N0].f[0] + un_kappa[N1].f[1] + un_kappa[N3].f[2] - un_kappa[N2].f[0] - un_kappa[N2].f[1] - un_kappa[N2].f[2] ) +
                                   topoC3 * (       un[N0].f[0] +       un[N1].f[1] +       un[N3].f[2] -       un[N2].f[0] -       un[N2].f[1] -       un[N2].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * 2.0 * ( un[N0].f[0] - un[N2].f[0] );
                resVec[N0].f[1] +=        topoC2 * (  un[N0].f[1] + un[N1].f[0] - un[N2].f[0] - un[N2].f[1] );
                resVec[N0].f[2] +=        topoC2 * (  un[N0].f[2] - un[N2].f[0] - un[N2].f[2] + un[N3].f[0] );
                /* ================ */
                resVec[N1].f[0] +=        topoC2 * (  un[N0].f[1] + un[N1].f[0] - un[N2].f[0] - un[N2].f[1] );
                resVec[N1].f[1] +=  prs + topoC2 * 2.0 * ( un[N1].f[1] - un[N2].f[1] );
                resVec[N1].f[2] +=        topoC2 * (  un[N1].f[2] - un[N2].f[1] - un[N2].f[2] + un[N3].f[1] );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( -2.0 * un[N0].f[0] - un[N0].f[1] - un[N0].f[2] - un[N1].f[0] + 4.0 * un[N2].f[0] + un[N2].f[1] + un[N2].f[2] - un[N3].f[0] );
                resVec[N2].f[1] += -prs + topoC2 * ( -un[N0].f[1] - un[N1].f[0] - 2.0 * un[N1].f[1] - un[N1].f[2] + un[N2].f[0] + 4.0 * un[N2].f[1] + un[N2].f[2] - un[N3].f[1] );
                resVec[N2].f[2] += -prs + topoC2 * ( -un[N0].f[2] - un[N1].f[2] + un[N2].f[0] + un[N2].f[1] + 4.0 * un[N2].f[2] - un[N3].f[0] - un[N3].f[1] - 2.0 * un[N3].f[2] );

                /* ================ */
                resVec[N3].f[0] +=        topoC2 * (  un[N0].f[2] - un[N2].f[0] - un[N2].f[2] + un[N3].f[0] );
                resVec[N3].f[1] +=        topoC2 * (  un[N1].f[2] - un[N2].f[1] - un[N2].f[2] + un[N3].f[1] );
                resVec[N3].f[2] +=  prs + topoC2 * 2.0 * ( un[N3].f[2] - un[N2].f[2] );
                break;

                case ( 4 ):
                N0 = 0;
                N1 = 6;
                N2 = 3;
                N3 = 5;

                topoC1 = topoC1 / 4.0;
                topoC2 = topoC2 / 4.0;
                topoC3 = topoC3 / 4.0;

                prs              = topoC1 * ( un_kappa[N0].f[0] + un_kappa[N0].f[1] + un_kappa[N0].f[2] + un_kappa[N1].f[0] - un_kappa[N1].f[1] - un_kappa[N1].f[2] - un_kappa[N2].f[0] - un_kappa[N2].f[1] + un_kappa[N2].f[2] - un_kappa[N3].f[0] + un_kappa[N3].f[1] - un_kappa[N3].f[2] ) +
                                   topoC3 * (       un[N0].f[0] +       un[N0].f[1] +       un[N0].f[2] +       un[N1].f[0] -       un[N1].f[1] -       un[N1].f[2] -       un[N2].f[0] -       un[N2].f[1] +       un[N2].f[2] -       un[N3].f[0] +       un[N3].f[1] -       un[N3].f[2] );

                resVec[N0].f[0] +=  prs + topoC2 * ( 4.0 * un[N0].f[0] + un[N0].f[1] + un[N0].f[2] + un[N1].f[1] + un[N1].f[2] - un[N2].f[1] - un[N2].f[2] - un[N3].f[1] - un[N3].f[2] - 2.0 * ( un[N2].f[0] + un[N3].f[0] ) );
                resVec[N0].f[1] +=  prs + topoC2 * (  un[N0].f[0] + 4.0 * un[N0].f[1] + un[N0].f[2] - un[N1].f[0] - un[N1].f[2] - un[N2].f[0] - un[N2].f[2] + un[N3].f[0] + un[N3].f[2] - 2.0 * ( un[N1].f[1] + un[N2].f[1] ) );
                resVec[N0].f[2] +=  prs + topoC2 * (  un[N0].f[0] + un[N0].f[1] + 4.0 * un[N0].f[2] - un[N1].f[0] - un[N1].f[1] + un[N2].f[0] + un[N2].f[1] - un[N3].f[0] - un[N3].f[1] - 2.0 * ( un[N1].f[2] + un[N3].f[2] ) );
                /* ================ */
                resVec[N1].f[0] +=  prs + topoC2 * ( -un[N0].f[1] - un[N0].f[2] + 4.0 * un[N1].f[0] -       un[N1].f[1] -       un[N1].f[2] + un[N2].f[1] + un[N2].f[2] + un[N3].f[1] + un[N3].f[2] - 2.0 * ( un[N3].f[0] + un[N2].f[0] ) );
                resVec[N1].f[1] += -prs + topoC2 * (  un[N0].f[0] - un[N0].f[2] -       un[N1].f[0] + 4.0 * un[N1].f[1] +       un[N1].f[2] - un[N2].f[0] + un[N2].f[2] + un[N3].f[0] - un[N3].f[2] - 2.0 * ( un[N3].f[1] + un[N0].f[1] ) );
                resVec[N1].f[2] += -prs + topoC2 * (  un[N0].f[0] - un[N0].f[1] -       un[N1].f[0] +       un[N1].f[1] + 4.0 * un[N1].f[2] +  un[N2].f[0] - un[N2].f[1] - un[N3].f[0] + un[N3].f[1] - 2.0 * ( un[N0].f[2] + un[N2].f[2] ) );
                /* ================ */
                resVec[N2].f[0] += -prs + topoC2 * ( 4.0 * un[N2].f[0] - un[N0].f[1] + un[N0].f[2] - un[N1].f[1] + un[N1].f[2] + un[N2].f[1] - un[N2].f[2] + un[N3].f[1] - un[N3].f[2] - 2.0 * ( un[N0].f[0] + un[N1].f[0] ) );
                resVec[N2].f[1] += -prs + topoC2 * ( 4.0 * un[N2].f[1] - un[N0].f[0] + un[N0].f[2] + un[N1].f[0] - un[N1].f[2] + un[N2].f[0] - un[N2].f[2] - un[N3].f[0] + un[N3].f[2] - 2.0 * ( un[N0].f[1] + un[N3].f[1] ) );
                resVec[N2].f[2] +=  prs + topoC2 * ( 4.0 * un[N2].f[2] - un[N0].f[0] - un[N0].f[1] + un[N1].f[0] + un[N1].f[1] - un[N2].f[0] - un[N2].f[1] + un[N3].f[0] + un[N3].f[1] - 2.0 * ( un[N3].f[2] + un[N1].f[2] ) );
                /* ================ */
                resVec[N3].f[0] += -prs + topoC2 * ( 4.0 * un[N3].f[0] + un[N0].f[1] - un[N0].f[2] + un[N1].f[1] - un[N1].f[2] - un[N2].f[1] + un[N2].f[2] - un[N3].f[1] + un[N3].f[2] - 2.0 * ( un[N0].f[0] + un[N1].f[0] ) );
                resVec[N3].f[1] +=  prs + topoC2 * ( 4.0 * un[N3].f[1] - un[N0].f[0] - un[N0].f[2] + un[N1].f[0] + un[N1].f[2] + un[N2].f[0] + un[N2].f[2] - un[N3].f[0] - un[N3].f[2] - 2.0 * ( un[N1].f[1] + un[N2].f[1] ) );
                resVec[N3].f[2] += -prs + topoC2 * ( 4.0 * un[N3].f[2] - un[N0].f[0] + un[N0].f[1] + un[N1].f[0] - un[N1].f[1] - un[N2].f[0] + un[N2].f[1] + un[N3].f[0] - un[N3].f[1] - 2.0 * ( un[N0].f[2] + un[N2].f[2] ) );
                break;

                }
            }
        }
    }
}

/* GPU port of the CPU topoBkt branch (damping.cc): for each VT cut cell
 * flagged edata->topoBkt, apply the fill-scaled 5-tet visco-elastic force
 * from the BKT damping vectors. shear/kappa vectors are read with the SAME
 * indexing kernelDampingCalcLocal uses (numcomp = lenum*8), so they are the
 * bit-identical damping_vector_shear/kappa the conv kernel already wrote.
 * Force sign mirrors the CPU: nodalForce -= localForce * dt^2.
 * CPML: defensive PML skip mirrors kernelDampingCalcConvBkt2/Local -- a
 * PML-flagged cell had its shear/kappa vectors skipped by the conv kernel,
 * so it must not get a force here (its attenuation is the elastic CPML
 * path's job). No-op under dashpot (pmlFlagArrayDevice == NULL). The builder
 * additionally excludes cpml_pml_flattened_element cells, so this is
 * belt-and-suspenders. */
__global__ static void kernelTopoBktForceEffective(
    int32_t lenum, gpu_data_t* gpuData, int32_t topoBktCount,
    const int32_t* topoBktMapper, const topo_bkt_const_t* topoBktConst,
    double dt2) {
  int32_t idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (idx >= topoBktCount) {
    return;
  }

  int32_t eindex = topoBktMapper[idx];

  if (gpuData->pmlFlagArrayDevice != NULL &&
      gpuData->pmlFlagArrayDevice[eindex] != 0) {
    return;
  }

  topo_bkt_const_t tc = topoBktConst[idx];
  int numcomp = lenum * 8;
  fvector_t un[8];
  fvector_t un_kappa[8];
  fvector_t localForce[8];
  int do_contrib = 0;

  solver_float* sVec = (solver_float*)gpuData->shearVectorDevice + eindex;
  solver_float* kVec = (solver_float*)gpuData->kappaVectorDevice + eindex;
  for (int i = 0; i < 8; i++) {
    un[i].f[0] = sVec[0];
    un[i].f[1] = sVec[numcomp];
    un[i].f[2] = sVec[2 * numcomp];
    un_kappa[i].f[0] = kVec[0];
    un_kappa[i].f[1] = kVec[numcomp];
    un_kappa[i].f[2] = kVec[2 * numcomp];
    localForce[i].f[0] = 0;
    localForce[i].f[1] = 0;
    localForce[i].f[2] = 0;
    do_contrib = (fabs(un[i].f[0]) > 1e-20) ? 1 : do_contrib;
    do_contrib = (fabs(un[i].f[1]) > 1e-20) ? 1 : do_contrib;
    do_contrib = (fabs(un[i].f[2]) > 1e-20) ? 1 : do_contrib;
    do_contrib = (fabs(un_kappa[i].f[0]) > 1e-20) ? 1 : do_contrib;
    do_contrib = (fabs(un_kappa[i].f[1]) > 1e-20) ? 1 : do_contrib;
    do_contrib = (fabs(un_kappa[i].f[2]) > 1e-20) ? 1 : do_contrib;
    sVec += lenum;
    kVec += lenum;
  }

  if (!do_contrib) {
    return;
  }

  topoTetraForcesBKT(un, un_kappa, localForce, tc.edgesize, tc.mu_bkt,
                     tc.kappa_bkt, tc.tetraVol, tc.cube_part);

  for (int i = 0; i < 8; i++) {
    int32_t lnid = gpuData->lnidArrayDevice[i * lenum + eindex];
    fvector_t* nodalForce = gpuData->forceDevice + lnid;
    atomicAdd((solver_float*)&(nodalForce->f[0]),
              (solver_float)(-localForce[i].f[0] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[1]),
              (solver_float)(-localForce[i].f[1] * dt2));
    atomicAdd((solver_float*)&(nodalForce->f[2]),
              (solver_float)(-localForce[i].f[2] * dt2));
  }
}

static char* kernel_handle(psolve_gpu_kernel_id_t kernel) {
  switch (kernel) {
    case PSOLVE_GPU_KERNEL_STIFFNESS:
      return (char*)kernelStiffnessCalcLocal;
    case PSOLVE_GPU_KERNEL_DAMPING_BKT2:
      return (char*)kernelDampingCalcConvBkt2;
    case PSOLVE_GPU_KERNEL_DAMPING_BKT2_FUSED:
      return (char*)kernelDampingBkt2Fused;
    case PSOLVE_GPU_KERNEL_DAMPING_BKT3:
      return (char*)kernelDampingCalcConvBkt3;
    case PSOLVE_GPU_KERNEL_DAMPING_GENERIC:
      return (char*)kernelDampingCalcConv;
    case PSOLVE_GPU_KERNEL_DAMPING_FORCE:
      return (char*)kernelDampingCalcLocal;
    case PSOLVE_GPU_KERNEL_DISPLACEMENT:
      return (char*)kernelDispCalc;
    case PSOLVE_GPU_KERNEL_STATION_TENSOR:
      return (char*)kernelStationVEStressGather;
  }
  return NULL;
}

char* psolve_gpu_kernel_handle(psolve_gpu_kernel_id_t kernel) {
  return kernel_handle(kernel);
}

int32_t psolve_gpu_kernel_blocksize(gpu_spec_t* gpuSpecs,
                                   psolve_gpu_kernel_id_t kernel,
                                   int32_t memPerThread) {
  char* handle = kernel_handle(kernel);
  return (handle == NULL) ? 0
                          : gpu_get_blocksize(gpuSpecs, handle, memPerThread);
}

void psolve_launch_station_tensor(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, int32_t outputStride, int32_t strainOffset,
    int32_t stressOffset, gpu_data_t* gpuData) {
  kernelStationVEStressGather<<<grid, block, shared, stream>>>(
      stationCount, outputStride, strainOffset, stressOffset, gpuData);
}

void psolve_launch_source_force_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const vector3D_t* srcForces,
    fvector_t* forceDevice, double dt2) {
  kernelSourceForceScatter<<<grid, block, shared, stream>>>(
      count, lnidList, srcForces, forceDevice, dt2);
}

void psolve_launch_sparse_force_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* srcForces,
    fvector_t* forceDevice) {
  kernelSparseForceScatter<<<grid, block, shared, stream>>>(
      count, lnidList, srcForces, forceDevice);
}

void psolve_launch_disp_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* dispDevice,
    fvector_t* bufferDevice) {
  kernelDispBoundaryGather<<<grid, block, shared, stream>>>(
      count, lnidList, dispDevice, bufferDevice);
}

void psolve_launch_disp_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dispDevice) {
  kernelDispBoundaryScatter<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dispDevice);
}

void psolve_launch_scalar_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* srcDevice,
    solver_float* bufferDevice) {
  kernelScalarBoundaryGather<<<grid, block, shared, stream>>>(
      count, lnidList, srcDevice, bufferDevice);
}

void psolve_launch_scalar_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice) {
  kernelScalarBoundaryScatter<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_scalar_boundary_scatter_serial(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice) {
  kernelScalarBoundaryScatterSerial<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_scalar_boundary_accumulate(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const solver_float* bufferDevice,
    solver_float* dstDevice) {
  kernelScalarBoundaryAccumulate<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_scalar_adjust_distribution(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t runCount,
    const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values) {
  kernelScalarAdjustDistribution<<<grid, block, shared, stream>>>(
      runCount, runOffsets, runParentBase3, edgeDn, danglingBase3, invDeps,
      values);
}

void psolve_launch_scalar_adjust_assignment(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t dnCount,
    const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, solver_float* values) {
  kernelScalarAdjustAssignment<<<grid, block, shared, stream>>>(
      dnCount, dnOrder, parentOffsets, parentBase3, danglingBase3, invDeps,
      values);
}

void psolve_launch_vector_boundary_gather(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* srcDevice,
    fvector_t* bufferDevice) {
  kernelVectorBoundaryGather<<<grid, block, shared, stream>>>(
      count, lnidList, srcDevice, bufferDevice);
}

void psolve_launch_vector_boundary_scatter(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice) {
  kernelVectorBoundaryScatter<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_vector_boundary_scatter_serial(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice) {
  kernelVectorBoundaryScatterSerial<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_vector_boundary_accumulate(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t count,
    const int32_t* lnidList, const fvector_t* bufferDevice,
    fvector_t* dstDevice) {
  kernelVectorBoundaryAccumulate<<<grid, block, shared, stream>>>(
      count, lnidList, bufferDevice, dstDevice);
}

void psolve_launch_vector_adjust_distribution(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t runCount,
    const int32_t* runOffsets, const int32_t* runParentBase3,
    const int32_t* edgeDn, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values) {
  kernelVectorAdjustDistribution<<<grid, block, shared, stream>>>(
      runCount, runOffsets, runParentBase3, edgeDn, danglingBase3, invDeps,
      values);
}

void psolve_launch_vector_adjust_assignment(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t dnCount,
    const int32_t* dnOrder, const int32_t* parentOffsets,
    const int32_t* parentBase3, const int32_t* danglingBase3,
    const solver_float* invDeps, fvector_t* values) {
  kernelVectorAdjustAssignment<<<grid, block, shared, stream>>>(
      dnCount, dnOrder, parentOffsets, parentBase3, danglingBase3, invDeps,
      values);
}

void psolve_launch_stiffness(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t linearElementCount,
    int32_t* linearElementMapper) {
  kernelStiffnessCalcLocal<<<grid, block, shared, stream>>>(
      lenum, gpuData, linearElementCount, linearElementMapper);
}

void psolve_launch_pml_elastic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t pmlElasticCount) {
  kernelPMLElasticStep<<<grid, block, shared, stream>>>(
      lenum, gpuData, pmlElasticCount);
}

void psolve_launch_damping_bkt2(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax,
    double shearCoef1, double shearCoef2, double shearExpCoef0,
    double shearCoef3, double shearCoef4, double shearExpCoef1,
    double kappaCoef1, double kappaCoef2, double kappaExpCoef0,
    double kappaCoef3, double kappaCoef4, double kappaExpCoef1,
    int32_t elementCount, int32_t* elementMapper) {
  kernelDampingCalcConvBkt2<<<grid, block, shared, stream>>>(
      lenum, numcomp, gpuData, rmax, shearCoef1, shearCoef2, shearExpCoef0,
      shearCoef3, shearCoef4, shearExpCoef1, kappaCoef1, kappaCoef2,
      kappaExpCoef0, kappaCoef3, kappaCoef4, kappaExpCoef1, elementCount,
      elementMapper);
}

void psolve_launch_damping_bkt3(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax,
    double shearCoef1, double shearCoef2, double shearExpCoef0,
    double shearCoef3, double shearCoef4, double shearExpCoef1,
    double shearCoef5, double shearCoef6, double shearExpCoef2,
    double kappaCoef1, double kappaCoef2, double kappaExpCoef0,
    double kappaCoef3, double kappaCoef4, double kappaExpCoef1,
    double kappaCoef5, double kappaCoef6, double kappaExpCoef2,
    int32_t elementCount, int32_t* elementMapper) {
  kernelDampingCalcConvBkt3<<<grid, block, shared, stream>>>(
      lenum, numcomp, gpuData, rmax, shearCoef1, shearCoef2, shearExpCoef0,
      shearCoef3, shearCoef4, shearExpCoef1, shearCoef5, shearCoef6,
      shearExpCoef2, kappaCoef1, kappaCoef2, kappaExpCoef0, kappaCoef3,
      kappaCoef4, kappaExpCoef1, kappaCoef5, kappaCoef6, kappaExpCoef2,
      elementCount, elementMapper);
}

void psolve_launch_damping_generic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    int32_t numcomp, gpu_data_t* gpuData, double rmax, int32_t dampingType,
    int32_t elementCount, int32_t* elementMapper) {
  kernelDampingCalcConv<<<grid, block, shared, stream>>>(
      lenum, numcomp, gpuData, rmax, dampingType, elementCount,
      elementMapper);
}

void psolve_launch_damping_bkt2_fused(
    int grid, int block, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, double rmax,
    double coef1, double coef2, double expCoef0,
    double coef3, double coef4, double expCoef1,
    int32_t elementCount, const int32_t* elementMapper, int saveVectors) {
  kernelDampingBkt2Fused<<<grid, block, 0, stream>>>(
      lenum, gpuData, rmax, coef1, coef2, expCoef0, coef3, coef4,
      expCoef1, elementCount, elementMapper, saveVectors);
}

void psolve_launch_damping_force(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, int32_t* elementMapper) {
  kernelDampingCalcLocal<<<grid, block, shared, stream>>>(
      lenum, gpuData, elementCount, elementMapper);
}

void psolve_launch_topography_bkt(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, const int32_t* elementMapper,
    const topo_bkt_const_t* constants, double dt2) {
  kernelTopoBktForceEffective<<<grid, block, shared, stream>>>(
      lenum, gpuData, elementCount, elementMapper, constants, dt2);
}

void psolve_launch_topography_effective(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t elementCount, const int32_t* elementMapper,
    const topoconstants_t* constants, double dt2) {
  kernelTopoForceEffective<<<grid, block, shared, stream>>>(
      lenum, gpuData, elementCount, elementMapper, constants, dt2);
}

void psolve_launch_acoustic_displacement_coupling(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt2, double dt) {
  kernelAcousticPressureDisplacementCoupling<<<grid, block, shared, stream>>>(
      faceCount, gpuData, dt2, dt);
}

void psolve_launch_acoustic_pressure_force(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, double dt2) {
  kernelAcousticPressureForce<<<grid, block, shared, stream>>>(
      lenum, gpuData, dt2);
}

void psolve_launch_pml_acoustic(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t pmlAcousticCount) {
  kernelPMLAcousticStep<<<grid, block, shared, stream>>>(
      lenum, gpuData, pmlAcousticCount);
}

void psolve_launch_pml_acoustic_gravity_surface(
    int grid, int block, size_t shared, cudaStream_t stream, int32_t lenum,
    gpu_data_t* gpuData, int32_t surfaceCount) {
  kernelPMLAcousticGravitySurfaceStep<<<grid, block, shared, stream>>>(
      lenum, gpuData, surfaceCount);
}

void psolve_launch_acoustic_interface_force(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt2, double dt) {
  kernelAcousticPressureInterfaceForce<<<grid, block, shared, stream>>>(
      faceCount, gpuData, dt2, dt);
}

void psolve_launch_acoustic_pressure_update(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData) {
  kernelAcousticPressureUpdate<<<grid, block, shared, stream>>>(
      nharbored, gpuData);
}

void psolve_launch_poroelastic_pressure_traction(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t elementCount, gpu_data_t* gpuData, double dt2,
    double biotAlpha) {
  kernelPoroelasticPressureTraction<<<grid, block, shared, stream>>>(
      elementCount, gpuData, dt2, biotAlpha);
}

void psolve_launch_acoustic_open_pore_target(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t faceCount, gpu_data_t* gpuData, double dt) {
  kernelAcousticPoroelasticOpenPoreTarget<<<grid, block, shared, stream>>>(
      faceCount, gpuData, dt);
}

void psolve_launch_poroelastic_pressure_force(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t elementCount, gpu_data_t* gpuData, double dt, double mobilityX,
    double mobilityY, double mobilityZ, int haveDiffusion, double biotAlpha) {
  kernelPoroelasticPressureForce<<<grid, block, shared, stream>>>(
      elementCount, gpuData, dt, mobilityX, mobilityY, mobilityZ,
      haveDiffusion, biotAlpha);
}

void psolve_launch_poroelastic_pressure_update(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData, int useDrainMass,
    int useOpenPore) {
  kernelPoroelasticPressureUpdate<<<grid, block, shared, stream>>>(
      nharbored, gpuData, useDrainMass, useOpenPore);
}

void psolve_launch_plane_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t pointCount, const int32_t* pointNodes, const double* pointCoords,
    const int32_t* pointModes, const fvector_t* tm1, double* pointOutput) {
  kernelPlaneDisplacementInterpolate<<<grid, block, shared, stream>>>(
      pointCount, pointNodes, pointCoords, pointModes, tm1, pointOutput);
}

void psolve_launch_acoustic_station_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const double* stationH,
    const solver_float* pressureTm1, const solver_float* pressureTm2,
    int32_t velocityPotential, double inverseDt, double* stationOutput) {
  kernelAcousticStationInterpolate<<<grid, block, shared, stream>>>(
      stationCount, stationNodes, stationCoords, stationH, pressureTm1,
      pressureTm2, velocityPotential, inverseDt, stationOutput);
}

void psolve_launch_pore_pressure_station_interpolate(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t stationCount, const int32_t* stationNodes,
    const double* stationCoords, const solver_float* porePressureTm1,
    const solver_float* porePressureMass, double* stationOutput) {
  kernelPorePressureStationInterpolate<<<grid, block, shared, stream>>>(
      stationCount, stationNodes, stationCoords, porePressureTm1,
      porePressureMass, stationOutput);
}

void psolve_launch_displacement(
    int grid, int block, size_t shared, cudaStream_t stream,
    int32_t nharbored, gpu_data_t* gpuData, noyesflag_t printAccel) {
  kernelDispCalc<<<grid, block, shared, stream>>>(
      nharbored, gpuData, printAccel);
}
