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

/*
 *
 * psolve.h: Generate an unstructured mesh and solve the linear system
 *           thereof derived and output the results.
 *
 * Input:    material database (cvm etree), physics.in, numerical.in.
 * Output:   mesh database (mesh.e) and 4D output.
 *
 *
 * Notes:
 *   For a history of changes see the ChangeLog file.  Also, if you
 *   perform any changes to the source code, please document the changes
 *   by adding the appropriate entry in the ChangeLog file and a descriptive
 *   message during CVS commit.  Thanks!
 */

#ifndef PSOLVE_H
#define PSOLVE_H

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "octor.h"

#define ERROR       HERC_ERROR

#ifndef DO_DEBUG
# ifdef  DEBUG
#  define DO_DEBUG        1
# else
#  define DO_DEBUG        0
# endif /* DEBUG */
#endif /* DO_DEBUG */

/* Maximum simultaneous CUDA streams */
#define MAX_CUDA_STREAM       2
#define CUDA_STREAM_MAIN      0
#define CUDA_STREAM_SECONDARY 1

/* Maximum number of kernels */
#define MAX_CUDA_KERNEL             5
#define CUDA_KERNEL_STIFFNESS_FORCE 0
#define CUDA_KERNEL_DAMPING_CONV    1
#define CUDA_KERNEL_DAMPING_FORCE   2
#define CUDA_KERNEL_DISPLACEMENT    3
#define CUDA_KERNEL_STATION_VE_STRESS 4  /* v3 station viscoelastic stress gather */

/* Kernel names */
#define CUDA_KERNEL_NAME_STIFFNESS_FORCE       "kernelStiffnessCalcLocal"
#define CUDA_KERNEL_NAME_DAMPING_CONV          "kernelDampingCalcConv"
#define CUDA_KERNEL_NAME_DAMPING_CONV_BKT2     "kernelDampingCalcConvBkt2"
#define CUDA_KERNEL_NAME_DAMPING_CONV_BKT3     "kernelDampingCalcConvBkt3"
#define CUDA_KERNEL_NAME_DAMPING_FORCE         "kernelDampingCalcLocal"
#define CUDA_KERNEL_NAME_DISPLACEMENT          "kernelDispCalc"
#define CUDA_KERNEL_NAME_STATION_VE_STRESS     "kernelStationVEStressGather"

/* Maximum simultaneous CUDA events */
#define MAX_CUDA_EVENT                   6
#define CUDA_EVENT_FORCE_UPLOAD_DONE     0
#define CUDA_EVENT_FORCE_STIFFNESS_DONE  1
#define CUDA_EVENT_FORCE_DAMPING_CONV_DONE 2
#define CUDA_EVENT_FORCE_DAMPING_LOCAL_DONE 3
#define CUDA_EVENT_DISP_KERNEL_DONE      4
#define CUDA_EVENT_FORCE_KERNEL_DONE     5

/* GPU device specifications */
typedef struct {
  int32_t device;
  int32_t max_threads;
  int32_t max_block_dim[3];
  int32_t max_grid_dim[3];
  int32_t regs_per_block;
  int32_t shared_per_block;
  int32_t warp_size;
  int32_t warp_allocation_size;
  int32_t register_allocation_size;

  /* Variables for tracking performance metrics */
  int32_t numdevices;
  int64_t numbytespci;
  int64_t numflops;
  int64_t numbytes;
} gpu_spec_t;


/* GPU kernel specifications */
typedef struct {
  char name[32];
  char *handle;
  int gridsize;
  int blocksize;
  int sharedmem;
} gpu_kernel_t;

extern MPI_Comm comm_solver;
extern MPI_Comm comm_IO;

/* Sets solver to use single or double precision in main variables */
#ifdef  SINGLE_PRECISION_SOLVER
   typedef float solver_float;
#else
   typedef double solver_float;
#endif

typedef int32_t     local_id_t;		/**< Local element and node id type  */
typedef local_id_t  lnid_t;		/**< Local node id type		     */
typedef local_id_t  leid_t;		/**< Local element id type	     */
typedef int64_t     global_id_t;	/**< Global element and node id type */
typedef global_id_t gnid_t;		/**< Global node id type	     */
typedef global_id_t geid_t;		/**< Global element id type	     */


struct corner_ref_t {
    lnid_t corner_lnid[8];              /**< Corners local node ids */
    double min_value, max_value;
};

typedef struct corner_ref_t corner_ref_t;

#define ACOUSTIC_FACE_X_MINUS 0x01u
#define ACOUSTIC_FACE_X_PLUS  0x02u
#define ACOUSTIC_FACE_Y_MINUS 0x04u
#define ACOUSTIC_FACE_Y_PLUS  0x08u
#define ACOUSTIC_FACE_Z_MINUS 0x10u
#define ACOUSTIC_FACE_Z_PLUS  0x20u
#define ACOUSTIC_INTERFACE_PML_ACOUSTIC_SIDE 0x01u
#define ACOUSTIC_INTERFACE_PML_ELASTIC_SIDE  0x02u
#define ACOUSTIC_INTERFACE_PROJECTION_POINT     0
#define ACOUSTIC_INTERFACE_PROJECTION_FACE_MASS 1
#define ACOUSTIC_VS_EPS 1.0e-6

/**
 * Mesh database (element) record payload.
 */
struct mdata_t {
    int64_t nid[8];
    float edgesize, Vp, Vs, rho,
          a0_shear, a1_shear, a2_shear,
          a0_kappa, a1_kappa, a2_kappa,
          g0_shear, g1_shear, g2_shear,
          g0_kappa, g1_kappa, g2_kappa,
          b_shear, b_kappa, sigma_0;
    int   topoBkt;
    int32_t topo_eindex;
    double nl_id;
};

typedef struct mdata_t mdata_t;


/**
 * Mesh element data fields
 */
typedef struct edata_t {
    float edgesize, Vp, Vs, rho,
		  a0_shear, a1_shear, a2_shear,
		  a0_kappa, a1_kappa, a2_kappa,
		  g0_shear, g1_shear, g2_shear,
		  g0_kappa, g1_kappa, g2_kappa,
		  b_shear, b_kappa, sigma_0;
    int     topoBkt;
    int32_t topo_eindex;
    double  nl_id;

} edata_t;

/* CPML per-element face mask bit definitions (over the 5 absorbing faces).
 * z_min is the free surface (never PML). Multiple bits set => corner/edge
 * PML element with multi-axis damping. Stored at runtime in a per-element
 * uint8_t array (Global.pmlFaceMask), not persisted in edata_t/mdata_t. */
#define PML_FACE_X_MIN 0x01u  /* south  */
#define PML_FACE_X_MAX 0x02u  /* north  */
#define PML_FACE_Y_MIN 0x04u  /* west   */
#define PML_FACE_Y_MAX 0x08u  /* east   */
#define PML_FACE_Z_MAX 0x10u  /* bottom */
#define PML_FACE_ANY   0x1Fu

/** 3-ary double vector */
struct fvector_t {
    solver_float f[3];
};

typedef struct fvector_t fvector_t;


union df3_t {
    struct {
        double x, y, z;
    };
    double f[3];
};


/**
 * out_hdr_t: 4D output file header.
 */
struct out_hdr_t {
    /** file type string identifier: "Hercules 4D output vnnn"  */
    char    file_type_str[29];
    int8_t  format_version;	/**< File format version		*/
    int8_t  endiannes;		/**< File endianess: 0=little, 1=big	*/

    /** Identifier of the platform where the file was generated */
    int8_t   platform_id;
    unsigned char ufid[16];	/**< "Unique" file identifier		*/
    int64_t  total_nodes;	/**< Node count.			*/
    int32_t  output_steps;	/**< Number of output time steps.	*/

    /** Number of components per (node) record, e.g., 1 vs. 3 */
    int32_t  scalar_count;
    int8_t   scalar_size;	/**< size (in bytes) of each scalar value  */

    /**<
     * Type of scalar, it can take one of the following values
     * - INVALID:		 0
     * - FLOAT32 (float):	 1
     * - FLOAT64 (double):	 2
     * - FLOAT128 (long double): 3
     * - INT8:			 4
     * - UINT8:			 5
     * - INT16:			 6
     * - UINT16:		 7
     * - INT32:			 8
     * - UINT32:		 9
     * - INT64:			 10
     * - UINT64:		 11
     */
    int8_t  scalar_type;

    /**<
     * scalar class, it can take one of the following values:
     * - INVALID:		 0
     * - FLOAT_CLASS		 1
     * - INT_CLASS		 2
     */
    int8_t  scalar_class;

    int8_t  quantity_type;	/* 0: unknown, 1: displacement, 2: velocity */


    /** Mesh parameters: extent of the simulated region in meters */
    double  domain_x, domain_y, domain_z;

    /**
     * Mesh parameter: tick size, factor for converting from domain units (m)
     * to etree units.
     */
    double  mesh_ticksize;

    /** Simulation parameter: delta t (time) */
    double  delta_t;

    /** Mesh parameter: total number of elements */
    int64_t total_elements;

    /** 4D Output parameter: how often is an output time step written out */
    int32_t output_rate;

    /** Simulation parameter: total number of simulation time steps */
    int32_t total_time_steps;

    int64_t generation_date;	/**< Time in seconds since the epoch	*/
};

typedef struct out_hdr_t out_hdr_t;


/*---------------- Solver data structures -----------------------------*/

/**
 * Constants initialized element structure.
 */
struct e_t {
    double c1, c2, c3, c4;
};

typedef struct e_t e_t;


/**
 * Constants initialized node structure.
 *
 * This structure has been changed as a result of the new algorithm
 * for the terms involved with the damping in the solution for the
 * next time step.
 */
struct n_t {
    solver_float mass_simple;
    solver_float mass2_minusaM[3];
    solver_float mass_minusaM[3];
};

typedef struct n_t n_t;

/**
 * Data type to represent stiffness matrices
 */
struct fmatrix_t {
    double f[3][3];
};

typedef struct fmatrix_t fmatrix_t;


/* Solver computation and communication routines */

typedef struct messenger_t messenger_t;

/**
 * A messenger keeps track of the data exchange.
 */
struct messenger_t {
    int32_t  procid;        /**< Remote PE id.     */

    int32_t  outsize;       /**< Outgoing record size. */
    int32_t  insize;        /**< Incoming record size. */

    void*    outdata;
    void*    indata;

    int32_t  nodecount;     /**< Number of mesh nodes involved.         */
    int32_t  nidx;          /**< Current index into the mapping table.  */
    int32_t* mapping;       /**< Array of local ids for involved nodes. */
    int32_t* mapping_itemoffset; /**< Cached mapping * itemsperentry. */
    int32_t  mapping_itemsperentry;
    int32_t  mapping_has_duplicates;
    int32_t  mapping_dup_reported;

    messenger_t* next;
};


/**
 * Communication schedule structure.
 */
struct schedule_t {
    /**
     * Number of PEs this PE contributes to or retrieves data from because
     * this PE harbors their nodes.
     */
    int32_t       c_count;
    messenger_t*  first_c;
    messenger_t** messenger_c;  /**< Fast lookup table to build c-list. */
    MPI_Request*  irecvreqs_c;  /**< control for non-blocking receives. */
    MPI_Status*   irecvstats_c;
    int*          waitsome_idxs_c; /**< scratch indices for MPI_Waitsome. */
    MPI_Request*  isendreqs_c;  /**< control for non-blocking sends.    */
    MPI_Status*   isendstats_c;

    /** Number of processors who share nodes owned by me. */
    int32_t       s_count;
    messenger_t*  first_s;
    messenger_t** messenger_s;  /* Fast lookup table to build s-list.      */
    MPI_Request*  irecvreqs_s;  /* controls for non-blocking MPI receives. */
    MPI_Status*   irecvstats_s;
    int*          waitsome_idxs_s; /* scratch indices for MPI_Waitsome.     */
    MPI_Request*  isendreqs_s;  /* controls for non-blocking MPI sends.    */
    MPI_Status*   isendstats_s;
};

typedef struct schedule_t schedule_t;


/* GPU copy of application data */
typedef struct gpu_data_t {
    int32_t      nharbored; // Number of harbored nodes
    int32_t      lenum; // Number of elements

    //elem_t*      elemTableDevice; // GPU copy of elemTable data structure
    //e_t*         eTableDevice; // GPU copy of eTable data structure
    //n_t*         nTableDevice; // GPU copy of nTable data structure
    fvector_t*   tm1Device; // GPU nodal tm1 displacements
    fvector_t*   tm2Device; // GPU nodal tm2 displacements
    fvector_t*   tm3Device; // GPU nodal tm3 displacements
    fvector_t*   forceDevice; // GPU nodal forces

    solver_float* pressureTm1Device;   // GPU nodal acoustic pressure at t
    solver_float* pressureTm2Device;   // GPU nodal acoustic pressure at t-dt/new
    solver_float* pressureForceDevice; // GPU nodal acoustic pressure RHS
    solver_float* pressureMassDevice;  // GPU lumped acoustic pressure mass
    solver_float* pressureDampDevice;  // GPU nodal acoustic CPML damping C
                                       // (Step B; NULL/0 = no damping ->
                                       // bit-identical to the legacy update)

    fvector_t*   conv_shear_1Device;  /* Approximate Convolution Calculation */
    fvector_t*   conv_shear_2Device;
    fvector_t*   conv_shear_3Device;
    fvector_t*   conv_kappa_1Device;
    fvector_t*   conv_kappa_2Device;
    fvector_t*   conv_kappa_3Device;

    double*      c1ArrayDevice;
    double*      c2ArrayDevice;
    double*      c3ArrayDevice;
    int32_t*     lnidArrayDevice;
    uint8_t*     acousticElementFlagsArrayDevice;
    uint8_t*     acousticFaceFlagsArrayDevice;
    int32_t      acousticInterfaceFaceCount;
    int32_t*     acousticInterfaceElementIdsDevice;
    uint8_t*     acousticInterfaceFaceIdsDevice;
    double*      acousticInterfaceFaceAreasDevice;
    double*      acousticInterfaceFaceCentersDevice;
    double*      acousticInterfaceFaceNormalsDevice;
    double*      acousticInterfaceFaceHalfWidthsDevice;
    int32_t*     acousticInterfaceAdjacentElementIndexesDevice;
    double*      acousticInterfaceAdjacentCentersDevice;
    double*      acousticInterfaceAdjacentHalfWidthsDevice;
    /* CPML-aware interface coupling correction. All NULL unless
     * abc_type=cpml AND at least one acoustic-elastic interface face has a
     * PML element on either side. Per-face sideMask bit 0 means the
     * acoustic side is in PML; bit 1 means the elastic side is in PML.
     * Each side stores its own face-normal-axis ADE coefficients. The
     * acoustic-side profile is used for the elastic->acoustic source term;
     * the elastic-side profile is used for the acoustic-pressure traction
     * on the solid. This avoids using water Vp/profile coefficients for an
     * elastic PML face, or vice versa. */
    uint8_t*     acousticInterfacePmlSideMaskDevice;
    double*      acousticInterfacePmlAcoustic_aDevice;
    double*      acousticInterfacePmlAcoustic_bDevice;
    double*      acousticInterfacePmlAcoustic_invKappaDevice;
    double*      acousticInterfacePmlElastic_aDevice;
    double*      acousticInterfacePmlElastic_bDevice;
    double*      acousticInterfacePmlElastic_invKappaDevice;
    solver_float* acousticInterfacePsiPressureDevice;
    solver_float* acousticInterfacePsiAccelDevice;
    uint8_t*     pressureFreeSurfaceNodeFlagsDevice;
    double*      acousticPressureCoeffArrayDevice;
    double*      acousticPressureBoundaryCoeffArrayDevice;
    double*      g0_shearArrayDevice;
    double*      g1_shearArrayDevice;
    double*      g2_shearArrayDevice;
    double*      g0_kappaArrayDevice;
    double*      g1_kappaArrayDevice;
    double*      g2_kappaArrayDevice;
    double*      b_shearArrayDevice;
    double*      a0_shearArrayDevice;
    double*      a1_shearArrayDevice;
    double*      a2_shearArrayDevice;
    double*      b_kappaArrayDevice;
    double*      a0_kappaArrayDevice;
    double*      a1_kappaArrayDevice;
    double*      a2_kappaArrayDevice;

    solver_float *mass_simpleArrayDevice;
    solver_float *mass2_minusaMArrayDevice[3];
    solver_float *mass_minusaMArrayDevice[3];

    fvector_t*   shearVectorDevice;
    fvector_t*   kappaVectorDevice;
    //edata_t*     matPropsDevice;

    /* CPML state on device, partitioned by material type. Each PML element
     * sits in exactly one mapper (elastic or acoustic); coefficient and
     * memory-variable arrays are sized by the corresponding per-type count.
     *
     * Profile coefficients (a_j, b_j, 1/kappa_j) are precomputed per element
     * from the local edata (Vp, edgesize) and the per-face damping profile
     * (d, kappa, alpha):
     *   pml_b_j      = exp(-(d_j/kappa_j + alpha_j) * dt);
     *   pml_a_j      = (d_j / (kappa_j * (d_j + kappa_j * alpha_j)))*(pml_b_j-1);
     *   pml_invKappa = 1 / kappa_j.
     * For axes outside the face mask the identity (b=1, a=0, invKappa=1) is
     * stored so the recurrence ψ ← b·ψ + a·∂u acts as a no-op on that axis.
     *
     * Elastic psi has 9 components: psi[i][j] is ψ_{i}^{(j)} = j-axis
     * convolution of ∂_j u_i. Acoustic psi has 3 components: psi[j] is the
     * j-axis convolution of ∂_j p (scalar pressure field).
     */
    /* --- Elastic PML state ---
     *
     * a_eff_j  = a_j / (4h),     storage of pml_elastic_a_jArrayDevice
     * invK_eff_j = (1/kappa_j) / (4h), storage of pml_elastic_invKappa_jArrayDevice
     *   (pre-scaling by 1/(4h) absorbs the hex8 1-pt B-matrix factor so the
     *    kernel just multiplies by the raw signed displacement sum.)
     *
     * muCoeff     = dt^2 * h^2 * mu / 4
     * lambdaCoeff = dt^2 * h^2 * lambda / 4
     *   (these combine the Galerkin volume factor h^2/4 with material and
     *    time-step constants so the kernel only needs one multiply per term.)
     */
    int32_t       pmlElasticElemCount;
    int32_t*      pmlElasticElemMapperDevice;
    double*       pml_elastic_b_xArrayDevice;
    double*       pml_elastic_b_yArrayDevice;
    double*       pml_elastic_b_zArrayDevice;
    double*       pml_elastic_a_xArrayDevice;        /* raw a (native norm) */
    double*       pml_elastic_a_yArrayDevice;
    double*       pml_elastic_a_zArrayDevice;
    double*       pml_elastic_invKappa_xArrayDevice; /* raw 1/kappa */
    double*       pml_elastic_invKappa_yArrayDevice;
    double*       pml_elastic_invKappa_zArrayDevice;
    double*       pml_elastic_muCoeffArrayDevice;    /* unused (cleanup) */
    double*       pml_elastic_lambdaCoeffArrayDevice;/* unused (cleanup) */
    solver_float* pml_elastic_psiArrayDevice[3][3];  /* [i_disp][j_axis] */
    /* Stage-2 viscoelastic CPML: separate gradient-side psi memory for the
     * BKT anelastic-displacement vectors (shear=deviatoric, kappa=bulk),
     * carried through the SAME stretch (a/b/invKappa) as the elastic psi.
     * Allocated/used only when cpmlPmlMaterialDamping==1. */
    solver_float* pml_elastic_psiShearArrayDevice[3][3];
    solver_float* pml_elastic_psiKappaArrayDevice[3][3];
    /* Divergence-side CPML, correct 2nd-order displacement form
     * (Matzen 2011; Xie et al. 2014). Multiplying the momentum equation
     * by s_j turns the divergence-side stretch into, per damping axis a:
     *   rho*d_a du_i/dt           (mass-prop. damping, in the mass arrays)
     *   + rho*d_a*alpha_a u_i     (CFS restoring)
     *   - rho*d_a*alpha_a^2 zeta_i^(a)   (CFS displacement memory)
     * The CFS restoring and memory are evaluated together in
     * kernelPMLElasticStep from the element-centroid displacement
     * ubar_i = (1/8) sum_n u_i^(n) so they cancel exactly at DC; the
     * per-node net CFS force is
     *   cfs_i = restoreCoeff*ubar_i - sum_a memForce_a*zeta_i^(a)
     * zeta_i^(a) advances by the exact exponential integrator
     *   zeta_new   = zetaDecay_a*zeta_old + zetaGain_a*ubar_i
     *   zetaDecay_a = exp(-alpha_a*dt)
     *   zetaGain_a  = (1 - exp(-alpha_a*dt)) / alpha_a
     *   memForce_a  = dt^2 * (rho*h^3/8) * d_a * alpha_a^2
     * divPsiArrayDevice is the zeta_i^(a) state (9 per element). */
    solver_float* pml_elastic_divPsiArrayDevice[3][3];/* zeta_i^(a) state */
    double*       pml_elastic_zetaDecayArrayDevice[3];/* exp(-alpha_a dt) */
    double*       pml_elastic_zetaGainArrayDevice[3]; /* (1-decay)/alpha_a */
    double*       pml_elastic_memForceArrayDevice[3]; /* dt^2 Mn d_a a_a^2 */
    /* CFS restoring coeff = dt^2 * M_node * sum_a(d_a*alpha_a). The
     * restoring +rho*d*alpha*u and the memory -rho*d*alpha^2*zeta are
     * BOTH evaluated in-kernel from the element-centroid displacement so
     * they cancel exactly at DC discretely (a per-node restoring vs a
     * centroid zeta did NOT cancel and blew up at corners in v10). */
    double*       pml_elastic_restoreCoeffArrayDevice;

    /* --- Acoustic PML state ---
     *
     * Same pre-scaling for a and invKappa as elastic above.
     * forceCoeff = dt^2 * h^2 / (4 * rho)
     *   (combines volume integration and 1/rho for the acoustic Laplacian.)
     */
    int32_t       pmlAcousticElemCount;
    int32_t*      pmlAcousticElemMapperDevice;
    double*       pml_acoustic_b_xArrayDevice;
    double*       pml_acoustic_b_yArrayDevice;
    double*       pml_acoustic_b_zArrayDevice;
    double*       pml_acoustic_a_xArrayDevice;        /* a / (4h) */
    double*       pml_acoustic_a_yArrayDevice;
    double*       pml_acoustic_a_zArrayDevice;
    double*       pml_acoustic_invKappa_xArrayDevice; /* (1/kappa) / (4h) */
    double*       pml_acoustic_invKappa_yArrayDevice;
    double*       pml_acoustic_invKappa_zArrayDevice;
    double*       pml_acoustic_forceCoeffArrayDevice;
    solver_float* pml_acoustic_psiArrayDevice[3];    /* [j_axis], scalar field */
    /* Step B -- acoustic divergence-side CFS-CPML (scalar mirror of the
     * validated elastic two-sided). zetaP is the pressure-convolution
     * memory (3 per element, one per stretch axis). restoreCoeff/memForce
     * use the acoustic nodal mass (volume/(8*bulk)) in place of rho*h^3/8.
     * The matched kappa-inertia and the mass-prop damping are applied to
     * the nodal pressure mass / pressureDampDevice respectively. */
    double*       pml_acoustic_zetaDecayArrayDevice[3];/* exp(-alpha_a dt) */
    double*       pml_acoustic_zetaGainArrayDevice[3]; /* (1-decay)/alpha_a */
    double*       pml_acoustic_memForceArrayDevice[3]; /* dt^2 Mp d_a a_a^2 */
    double*       pml_acoustic_restoreCoeffArrayDevice;/* dt^2 Mp sum d_a a_a */
    solver_float* pml_acoustic_zetaPArrayDevice[3];    /* [j_axis] zeta_p */

    /* Per-element (sized [lenum]) PML flag. NULL when abc_type != cpml; when
     * non-NULL, has value 1 for any element that lives in either PML pipeline
     * (elastic or acoustic) and 0 otherwise. Lets the existing non-PML
     * kernels (kernelStiffnessCalcLocal, kernelAcousticPressureForce, BKT
     * damping kernels) early-exit on PML elements so we don't double-count. */
    uint8_t*      pmlFlagArrayDevice;

    /* Per-node (sized [nharbored]) Dirichlet clamp mask for the outermost
     * PML face. Set to 1 only when cpml_clamp_outer = yes. Nodes that lie
     * on a domain-extent face touched by a PML element get their
     * displacement forced to zero each step in kernelDispCalc. z_min is the
     * free surface and is never clamped. NULL otherwise. */
    uint8_t*      pmlClampNodeMaskDevice;

    /* Elastic PML stiffness mode (set from Param.theCpmlNativeStiffness
     * at CPML setup). 1 = native effective stiffness + CPML correction
     * (kernelStiffnessCalcLocal processes PML elements; kernelPMLElastic-
     * Step adds -(nodalPml-nodalRaw)). 0 = standalone 1-pt reduced
     * (kernelStiffnessCalcLocal skips PML; kernelPMLElasticStep adds the
     * full -nodalPml). Irrelevant in dashpot mode. */
    int           pmlNativeStiffness;

    /* Acoustic PML stiffness mode (set from
     * Param.theCpmlNativeStiffnessAcoustic at CPML setup). Independent of
     * pmlNativeStiffness so the acoustic PML can run native+correction
     * while the elastic PML stays reduced (or vice versa). 1 = native
     * (kernelAcousticPressureForce processes PML acoustic elements;
     * kernelPMLAcousticStep adds -(P_pml-P_raw)+cfs). 0 = standalone
     * 1-pt reduced (kernelAcousticPressureForce skips PML acoustic;
     * kernelPMLAcousticStep adds full -P_pml+cfs). Irrelevant without
     * acoustic CPML. */
    int           pmlNativeStiffnessAcoustic;
    /* Opt-in material damping inside the elastic PML zone (viscoelastic
     * CFS-CPML). 0 none (default; 5d91b09 skip, bit-identical), 1 bkt,
     * 2 rayleigh. Consumed by the CPML kernels from Stage 2 on. */
    int           cpmlPmlMaterialDamping;

    /* Velocity-potential acoustic-elastic coupling. 0 = pressure;
     * 1 = the acoustic buffers hold a velocity potential phi, the interface
     * exchanges solid normal VELOCITY (1st time-difference, ~omega) and
     * the solid feels p = -dphi/dt, and station pressure output = -dphi/dt. */
    int           acousticVelocityPotential;
    /* Acoustic interface projection. 0 = legacy one-point fragment coupling;
     * 1 = 2x2 face-mass projection over each interface fragment. */
    int           acousticInterfaceProjectionMode;
    double        acousticStiffnessDampingBeta; /* Kelvin-Voigt beta; 0=lossless */

    /* v3 station viscoelastic stress gather (print_station_strains /
     * print_station_stresses under BKT-family damping). Populated only
     * when the v3 path is enabled; otherwise all NULL / zero. The
     * gather kernel computes per-station tensor strain from the
     * modified-displacement shearVectorDevice / kappaVectorDevice
     * (deviatoric channel / volumetric channel) and forms the
     * viscoelastic stress sigma_visco = 2*mu*dev(eps^s) + K*tr(eps^k)*I
     * (K = lambda + 2*mu/3); writes 12 doubles (6 eps^s + 6 sigma)
     * per station into stationVEStressArrayDevice for a single
     * cudaMemcpyAsync back to host each station-output step. */
    int32_t       stationGatherCount;        /* # elastic non-acoustic stations on this rank */
    int32_t*      stationGatherEindexDevice; /* owning element id, per station */
    double*       stationGatherCsiDevice;    /* local coord x (cube case), per station */
    double*       stationGatherEtaDevice;
    double*       stationGatherDzetaDevice;
    double*       stationGatherHDevice;      /* owning cube edge size, per station */
    double*       stationGatherMuDevice;     /* per-station Lame mu */
    double*       stationGatherLambdaDevice; /* per-station Lame lambda */
    int*          stationGatherTethIDevice;  /* tet index (0..4) or -1 for regular cube */
    int*          stationGatherCubePartDevice; /* 1 or 2 (VT cut pattern); 0 ignored */
    double*       stationVEStressArrayDevice;  /* 12 * stationGatherCount: [eps,sig]*6 */

} gpu_data_t;


/**
 * Solver data structure.
 */
struct solver_t {
    e_t *eTable;            /* Element computation-invariant table */
    n_t *nTable;            /* Node computation-invariant table */

    fvector_t *tm1;         /* Displacements at timestep t - 1 */
    fvector_t *tm2;         /* Displacements at timestep t - 2 */
    fvector_t *force;       /* Force accumulation at timestep t */
};

typedef struct solver_t solver_t;


/**
 * \todo Document this struct.
 */
struct mysolver_t {
    e_t* eTable;              /* Element computation-invariant table  */
    n_t* nTable;              /* Node computation-invariant table     */

    /* \TODO Explain better what is tm1,2,3 */
    fvector_t* tm1;           /* Displacements at time-step t - 1      */
    fvector_t* tm2;           /* Displacements at time-step t - 2      */
    fvector_t* tm3;           /* Displacements at time-step t - 3      */
    fvector_t* force;         /* Force accumulation at time-step t     */

    solver_float* pressure_tm1;   /* Acoustic pressure at time-step t - 1 */
    solver_float* pressure_tm2;   /* Acoustic pressure at time-step t - 2 */
    solver_float* pressure_force; /* Acoustic pressure RHS accumulation */
    solver_float* pressure_mass;  /* Acoustic pressure lumped mass */
    solver_float* pressure_damp;  /* Acoustic CPML nodal damping C (Step B;
                                   * 0 everywhere unless acoustic PML) */

    schedule_t* dn_sched;     /* Dangling node communication schedule */
    schedule_t* an_sched;     /* Anchored node communication schedule */

    fvector_t* conv_shear_1;  /* Approximate convolution Calculation */
    fvector_t* conv_shear_2;
    fvector_t* conv_shear_3;

    fvector_t* conv_kappa_1;
    fvector_t* conv_kappa_2;
    fvector_t* conv_kappa_3;

    /* CPML host-side state (used for checkpoint I/O only). NULL when
     * abc_type != cpml or the respective count is zero. The psi host mirrors
     * shadow the device arrays of the same name; the elem_list arrays are
     * host copies of the device mappers so checkpoint writes don't need a
     * device->host pull just to know which eindex each slot belongs to. */
    int32_t       pmlElasticElemCount;
    int32_t*      pmlElasticElemList;
    solver_float* pml_elastic_psi[3][3];

    int32_t       pmlAcousticElemCount;
    int32_t*      pmlAcousticElemList;
    solver_float* pml_acoustic_psi[3];

    /* GPU device data structures - should be declared in a separate struct */
    gpu_spec_t*  gpu_spec; // GPU device specifications
    gpu_kernel_t*  gpu_kernel; // GPU kernel specifications
    gpu_data_t*  gpuData;
    gpu_data_t*  gpuDataDevice;

    cudaStream_t streams[MAX_CUDA_STREAM]; // Streams asynchronous operation
    cudaEvent_t  events[MAX_CUDA_EVENT];   // Events for asynch operation
};

typedef struct mysolver_t mysolver_t;


/* This structures was copied from geometrics.h
 * \todo We need to unify this structure because is used in several
 *       parts of the code and with different names, see above, for
 *       example, fvector_t
 */
struct vector3D_t {
    double x[3];
};

typedef struct vector3D_t vector3D_t;


/* ------------------------------------------------------------------------- *
 *               Output stations data structures
 * ------------------------------------------------------------------------- */

struct station_t {
    int32_t    id, nodestointerpolate[8];
    int32_t    pressure_is_defined;
    int32_t    acoustic_kinematics_ready;
    double*    displacementsX;
    double*    displacementsY;
    double*    displacementsZ;
    double     prev_disp[3];
    double     acoustic_disp[3];
    double     acoustic_vel[3];
    double     acoustic_acc[3];
    double     acoustic_prev_acc[3];
    double     h; /* dorian: added size edge to compute station's strains */
    double     rho;
    double     mu;
    double     lambda;
    int32_t    owningElementIndex; /* eindex of owning octree element (v3 gather) */

    vector3D_t coords;        /* cartesian */
    vector3D_t localcoords;   /* csi, eta, dzeta in (-1,1) x (-1,1) x (-1,1) */
    FILE*      fpoutputfile;
};

typedef struct station_t station_t;



#ifdef __cplusplus
extern "C" {
#endif

extern int  solver_abort (const char* function_name, const char* error_msg,
			  const char* format, ...);
extern void solver_output_seq (void);
extern int  parsetext (FILE* fp, const char* querystring, const char type,
		       void* result);



/**
 * mu_and_lambda: Calculates mu and lambda according to the element values
 *                of Vp, Vs, and Rho and verifies/applies some rules defined
 *                by Jacobo, Leonardo and Ricardo.  It was originally within
 *                solver_init but was moved out because it needs to be used
 *                in other places as well (nonlinear)
 */
void mu_and_lambda(double *theMu, double *theLambda, edata_t *edata, int32_t eindex);

/**
 * Search a point in the domain of the local mesh.
 *
 *   input: coordinates
 *  output: 0 fail 1 success
 */
int32_t search_point( vector3D_t point, octant_t **octant );
int32_t solver_get_octant_material( octant_t *octant, edata_t *material );

/* SPECFEM-style flat/regular PML: returns 1 iff element `eindex` is a
 * VT-topography CPML PML solid element whose body-fitted deformation
 * must be suppressed (abc_type=cpml AND include_topography=yes AND
 * computation_method=VT AND pml_face_mask!=0 AND not an air element).
 * VT body-fits/deforms surface elements (-> curvilinear topo kernel),
 * which mismatches the cube CPML delta; this re-routes them to the
 * regular cube path. FEM topography is a regular Cartesian mesh + air
 * elements (no element deformation), so this is a strict NO-OP for FEM
 * -- and 0 whenever CPML or topography is off, so every flat-surface
 * CPML run (all Stage 1-4 validation) and every FEM/dashpot run is
 * byte-identical. Air elements (Vp=-1) always return 0 (must not enter
 * the solid stiffness mapper). Single source of truth for the
 * linear/topo mapper partition (stiffness.cu), the topo GPU mapper
 * filter and the toponodes_mass gate. g_cpmlBuild.face_mask is filled
 * by solver_init() before stiffness_init()/the topo mapper, so it is
 * valid at all call sites. */
int cpml_pml_flattened_element( int32_t eindex );


#ifdef __cplusplus
}
#endif
#endif /* PSOLVE_H */
