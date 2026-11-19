#ifndef HERCULES_POROELASTIC_MASS_H
#define HERCULES_POROELASTIC_MASS_H

/* Workspace for the optional projected, distributed consistent-mass solve.
 * Buffers contain contiguous u[3*N], w[3*N], p[N] segments. Pressure and
 * mechanics have disjoint segments, including reduction slots. The optional
 * batched solve shares scheduling, never scalar/mechanical Krylov state.
 * No allocation in a time step. */
struct biot_mass_t {
  int nodes, elements, full_nodes, diagonal_nodes;
  int32_t *node, *index, *diagonal; // compact slot -> full node, inverse, direct-only nodes
  int32_t *element, *boundary;
  unsigned char *owned; /* bit 0: independent solid; bit 1: independent pore */
  double *rho, *remainder; /* [node][rho,rho_f,m,capacity] lumped remainder */
  double rho_f, fluid_mass, dt;
  fvector_t* solid_damping; /* physical porous-cell mass damping/dashpots, kg/s */
  solver_float *x, *r, *z, *direction, *product, *rhs;
  double *dots;
  solver_float *precondition, *applied, *history[4][3];
  solver_float *history_product[4][3], *history_preconditioned_product[4][3];
  int rows, uniform_rows;
  int32_t *row_node, *neighbor;
  unsigned char *weight;
  double *row_mass_weight, *row_storage_weight; /* weighted CPML-adjacent rows */
  double row_volume, row_rho, row_capacity;
  int owner_rows;
  solver_float* remote_values;
};

void psolve_biot_mass_elements(cudaStream_t, gpu_data_t*, biot_mass_t, int,
    const solver_float*, solver_float*);
void psolve_biot_mass_operator(cudaStream_t, gpu_data_t*, biot_mass_t, int, int,
    const solver_float*, solver_float*);
void psolve_biot_mass_nodes(cudaStream_t, gpu_data_t*, biot_mass_t, int, int,
    double, double, const solver_float*, solver_float*);
void psolve_biot_mass_dot(cudaStream_t, biot_mass_t, int,
    const solver_float*, const solver_float*);
void psolve_biot_mass_dot_pair(cudaStream_t, biot_mass_t, int,
    const solver_float*, const solver_float*, const solver_float*);
void psolve_biot_mass_finish_product(cudaStream_t, gpu_data_t*, biot_mass_t,
    int, double, const solver_float*, const solver_float*, solver_float*,
    int, const solver_float*, int);
void psolve_biot_mass_local(cudaStream_t, gpu_data_t*, biot_mass_t,
    int, int, double, double, const solver_float*, solver_float*, int);
void psolve_biot_mass_diagonal(cudaStream_t, gpu_data_t*, biot_mass_t, int, int, solver_float*);
void psolve_biot_mass_scalar_io(cudaStream_t, biot_mass_t, solver_float*, int);
void psolve_biot_mass_state(cudaStream_t, gpu_data_t*, biot_mass_t, int, double);
void psolve_biot_mass_predict(cudaStream_t, biot_mass_t, int, int, int);
void psolve_biot_mass_cached_predict(cudaStream_t, gpu_data_t*, biot_mass_t,
    int, int, int, int);
void psolve_biot_mass_cache_store(cudaStream_t, biot_mass_t, int, int, int);
void psolve_biot_mass_input_pack(cudaStream_t,int,int,int,const int32_t*,
    const solver_float*,solver_float*);

/* Entries refer to existing peer packet buffers. Receive entries are stably
 * grouped by destination node, preserving peer and duplicate-entry order. */
struct biot_mass_packet_entry_t {
  solver_float* packet;
  int32_t node, entry, count;
};
void psolve_biot_mass_pack(cudaStream_t, int, int, int, int, int,
    const biot_mass_packet_entry_t*, const solver_float*);
void psolve_biot_mass_unpack(cudaStream_t, int, int, int, int, int, int,
    const int32_t*, const biot_mass_packet_entry_t*, solver_float*);

#endif
