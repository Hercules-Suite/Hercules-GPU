#ifndef HERCULES_POROELASTIC_DRAIN_CUH
#define HERCULES_POROELASTIC_DRAIN_CUH

/* Optional diagnostics only. For consistent storage, C here is its assembled
 * row sum: summing C_i*dp_i equals 1^T C dp, not the pointwise residual.
 * rhs is already F_step - 2 B p_old, with B=dt R/2. */
__device__ static inline void biot_drain_record(gpu_data_t* data,int n,
    double oldp,double dp,double rhs) {
  if (!data->biotDrainDiagnosticsDevice || !data->biotDrainOwnedDevice[n]) return;
  const double C=data->porePressureMassDevice[n];
  const double B=data->biotDrainCapacityDevice[n],mid=oldp+0.5*dp;
  const double values[6]={C*dp+B*dp-rhs,
      fabs(C*dp)+fabs(B*dp)+fabs(rhs),2*B*mid,2*B*mid*mid,
      C*dp,rhs+2*B*oldp};
  for (int j=0;j<6;++j) atomicAdd(data->biotDrainDiagnosticsDevice+j,values[j]);
}

#endif
