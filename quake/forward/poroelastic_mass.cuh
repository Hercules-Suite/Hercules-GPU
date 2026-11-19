#ifndef HERCULES_POROELASTIC_MASS_CUH
#define HERCULES_POROELASTIC_MASS_CUH

#include "poroelastic_biot.h"
#include "poroelastic_drain.cuh"

__device__ static double bm_surface_capacity(gpu_data_t* data,int n,double scale) {
  return scale && data->biotDrainCapacityDevice?
      scale*data->biotDrainCapacityDevice[n]:0;
}

__device__ static double bm_solid_damping(biot_mass_t m,int n,int d) {
  return m.solid_damping?m.solid_damping[n].f[d]:0;
}

__device__ static size_t bm_index(biot_mass_t m,int n,int c) {
  if (m.index) n=m.index[n];
  return c<3?3*(size_t)n+c:(c<6?3*(size_t)m.nodes+3*(size_t)n+c-3:6*(size_t)m.nodes+n);
}

/* Exact Cartesian Q1 or native retained-volume P1 integral of Ni*Nj.
 * The native VT fractions describe the same approximation used by stiffness. */
__global__ static void bm_elements(gpu_data_t* data,biot_mass_t m,int scalar,int patch,
    const solver_float* in,solver_float* out) {
  const int e=blockIdx.x*blockDim.x+threadIdx.x;
  if (e>=m.elements) return;
  const int pi=m.element[e], ei=data->poroelasticElemMapperDevice[pi];
  const double h=data->poroelasticElemHDevice[pi], volume=h*h*h;
  const int cpml=data->poroCpmlIndexDevice?data->poroCpmlIndexDevice[pi]:-1;
  const double inertia=cpml>=0?data->poroCpmlCoeffDevice[cpml].inertia.mass:1;
  const double capacity=cpml>=0?poro_cpml_capacity(data->poroCpmlCoeffDevice[cpml]):1;
  int ids[8]; double mass[8][8];
  for (int i=0;i<8;++i) ids[i]=data->lnidArrayDevice[ei+i*data->lenum];
  const int ci=data->biotCutIndexDevice?data->biotCutIndexDevice[pi]:-1;
  if (ci>=0) {
    for (int i=0;i<8;++i) for (int j=0;j<8;++j) mass[i][j]=0;
    const biot_cut_t cut=data->biotCutsDevice[ci];
    for (int t=0;t<5;++t) {
      int nodes[4];double grad[4][3];
      const double v=biot_tet_geometry(cut.partition,t,nodes,grad)*volume*cut.fraction[t]/20;
      for (int a=0;a<4;++a) for (int b=0;b<4;++b)
        mass[nodes[a]][nodes[b]]+=v*(patch?(a==b?20:-5):(a==b?2:1));
    }
  }
  for (int i=0;i<8;++i) {
    double sum[7]={};
    for (int j=0;j<8;++j) {
      double v;
      if (ci>=0) v=mass[i][j];
      else {
        v=volume/216;
        for (int d=0;d<3;++d) if (((i^j)&(1<<d))==0) v*=2;
        if (patch) v*=27*((__popc((unsigned)(i^j))&1)?-1:1);
      }
      if (scalar!=0) sum[6]+=v*capacity*in[bm_index(m,ids[j],6)]/data->poroelasticStorageMDevice[pi];
      if (scalar!=1) for (int d=0;d<3;++d) {
        const double u=in[bm_index(m,ids[j],d)], w=in[bm_index(m,ids[j],d+3)];
        sum[d]+=v*inertia*(m.rho[e]*u+m.rho_f*w);
        sum[d+3]+=v*inertia*(m.rho_f*u+m.fluid_mass*w);
      }
    }
    for (int c=scalar==1?6:0;c<(scalar!=0?7:6);++c)
      atomicAdd(out+bm_index(m,ids[i],c),(solver_float)sum[c]);
  }
}

/* Locally assembled Cartesian Q1 rows. Each row has one writer; shared-node
 * partial rows are still accumulated by the established MPI contribution
 * exchange. Opposite-corner signs also give the element-inverse preconditioner
 * H=sum Le Me^-1 Le, without assembling a second stencil. */
template<int scalar,bool weighted=false>
__global__ static void bm_rows(biot_mass_t m,int patch,
    const solver_float* in,solver_float* out) {
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  const int count=weighted?m.rows-m.uniform_rows:m.uniform_rows;
  if (slot>=count) return;
  const int row=slot+(weighted?m.uniform_rows:0);
  double sum[7]={};
  for (int j=0;j<27;++j) {
    const size_t index=(size_t)j*m.rows+row;
    const int weight=m.weight[index];
    if (!weight) continue;
    const int node=m.neighbor[index];
    const size_t wi=(size_t)j*count+slot;
    double factor=weighted?m.row_mass_weight[wi]:weight;
    double capacity=weighted?m.row_storage_weight[wi]:weight;
    if (patch) {
      const int odd=(j%3!=1)+(j/3%3!=1)+(j/9!=1);
      factor*=27*((odd&1)?-1:1);
      capacity*=27*((odd&1)?-1:1);
    }
    #pragma unroll
    for (int c=scalar==1?6:0;c<(scalar!=0?7:6);++c) {
      const int width=scalar==1?1:(scalar==2?7:6);
      const double value=m.owner_rows && node<0?
          m.remote_values[width*(size_t)(-1-node)+(scalar==1?0:c)]:in[bm_index(m,node,c)];
      sum[c]+=(c==6?capacity:factor)*value;
    }
  }
  const int n=m.row_node[row];
  const double factor=m.row_volume/216;
  if (scalar!=0) out[bm_index(m,n,6)]=factor*m.row_capacity*sum[6];
  if (scalar!=1) {
    #pragma unroll
    for (int d=0;d<3;++d) {
    out[bm_index(m,n,d)]=factor*(m.row_rho*sum[d]+m.rho_f*sum[d+3]);
    out[bm_index(m,n,d+3)]=factor*(m.rho_f*sum[d]+m.fluid_mass*sum[d+3]);
    }
  }
}

__device__ static void bm_project(gpu_data_t* data,biot_mass_t m,int n,
    int scalar,double v[7]) {
  if (scalar) {
    if (!(m.owned[n]&2) || !(data->porePressureMassDevice[n]>0)) v[6]=0;
    return;
  }
  if (!(m.owned[n]&1) || !(data->mass_simpleArrayDevice[n]>0) ||
      (data->pmlClampNodeMaskDevice && data->pmlClampNodeMaskDevice[n])) {
    for (int c=0;c<6;++c) v[c]=0;
    return;
  }
  if (m.boundary[n]>=0) {
    const biot_boundary_t bc=data->biotBoundaryDevice[m.boundary[n]];
    double w[3]={};
    for (int a=0;a<bc.modes;++a) {
      double dot=0;for (int d=0;d<3;++d) dot+=bc.basis[a][d]*v[d+3];
      for (int d=0;d<3;++d) w[d]+=bc.basis[a][d]*dot;
    }
    for (int d=0;d<3;++d) v[d+3]=w[d];
  } else for (int d=0;d<3;++d)
    if (!(data->biotMassDevice[n].f[1]>0) || data->biotDragDevice[n].f[d]<0) v[d+3]=0;
}

__device__ static void bm_precondition(gpu_data_t* data,biot_mass_t m,int n,
    int scalar,double drag_dt,double v[7]) {
  bm_project(data,m,n,scalar,v);
  if (scalar) {
    const double C=data->porePressureMassDevice[n]+bm_surface_capacity(data,n,drag_dt);
    if (C>0) v[6]/=C;
    return;
  }
  const double ms=data->mass_simpleArrayDevice[n];
  if (!(ms>0)) return;
  const double mf=data->biotMassDevice[n].f[0],mw=data->biotMassDevice[n].f[1];
  double w[3]={},force[3],solid[3],schur[3];
  for (int d=0;d<3;++d) {
    solid[d]=ms+drag_dt*bm_solid_damping(m,n,d);
    force[d]=v[d+3]-mf/solid[d]*v[d];
    schur[d]=mw+drag_dt*fmax(0.,(double)data->biotDragDevice[n].f[d])-mf*mf/solid[d];
  }
  if (mw>0) {
    if (m.boundary[n]>=0) {
      const biot_boundary_t bc=data->biotBoundaryDevice[m.boundary[n]];
      // The anisotropic solid damping changes the fluid Schur complement;
      // project the full block instead of reusing the Darcy-only eigenvalues.
      double A[2][2]={},f[2]={},z[2]={};
      for (int a=0;a<bc.modes;++a) for (int d=0;d<3;++d) {
        f[a]+=bc.basis[a][d]*force[d];
        for (int b=0;b<bc.modes;++b) A[a][b]+=bc.basis[a][d]*schur[d]*bc.basis[b][d];
      }
      if (bc.modes==1) z[0]=f[0]/A[0][0];
      if (bc.modes==2) {
        const double det=A[0][0]*A[1][1]-A[0][1]*A[1][0];
        z[0]=(A[1][1]*f[0]-A[0][1]*f[1])/det;
        z[1]=(A[0][0]*f[1]-A[1][0]*f[0])/det;
      }
      for (int a=0;a<bc.modes;++a) for (int d=0;d<3;++d) w[d]+=bc.basis[a][d]*z[a];
    } else for (int d=0;d<3;++d) if (data->biotDragDevice[n].f[d]>=0)
      w[d]=force[d]/schur[d];
  }
  for (int d=0;d<3;++d) {v[d]=(v[d]-mf*w[d])/solid[d];v[d+3]=w[d];}
  bm_project(data,m,n,scalar,v);
}

/* op 0: add owned lumped remainder + drag to assembled element product;
 * 1: project/copy; 2: block-Jacobi; 3: out=a*in+b*out; 4: rhs-in.
 * All Krylov vectors are zero on ghosts/constrained slaves until sharing. */
__global__ static void bm_nodes(gpu_data_t* data,biot_mass_t m,int scalar,int op,
    double a,double b,const solver_float* in,solver_float* out) {
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  if (slot>=m.nodes) return;
  const int n=m.node?m.node[slot]:slot;
  double v[7]={};
  const int start=scalar?6:0,end=scalar?7:6;
  for (int c=start;c<end;++c) {
    const size_t i=bm_index(m,n,c);
    if (op==3) v[c]=a*in[i]+b*out[i];
    else if (op==4) v[c]=m.rhs[i]-in[i];
    else v[c]=in[i];
  }
  if (op==0) {
    const double* R=m.remainder+4*(size_t)n;
    if (scalar) v[6]=out[bm_index(m,n,6)]+(R[3]+bm_surface_capacity(data,n,a))*v[6];
    else for (int d=0;d<3;++d) {
      const double u=v[d], w=v[d+3], D=fmax(0.,(double)data->biotDragDevice[n].f[d]);
      v[d]=out[bm_index(m,n,d)]+(R[0]+a*bm_solid_damping(m,n,d))*u+R[1]*w;
      v[d+3]=out[bm_index(m,n,d+3)]+R[1]*u+(R[2]+a*D)*w;
    }
  }
  if (op==2) bm_precondition(data,m,n,scalar,a,v);
  else bm_project(data,m,n,scalar,v);
  for (int c=start;c<end;++c) out[bm_index(m,n,c)]=(solver_float)v[c];
}

__global__ static void bm_dot(biot_mass_t m,int scalar,
    const solver_float* x,const solver_float* y) {
  __shared__ double sum[256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  const int n=slot<m.nodes && m.node?m.node[slot]:slot;
  double v=0;
  if (slot<m.nodes && (m.owned[n]&(scalar?2:1)))
    for (int c=scalar?6:0;c<(scalar?7:6);++c) {
      const size_t i=bm_index(m,n,c);v+=(double)x[i]*y[i];
    }
  sum[threadIdx.x]=v;__syncthreads();
  for (int k=128;k;k>>=1) {
    if (threadIdx.x<k) sum[threadIdx.x]+=sum[threadIdx.x+k];
    __syncthreads();
  }
  if (threadIdx.x==0) atomicAdd(m.dots,sum[0]);
}

__global__ static void bm_dot_pair(biot_mass_t m,int scalar,
    const solver_float* r,const solver_float* z,const solver_float* az) {
  __shared__ double sum[2][256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  const int n=slot<m.nodes && m.node?m.node[slot]:slot;
  double rz=0,zaz=0;
  if (slot<m.nodes && (m.owned[n]&(scalar?2:1)))
    for (int c=scalar?6:0;c<(scalar?7:6);++c) {
      const size_t i=bm_index(m,n,c);
      rz+=(double)r[i]*z[i];zaz+=(double)z[i]*az[i];
    }
  sum[0][threadIdx.x]=rz;sum[1][threadIdx.x]=zaz;__syncthreads();
  for (int k=128;k;k>>=1) {
    if (threadIdx.x<k) for (int j=0;j<2;++j)
      sum[j][threadIdx.x]+=sum[j][threadIdx.x+k];
    __syncthreads();
  }
  if (!threadIdx.x) for (int j=0;j<2;++j) atomicAdd(m.dots+j,sum[j][0]);
}

template<int count>
__device__ static void bm_fused_sum(biot_mass_t m,double a,double b,
    double (&sum)[count?count:1][256]) {
  if (count) {
    sum[0][threadIdx.x]=a;
    if (count==2) sum[1][threadIdx.x]=b;
    __syncthreads();
    for (int k=128;k;k>>=1) {
      if (threadIdx.x<k) for (int j=0;j<count;++j)
        sum[j][threadIdx.x]+=sum[j][threadIdx.x+k];
      __syncthreads();
    }
    if (!threadIdx.x) for (int j=0;j<count;++j) atomicAdd(m.dots+j,sum[j][0]);
  }
}

/* Finish only AFTER MPI contribution assembly. Optionally apply the final
 * nodal inverse of L^-1 H L^-1 and reduce on independent owned DOFs. The
 * explicit solver_float roundings retain the separate-kernel storage points.
 * dot_mode 1: left.out; 2: (left.in,in.out), for single-reduction CG. */
template<bool scalar,int dot_mode>
__global__ static void bm_finish_product(gpu_data_t* data,biot_mass_t m,double drag_dt,
    const solver_float* in,const solver_float* assembled,solver_float* out,
    int precondition,const solver_float* left) {
  __shared__ double sum[dot_mode?dot_mode:1][256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  const int n=slot<m.nodes && m.node?m.node[slot]:slot;
  double a=0,b=0;
  if (slot<m.nodes) {
    double v[7]={};
    const double* R=m.remainder+4*(size_t)n;
    if (scalar) {
      const size_t i=bm_index(m,n,6);
      v[6]=assembled[i]+(R[3]+bm_surface_capacity(data,n,drag_dt))*in[i];
    } else for (int d=0;d<3;++d) {
      const size_t iu=bm_index(m,n,d),iw=bm_index(m,n,d+3);
      const double u=in[iu],w=in[iw],D=fmax(0.,(double)data->biotDragDevice[n].f[d]);
      v[d]=assembled[iu]+(R[0]+drag_dt*bm_solid_damping(m,n,d))*u+R[1]*w;
      v[d+3]=assembled[iw]+R[1]*u+(R[2]+drag_dt*D)*w;
    }
    bm_project(data,m,n,scalar,v);
    if (precondition) {
      for (int c=scalar?6:0;c<(scalar?7:6);++c) v[c]=(solver_float)v[c];
      bm_precondition(data,m,n,scalar,drag_dt,v);
    }
    for (int c=scalar?6:0;c<(scalar?7:6);++c) {
      const size_t i=bm_index(m,n,c);
      const solver_float value=(solver_float)v[c];
      if (dot_mode && (m.owned[n]&(scalar?2:1))) {
        if (dot_mode==1) a+=(double)left[i]*value;
        else {a+=(double)left[i]*in[i];b+=(double)in[i]*value;}
      }
      out[i]=value;
    }
  }
  bm_fused_sum<dot_mode>(m,a,b,sum);
}

/* mode 0: nodal precondition only; 1: x/r update + nodal precondition;
 * 2: true residual + nodal precondition; 3: projected RHS + precondition.
 * For element-inverse PCG this is ONLY the first L^-1 stage. Its H application
 * and final L^-1 stage still follow assembly/communication in finish_product. */
template<bool scalar,int dot>
__global__ static void bm_local(gpu_data_t* data,biot_mass_t m,int mode,
    double alpha,double drag_dt,const solver_float* in,solver_float* out) {
  __shared__ double sum[dot?dot:1][256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  const int n=slot<m.nodes && m.node?m.node[slot]:slot;
  double rz=0;
  if (slot<m.nodes) {
    const int first=scalar?6:0,last=scalar?7:6;
    double v[7]={},r[7]={};
    if (mode==1) {
      for (int c=first;c<last;++c) {
        const size_t i=bm_index(m,n,c);v[c]=alpha*m.direction[i]+m.x[i];
      }
      bm_project(data,m,n,scalar,v);
      for (int c=first;c<last;++c) m.x[bm_index(m,n,c)]=(solver_float)v[c];
    }
    for (int c=first;c<last;++c) {
      const size_t i=bm_index(m,n,c);
      v[c]=mode==1?-alpha*m.product[i]+m.r[i]:
          (mode==2?(double)m.rhs[i]-m.product[i]:(mode==3?m.rhs[i]:in[i]));
    }
    if (mode) bm_project(data,m,n,scalar,v);
    for (int c=first;c<last;++c) {
      r[c]=v[c]=(solver_float)v[c];
      if (mode) m.r[bm_index(m,n,c)]=(solver_float)v[c];
    }
    bm_precondition(data,m,n,scalar,drag_dt,v);
    for (int c=first;c<last;++c) {
      const size_t i=bm_index(m,n,c);const solver_float value=(solver_float)v[c];
      out[i]=value;
      if (dot && (m.owned[n]&(scalar?2:1))) rz+=r[c]*value;
    }
  }
  bm_fused_sum<dot>(m,rz,0,sum);
}

/* op 0: old increments -> x; 1: rhs += a*dt^2*force;
 * 2: rhs -= a*D*x (SDIRK stage 2); 3/4: finish central-difference/SDIRK;
 * 5: pressure update. Ordinary displacement launch already saved u_old in tm3. */
__global__ static void bm_state(gpu_data_t* data,biot_mass_t m,int op,double a) {
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  if (slot>=m.nodes) return;
  const int n=m.node?m.node[slot]:slot;
  if (op==5) {
    data->porePressureTm2Device[n]=data->porePressureTm1Device[n]+m.x[bm_index(m,n,6)];
    biot_drain_record(data,n,data->porePressureTm1Device[n],
        m.x[bm_index(m,n,6)],data->porePressureForceDevice[n]);
    data->porePressureForceDevice[n]=0;
    return;
  }
  double v[7]={};
  for (int d=0;d<3;++d) {
    const size_t iu=bm_index(m,n,d),iw=bm_index(m,n,d+3);
    const double olddu=(double)data->tm1Device[n].f[d]-data->tm3Device[n].f[d];
    const double olddw=(double)data->biotW1Device[n].f[d]-data->biotW2Device[n].f[d];
    if (op==0) {v[d]=olddu;v[d+3]=olddw;}
    else if (op==1) {
      v[d]=m.rhs[iu]+a*(data->forceDevice[n].f[d]-
          (data->mass_simpleArrayDevice[n]-data->mass_minusaMArrayDevice[d][n]-
           m.dt*bm_solid_damping(m,n,d))*olddu);
      v[d+3]=m.rhs[iw]+a*data->biotForceDevice[n].f[d];
    } else if (op==2) {
      v[d]=m.rhs[iu]-a*bm_solid_damping(m,n,d)*m.x[iu];
      v[d+3]=m.rhs[iw]-a*fmax(0.,(double)data->biotDragDevice[n].f[d])*m.x[iw];
    } else {
      v[d]=m.x[iu]+(op==3?olddu:0);
      v[d+3]=m.x[iw]+(op==3?olddw:0);
    }
  }
  bm_project(data,m,n,0,v);
  for (int d=0;d<3;++d) {
    const size_t iu=bm_index(m,n,d),iw=bm_index(m,n,d+3);
    if (op==0) {m.x[iu]=v[d];m.x[iw]=v[d+3];}
    else if (op<=2) {m.rhs[iu]=v[d];m.rhs[iw]=v[d+3];}
    else {
      data->biotW3Device[n].f[d]=data->biotW2Device[n].f[d];
      data->tm2Device[n].f[d]=data->tm1Device[n].f[d]+v[d];
      data->biotW2Device[n].f[d]=data->biotW1Device[n].f[d]+v[d+3];
    }
  }
}

/* Independent nodal blocks have no consistent-cell/constraint coupling.
 * Advance them directly once per physical update, retaining the same nodal
 * damping stages and the full-RHS normalization used by compact PCG. */
__global__ static void bm_diagonal(gpu_data_t* data,biot_mass_t m,int scalar,
    int mode,solver_float* rhs) {
  __shared__ double sum[2][256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  double norm[2]={};
  if (slot<m.diagonal_nodes) {
    const int n=m.diagonal[slot];double r[7]={},z[7]={};
    if (scalar) {
      r[6]=rhs[n];bm_project(data,m,n,1,r);z[6]=r[6];
      bm_precondition(data,m,n,1,mode?1:0,z);z[6]=(solver_float)z[6];
      norm[0]=r[6]*z[6];
      if (mode) {
        data->porePressureTm2Device[n]=data->porePressureTm1Device[n]+z[6];
        biot_drain_record(data,n,data->porePressureTm1Device[n],z[6],rhs[n]);
        data->porePressureForceDevice[n]=0;
      } else rhs[n]=z[6];
    } else {
      const double gamma=1-1/sqrt(2.0),dt=m.dt;
      double old[7]={},force[7]={};
      for (int d=0;d<3;++d) {
        old[d]=(double)data->tm1Device[n].f[d]-data->tm3Device[n].f[d];
        old[d+3]=(double)data->biotW1Device[n].f[d]-data->biotW2Device[n].f[d];
        force[d]=data->forceDevice[n].f[d]-(data->mass_simpleArrayDevice[n]-
            data->mass_minusaMArrayDevice[d][n]-dt*bm_solid_damping(m,n,d))*old[d];
        force[d+3]=data->biotForceDevice[n].f[d];
      }
      if (mode) {
        double v[7];for (int c=0;c<6;++c)v[c]=old[c];
        bm_project(data,m,n,0,v);
        for (int c=0;c<6;++c)v[c]=(solver_float)v[c];
        const double* R=m.remainder+4*(size_t)n;
        for (int d=0;d<3;++d) {r[d]=R[0]*v[d]+R[1]*v[d+3];r[d+3]=R[1]*v[d]+R[2]*v[d+3];}
        bm_project(data,m,n,0,r);
        for (int c=0;c<6;++c)r[c]=(solver_float)r[c];
      }
      for (int c=0;c<6;++c)r[c]+=(mode?gamma:1)*force[c];
      bm_project(data,m,n,0,r);
      for (int c=0;c<6;++c)z[c]=r[c]=(solver_float)r[c];
      bm_precondition(data,m,n,0,mode?gamma*dt:0,z);
      for (int c=0;c<6;++c) {z[c]=(solver_float)z[c];norm[0]+=r[c]*z[c];}
      if (mode) {
        for (int c=0;c<6;++c)r[c]+=(1-gamma)*force[c];
        bm_project(data,m,n,0,r);
        for (int c=0;c<6;++c)r[c]=(solver_float)r[c];
        for (int d=0;d<3;++d) {
          r[d]-=(1-gamma)*dt*bm_solid_damping(m,n,d)*z[d];
          r[d+3]-=(1-gamma)*dt*fmax(0.,(double)data->biotDragDevice[n].f[d])*z[d+3];
        }
        bm_project(data,m,n,0,r);
        for (int c=0;c<6;++c)z[c]=r[c]=(solver_float)r[c];
        bm_precondition(data,m,n,0,gamma*dt,z);
        for (int c=0;c<6;++c) {z[c]=(solver_float)z[c];norm[1]+=r[c]*z[c];}
      } else for (int c=0;c<6;++c)z[c]+=old[c];
      bm_project(data,m,n,0,z);
      for (int d=0;d<3;++d) {
        data->biotW3Device[n].f[d]=data->biotW2Device[n].f[d];
        data->tm2Device[n].f[d]=data->tm1Device[n].f[d]+z[d];
        data->biotW2Device[n].f[d]=data->biotW1Device[n].f[d]+z[d+3];
      }
    }
  }
  m.dots+=2;bm_fused_sum<2>(m,norm[0],norm[1],sum);
}
__global__ static void bm_scalar_io(biot_mass_t m,solver_float* full,int scatter) {
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  if (slot>=m.nodes)return;
  const int n=m.node?m.node[slot]:slot;
  if (scatter) full[n]=m.x[6*(size_t)m.nodes+slot];
  else m.rhs[6*(size_t)m.nodes+slot]=full[n];
}
void psolve_biot_mass_diagonal(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,int mode,solver_float* rhs) {
  cudaMemsetAsync(m.dots+2,0,2*sizeof(double),stream);
  if (m.diagonal_nodes) bm_diagonal<<<(m.diagonal_nodes+255)/256,256,0,stream>>>(data,m,scalar,mode,rhs);
}
void psolve_biot_mass_scalar_io(cudaStream_t stream,biot_mass_t m,solver_float* full,int scatter) {
  if (m.nodes)bm_scalar_io<<<(m.nodes+255)/256,256,0,stream>>>(m,full,scatter);
}

void psolve_biot_mass_elements(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,const solver_float* in,solver_float* out) {
  if (m.elements) bm_elements<<<(m.elements+63)/64,64,0,stream>>>(data,m,scalar,0,in,out);
}
void psolve_biot_mass_operator(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,int patch,const solver_float* in,solver_float* out) {
  if (m.owner_rows || m.neighbor) {
    if (m.uniform_rows) {
      const int blocks=(m.uniform_rows+127)/128;
      if (scalar==2) bm_rows<2><<<blocks,128,0,stream>>>(m,patch,in,out);
      if (scalar==1) bm_rows<1><<<blocks,128,0,stream>>>(m,patch,in,out);
      if (!scalar) bm_rows<0><<<blocks,128,0,stream>>>(m,patch,in,out);
    }
    if (m.rows>m.uniform_rows) {
      const int blocks=(m.rows-m.uniform_rows+127)/128;
      if (scalar==2) bm_rows<2,true><<<blocks,128,0,stream>>>(m,patch,in,out);
      if (scalar==1) bm_rows<1,true><<<blocks,128,0,stream>>>(m,patch,in,out);
      if (!scalar) bm_rows<0,true><<<blocks,128,0,stream>>>(m,patch,in,out);
    }
  } else if (m.elements)
    bm_elements<<<(m.elements+63)/64,64,0,stream>>>(data,m,scalar,patch,in,out);
}
void psolve_biot_mass_nodes(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,int op,double a,double b,const solver_float* in,solver_float* out) {
  if (m.nodes) bm_nodes<<<(m.nodes+255)/256,256,0,stream>>>(data,m,scalar,op,a,b,in,out);
}

/* Owner-row input halo: no field staging on CUDA-aware MPI and no output
 * assembly. Each remote node occupies width contiguous packet components. */
__global__ static void bm_input_pack(int nodes,int scalar,int count,
    const int32_t* map,const solver_float* in,solver_float* out) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i>=count)return;const int n=map[i],width=scalar==1?1:(scalar==2?7:6);
  for(int j=0;j<width;++j) {
    const int c=scalar==1?6:j;
    const size_t k=c<6?(c/3)*3*(size_t)nodes+3*(size_t)n+c%3:6*(size_t)nodes+n;
    out[width*(size_t)i+j]=in[k];
  }
}
void psolve_biot_mass_input_pack(cudaStream_t stream,int nodes,int scalar,int count,
    const int32_t* map,const solver_float* in,solver_float* out) {
  if(count)bm_input_pack<<<(count+255)/256,256,0,stream>>>(nodes,scalar,count,map,in,out);
}
void psolve_biot_mass_dot(cudaStream_t stream,biot_mass_t m,int scalar,
    const solver_float* x,const solver_float* y) {
  cudaMemsetAsync(m.dots,0,sizeof(double),stream);
  if (m.nodes) bm_dot<<<(m.nodes+255)/256,256,0,stream>>>(m,scalar,x,y);
}
void psolve_biot_mass_dot_pair(cudaStream_t stream,biot_mass_t m,int scalar,
    const solver_float* r,const solver_float* z,const solver_float* az) {
  cudaMemsetAsync(m.dots,0,2*sizeof(double),stream);
  if (m.nodes) bm_dot_pair<<<(m.nodes+255)/256,256,0,stream>>>(m,scalar,r,z,az);
}
void psolve_biot_mass_finish_product(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,double drag_dt,const solver_float* in,const solver_float* assembled,
    solver_float* out,int precondition,const solver_float* left,int dot_mode) {
  if (dot_mode) cudaMemsetAsync(m.dots,0,dot_mode*sizeof(double),stream);
  if (!m.nodes) return;
  const int blocks=(m.nodes+255)/256;
  #define BM_FINISH(S,D) bm_finish_product<S,D><<<blocks,256,0,stream>>>(data,m,drag_dt,in,assembled,out,precondition,left)
  if (scalar) {
    if (dot_mode==2) BM_FINISH(true,2);
    else if (dot_mode==1) BM_FINISH(true,1);
    else BM_FINISH(true,0);
  } else {
    if (dot_mode==2) BM_FINISH(false,2);
    else if (dot_mode==1) BM_FINISH(false,1);
    else BM_FINISH(false,0);
  }
  #undef BM_FINISH
}
void psolve_biot_mass_local(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int scalar,int mode,double alpha,double drag_dt,const solver_float* in,
    solver_float* out,int dot) {
  if (dot) cudaMemsetAsync(m.dots,0,sizeof(double),stream);
  if (!m.nodes) return;
  const int blocks=(m.nodes+255)/256;
  #define BM_LOCAL(S,D) bm_local<S,D><<<blocks,256,0,stream>>>(data,m,mode,alpha,drag_dt,in,out)
  if (scalar) {
    if (dot) BM_LOCAL(true,1);else BM_LOCAL(true,0);
  } else {
    if (dot) BM_LOCAL(false,1);else BM_LOCAL(false,0);
  }
  #undef BM_LOCAL
}
template<int scalar>
__global__ static void bm_pack(int nodes,int field,int fields,int count,
    const biot_mass_packet_entry_t* entries,const solver_float* values) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i>=count)return;
  const auto e=entries[i];const int width=scalar?1:3;
  if(scalar==2) {
    for(int c=0;c<7;++c) {
      const size_t source=c<6?(c/3)*3*(size_t)nodes+3*(size_t)e.node+c%3:6*(size_t)nodes+e.node;
      e.packet[7*(size_t)e.entry+c]=values[source];
    }
    return;
  }
  for(int f=0;f<fields;++f)for(int d=0;d<width;++d)
    e.packet[(size_t)width*((size_t)f*e.count+e.entry)+d]=
        values[(scalar?6:3*(field+f))*(size_t)nodes+(size_t)width*e.node+d];
}
template<int scalar>
__global__ static void bm_unpack(int nodes,int field,int fields,int add,int count,
    const int32_t* offsets,const biot_mass_packet_entry_t* entries,solver_float* values) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i>=count)return;
  const int begin=offsets[i],end=offsets[i+1],node=entries[begin].node;
  const int width=scalar?1:3;
  if(scalar==2) {
    for(int c=0;c<7;++c) {
      const size_t dest=c<6?(c/3)*3*(size_t)nodes+3*(size_t)node+c%3:6*(size_t)nodes+node;
      solver_float value=values[dest];
      for(int j=add?begin:end-1;j<end;++j) {
        const auto e=entries[j];const solver_float v=e.packet[7*(size_t)e.entry+c];
        value=add?(solver_float)(value+v):v;
      }
      values[dest]=value;
    }
    return;
  }
  for(int f=0;f<fields;++f)for(int d=0;d<width;++d) {
    const size_t dest=(scalar?6:3*(field+f))*(size_t)nodes+(size_t)width*node+d;
    solver_float value=values[dest];
    // A node owns its reduction. Round after every contribution in the same
    // peer order as separate launches; assignment uses the final entry.
    for(int j=add?begin:end-1;j<end;++j) {
      const auto e=entries[j];
      const solver_float v=e.packet[(size_t)width*((size_t)f*e.count+e.entry)+d];
      value=add?(solver_float)(value+v):v;
    }
    values[dest]=value;
  }
}
void psolve_biot_mass_pack(cudaStream_t stream,int nodes,int scalar,int field,int fields,
    int count,const biot_mass_packet_entry_t* entries,const solver_float* values) {
  if(!count)return;
  const size_t blocks=((size_t)count+255)/256;
  if(scalar==2)bm_pack<2><<<blocks,256,0,stream>>>(nodes,field,fields,count,entries,values);
  else if(scalar)bm_pack<true><<<blocks,256,0,stream>>>(nodes,field,fields,count,entries,values);
  else bm_pack<false><<<blocks,256,0,stream>>>(nodes,field,fields,count,entries,values);
}
void psolve_biot_mass_unpack(cudaStream_t stream,int nodes,int scalar,int field,int fields,
    int add,int count,const int32_t* offsets,const biot_mass_packet_entry_t* entries,solver_float* values) {
  if(!count)return;
  const size_t blocks=((size_t)count+255)/256;
  if(scalar==2)bm_unpack<2><<<blocks,256,0,stream>>>(nodes,field,fields,add,count,offsets,entries,values);
  else if(scalar)bm_unpack<true><<<blocks,256,0,stream>>>(nodes,field,fields,add,count,offsets,entries,values);
  else bm_unpack<false><<<blocks,256,0,stream>>>(nodes,field,fields,add,count,offsets,entries,values);
}

/* Predict each solved increment independently at fixed timestep; SDIRK stages
 * and the two storage right-hand sides never share histories. */
__global__ static void bm_predict(biot_mass_t m,int scalar,int history,int order) {
  const size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
  const size_t count=(scalar?1:6)*(size_t)m.nodes;
  if(i>=count)return;
  const auto* a=m.history[history][0];const auto* b=m.history[history][1];
  double value=2*(double)a[i]-b[i];
  if(order==3)value=3*(double)a[i]-3*(double)b[i]+m.history[history][2][i];
  m.x[(scalar?6*(size_t)m.nodes:0)+i]=(solver_float)value;
}
void psolve_biot_mass_predict(cudaStream_t stream,biot_mass_t m,int scalar,int history,int order) {
  const size_t count=(scalar?1:6)*(size_t)m.nodes;
  if(count)bm_predict<<<(count+255)/256,256,0,stream>>>(m,scalar,history,order);
}

__device__ static double bm_history_value(solver_float* const h[3],size_t i,int order) {
  if(order==1)return h[0][i];
  if(order==2)return 2*(double)h[0][i]-h[1][i];
  return 3*(double)h[0][i]-3*(double)h[1][i]+h[2][i];
}

/* The final verified products are linear under a fixed A and B. Reuse them
 * for the starting r=b-A*x and z=B*b-B*A*x, not for final acceptance. The
 * oldest B*A*x slot is no longer needed after this thread reads its history;
 * use it to retain B*b until cache_store. No extra scratch vector or exchange. */
template<bool scalar>
__global__ static void bm_cached_predict(gpu_data_t* data,biot_mass_t m,
    int history,int order,int depth) {
  __shared__ double sum[1][256];
  const int slot=blockIdx.x*blockDim.x+threadIdx.x;
  double rz=0;
  if(slot<m.nodes) {
    const int n=m.node?m.node[slot]:slot;
    const size_t offset=scalar?6*(size_t)m.nodes:0;
    double r[7]={},z[7]={};
    for(int c=scalar?6:0;c<(scalar?7:6);++c) {
      const size_t i=bm_index(m,n,c),q=i-offset;
      const double product=bm_history_value(m.history_product[history],q,order);
      const double preconditioned=bm_history_value(m.history_preconditioned_product[history],q,order);
      const solver_float bb=m.z[i];
      m.x[i]=(solver_float)bm_history_value(m.history[history],q,order);
      r[c]=(double)m.rhs[i]-product;
      z[c]=(double)bb-preconditioned;
      m.history_preconditioned_product[history][depth-1][q]=bb;
    }
    bm_project(data,m,n,scalar,r);
    bm_project(data,m,n,scalar,z);
    for(int c=scalar?6:0;c<(scalar?7:6);++c) {
      const size_t i=bm_index(m,n,c);
      const solver_float rr=(solver_float)r[c],zz=(solver_float)z[c];
      m.r[i]=rr;m.z[i]=zz;
      if(m.owned[n]&(scalar?2:1))rz+=(double)rr*zz;
    }
  }
  bm_fused_sum<1>(m,rz,0,sum);
}
void psolve_biot_mass_cached_predict(cudaStream_t stream,gpu_data_t* data,
    biot_mass_t m,int scalar,int history,int order,int depth) {
  cudaMemsetAsync(m.dots,0,sizeof(double),stream);
  if(!m.nodes)return;
  const int blocks=(m.nodes+255)/256;
  if(scalar)bm_cached_predict<true><<<blocks,256,0,stream>>>(data,m,history,order,depth);
  else bm_cached_predict<false><<<blocks,256,0,stream>>>(data,m,history,order,depth);
}

/* Histories have just been rotated: newest B*A*x still holds this solve's
 * saved B*b. Store A*x from fresh verification and B*A*x = B*b - B*r.
 * An accepted zero start has A*0=B*A*0=0 without an operator application. */
__global__ static void bm_cache_store(biot_mass_t m,int scalar,int history,int zero) {
  const size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
  if(i>=(scalar?1:6)*(size_t)m.nodes)return;
  const size_t q=(scalar?6*(size_t)m.nodes:0)+i;
  auto* ba=m.history_preconditioned_product[history][0];
  m.history[history][0][i]=m.x[q];
  m.history_product[history][0][i]=zero?0:m.product[q];
  ba[i]=zero?0:(solver_float)((double)ba[i]-m.z[q]);
}
void psolve_biot_mass_cache_store(cudaStream_t stream,biot_mass_t m,
    int scalar,int history,int zero) {
  const size_t count=(scalar?1:6)*(size_t)m.nodes;
  if(count)bm_cache_store<<<(count+255)/256,256,0,stream>>>(m,scalar,history,zero);
}
void psolve_biot_mass_state(cudaStream_t stream,gpu_data_t* data,biot_mass_t m,
    int op,double a) {
  if (m.nodes) bm_state<<<(m.nodes+255)/256,256,0,stream>>>(data,m,op,a);
}
#endif
