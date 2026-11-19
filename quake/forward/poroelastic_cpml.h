/* Mixed displacement/pressure CFS-CPML. Element-constant Cartesian stretches.
 * With J=Sx Sy Sz the weak operators are J R dtt, J B dt,
 * K(J/(Sj Sl)), Q(J/Sj), and J/M storage. The elastic kernel supplies
 * the drained K and solid diagonal inertia. This file supplies the mixed
 * terms without converting a porous layer into an elastic exterior.
 */
#ifndef HERCULES_POROELASTIC_CPML_H
#define HERCULES_POROELASTIC_CPML_H
#include "cpml.h"
#ifdef __CUDACC__
#define PORO_CPML_HD __host__ __device__
#else
#define PORO_CPML_HD
#endif

typedef struct poro_cpml_t {
  cpml_element_t stretch;
  cpml_inertia_t inertia, darcy; /* J*s^2 and J*s (leading drag handled locally) */
  double h, volume, storage, dt, rho, rho_f, fluid_mass, drag[3];
  int consistent; /* exact Q1 inertia/storage with matching history products */
  double drain[8]; /* Unstretched BrickUP drainage compliance; zero for full Biot. */
} poro_cpml_t;

/* 48 traction + 48 continuity + 24 storage + 36 reduced diffusion,
 * followed by 72 solid + 72 relative-fluid + 72 drag histories for full Biot.
 * All arrays are compact PML-only, SoA by channel. */
enum { PORO_CPML_CONTINUITY=48, PORO_CPML_STORAGE=96,
       PORO_CPML_DIFFUSION=120, PORO_CPML_INERTIA=156,
       PORO_CPML_REDUCED_STATES=156, PORO_CPML_BIOT_STATES=372 };

#ifdef __cplusplus
static inline cpml_inertia_t poro_cpml_darcy(const double d[3],
    const double alpha[3], const double kappa[3], double dt) {
  const cpml_inertia_t a=cpml_inertia(d,alpha,kappa,3,0,dt);
  cpml_inertia_t b=cpml_inertia(d,alpha,kappa,3,1,dt);
  b.mass-=a.mass; b.damping-=a.damping; b.stiffness-=a.stiffness;
  for (int i=0;i<3;++i) b.residue[i]-=a.residue[i];
  return b;
}

PORO_CPML_HD static inline double poro_cpml_capacity(const poro_cpml_t& c) {
  return c.stretch.forward[0].direct*c.stretch.forward[1].direct*
         c.stretch.forward[2].direct;
}

/* Ni*Nj integral divided by the row-lumped volume h^3/8. Keeping the
 * element-constant stretch outside this matrix pairs leading and memory
 * operators exactly. This separable three-axis application costs 24 pairs,
 * rather than 64 dense products, and preserves constants. */
template<int components>
PORO_CPML_HD static inline void poro_cpml_q1_average(double (&v)[8][components]) {
  for (int d=0;d<3;++d) for (int i=0;i<8;++i) if (!(i&(1<<d))) {
    const int j=i|(1<<d);
    for (int c=0;c<components;++c) {
      const double a=v[i][c],b=v[j][c];
      v[i][c]=(2*a+b)/3; v[j][c]=(a+2*b)/3;
    }
  }
}

/* Exact cube integral int dNi/dx_d Nj dV. */
PORO_CPML_HD static inline double poro_cpml_q(int i,int j,int d,double h) {
  const int a=(d+1)%3,b=(d+2)%3;
  return (((i>>d)&1)?1:-1)*h*h/72*
      ((((i^j)>>a)&1)?1:2)*((((i^j)>>b)&1)?1:2);
}

template<class T>
PORO_CPML_HD static inline double poro_cpml_cofactor(const poro_cpml_t& c,
    int d,double v,T* state,size_t stride) {
  v=cpml_apply(c.stretch.forward[(d+1)%3],v,state);
  return cpml_apply(c.stretch.forward[(d+2)%3],v,state+stride);
}

template<class T>
PORO_CPML_HD static inline void poro_cpml_traction(const poro_cpml_t& c,
    const double p[8],double force[8][3],T* state,size_t stride) {
  for (int i=0;i<8;++i) for (int d=0;d<3;++d) {
    double f=0;
    for (int j=0;j<8;++j) f+=poro_cpml_q(i,j,d,c.h)*p[j];
    force[i][d]=poro_cpml_cofactor(c,d,f,state+2*(3*i+d)*stride,stride);
  }
}

/* Solve increments of C J p = -Q(J/S)^T (alpha*u+w). A cascade has
 * y[n]=D*p[n]+tail[n]. Update its state with p[n], then subtract
 * C*(tail[n+1]-tail[n]); assembled C*D is the pressure update denominator.
 * This retains a nodal pressure shared by cells with different profiles.
 * BrickUP's per-update drainage lift adds D*J*p[n+1]. Its direct part is
 * assembled with the surface capacity; its tail uses this same J cascade.
 * No additional history, nodal exchange or implicit global solve is needed. */
template<bool reducedDrain=false,class T>
PORO_CPML_HD static inline void poro_cpml_continuity(const poro_cpml_t& c,
    const double delta[8][3],const double p[8],double rhs[8],T* state,size_t stride) {
  double storage_tail[8][1];
  for (int j=0;j<8;++j) {
    rhs[j]=0;
    for (int d=0;d<3;++d) {
      double v=0;
      for (int i=0;i<8;++i) v+=poro_cpml_q(i,j,d,c.h)*delta[i][d];
      rhs[j]-=poro_cpml_cofactor(c,d,v,
          state+(PORO_CPML_CONTINUITY+2*(3*j+d))*stride,stride);
    }
    T* s=state+(PORO_CPML_STORAGE+3*j)*stride;
    double old_tail=0,new_tail=0,x=p[j];
    for (int d=0;d<3;++d) {
      if (c.stretch.forward[d].feed==0.0) {
        old_tail=c.stretch.forward[d].direct*old_tail;
        x=c.stretch.forward[d].direct*x;
        new_tail=c.stretch.forward[d].direct*new_tail;
        continue;
      }
      old_tail=c.stretch.forward[d].direct*old_tail+(double)s[d*stride];
      x=cpml_apply(c.stretch.forward[d],x,s+d*stride);
      new_tail=c.stretch.forward[d].direct*new_tail+(double)s[d*stride];
    }
    if (c.consistent) storage_tail[j][0]=c.volume/c.storage*(new_tail-old_tail);
    else rhs[j]-=c.volume/c.storage*(new_tail-old_tail);
    if (reducedDrain) rhs[j]-=c.drain[j]*new_tail;
  }
  if (c.consistent) {
    poro_cpml_q1_average(storage_tail);
    for (int j=0;j<8;++j) rhs[j]-=storage_tail[j][0];
  }
}

template<class T>
PORO_CPML_HD static inline void poro_cpml_diffusion(const poro_cpml_t& c,
    const double mobility[3],const double p[8],double rhs[8],T* state,size_t stride) {
  double modes[8]={},fm[8]={};
  for (int m=1;m<8;++m) for (int j=0;j<8;++j)
    modes[m]+=0.125*cpml_sign(j,m)*p[j];
  for (int d=0;d<3;++d) for (int t=0;t<4;++t) {
    const int a=(d+1)%3,b=(d+2)%3;
    const int m=(1<<d)|((t&1)?1<<a:0)|((t&2)?1<<b:0);
    const double weight=0.5/((t&1)?3:1)/((t&2)?3:1);
    fm[m]+=weight*c.dt*c.h*mobility[d]*cpml_diagonal(c.stretch,d,modes[m],
        state+(PORO_CPML_DIFFUSION+3*(4*d+t))*stride,stride);
  }
  for (int j=0;j<8;++j) for (int m=1;m<8;++m) rhs[j]-=cpml_sign(j,m)*fm[m];
}

template<class T>
PORO_CPML_HD static inline void poro_cpml_fluid_inertia(const poro_cpml_t& c,
    const double u[8][3],const double w[8][3],const double du[8][3],
    const double dw[8][3],double fu[8][3],double fw[8][3],T* state,size_t stride) {
  for (int j=0;j<8;++j) for (int d=0;d<3;++d) {
    const int channel=3*(3*j+d);
    const double hu=cpml_inertia_force(c.inertia,u[j][d],
        state+(PORO_CPML_INERTIA+channel)*stride,stride);
    const double hw=cpml_inertia_force(c.inertia,w[j][d],
        state+(PORO_CPML_INERTIA+72+channel)*stride,stride);
    const double drag=c.drag[d]==0.0?0.0:cpml_inertia_force(c.darcy,w[j][d],
        state+(PORO_CPML_INERTIA+144+channel)*stride,stride);
    const double U=c.dt*c.dt*hu-c.dt*c.inertia.damping*du[j][d];
    const double W=c.dt*c.dt*hw-c.dt*c.inertia.damping*dw[j][d];
    fu[j][d]+=c.volume*c.rho_f*W;
    fw[j][d]+=c.volume*(c.rho_f*U+c.fluid_mass*W+c.dt*c.dt*c.drag[d]*drag);
  }
}

/* Consistent full-Biot CPML: spatially integrate ALL polynomial/restoring/
 * memory parts of J*R*s^2 with the same Q1 matrix as leading inertia. Darcy
 * drag and physical solid mass damping remain row-lumped, as in the interior;
 * their complete J histories retain that same lumping. Base and material
 * solid inertia share poles, so one u history suffices for both operators.
 * The elastic CPML kernel supplies stiffness only on these elements. */
template<class T>
PORO_CPML_HD static inline void poro_cpml_consistent_inertia(const poro_cpml_t& c,
    const double u[8][3],const double w[8][3],const double du[8][3],
    const double dw[8][3],double fu[8][3],double fw[8][3],T* state,size_t stride) {
  double coupled[8][6],solid_damping[8][3],fluid_drag[8][3];
  for (int j=0;j<8;++j) for (int d=0;d<3;++d) {
    const int channel=3*(3*j+d);
    T* us=state+(PORO_CPML_INERTIA+channel)*stride;
    double x=u[j][d],hu=-c.inertia.stiffness*x;
    double ha=-(c.stretch.inertia.stiffness-c.inertia.stiffness)*x;
    for (int k=0;k<c.inertia.order;++k) {
      x=cpml_apply(c.inertia.pole[k],x,us+k*stride);
      hu-=c.inertia.residue[k]*x;
      ha-=(c.stretch.inertia.residue[k]-c.inertia.residue[k])*x;
    }
    const double hw=cpml_inertia_force(c.inertia,w[j][d],
        state+(PORO_CPML_INERTIA+72+channel)*stride,stride);
    const double drag=c.drag[d]==0?0:cpml_inertia_force(c.darcy,w[j][d],
        state+(PORO_CPML_INERTIA+144+channel)*stride,stride);
    const double U=c.dt*c.dt*hu-c.dt*c.inertia.damping*du[j][d];
    const double W=c.dt*c.dt*hw-c.dt*c.inertia.damping*dw[j][d];
    coupled[j][d]=c.volume*(c.rho*U+c.rho_f*W);
    coupled[j][d+3]=c.volume*(c.rho_f*U+c.fluid_mass*W);
    solid_damping[j][d]=c.volume*c.rho*c.dt*c.dt*ha;
    fluid_drag[j][d]=c.volume*c.dt*c.dt*c.drag[d]*drag;
  }
  poro_cpml_q1_average(coupled);
  for (int j=0;j<8;++j) for (int d=0;d<3;++d) {
    fu[j][d]+=coupled[j][d]+solid_damping[j][d];
    fw[j][d]+=coupled[j][d+3]+fluid_drag[j][d];
  }
}
#endif /* __cplusplus */
#undef PORO_CPML_HD
#endif
