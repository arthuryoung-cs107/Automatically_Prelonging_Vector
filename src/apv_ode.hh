#ifndef APV_ODE_HH
#define APV_ODE_HH

#include "apv_util.hh"

struct ode_solspc_meta
{
  ode_solspc_meta(int eor_,int ndep_): eor(eor_), ndep(ndep_) {}
  ~ode_solspc_meta() {}

  const int eor,
            ndep,
            ndim = 1 + ndep*(eor+1),
            nvar = 1 + ndep;
};

struct ode_solspc_element
{
  ode_solspc_element(ode_solspc_meta &meta_): meta(meta_) {}
  ~ode_solspc_element() {}

  ode_solspc_meta &meta;

  const int &eor = meta.eor,
            &ndep = meta.ndep,
            &ndim = meta.ndim,
            &nvar = meta.nvar;
};

struct ode_solution: public ode_solspc_element
{
  ode_solution(ode_solspc_meta &meta_, double *pts_): ode_solspc_element(meta_), pts(pts_) {}
  ~ode_solution() {}

  double  * const pts,
          &x = pts[0],
          * const u = pts + 1,
          * const dxu = u + ndep,
          * const dnxu = u + ndep*eor;

  void print_sol();
};

struct ode_system: public ode_solspc_meta
{
  ode_system(int eor_, int ndep_): ode_solspc_meta(eor_,ndep_) {}
  ~ode_system() {}

  virtual void dudx_eval(double x_, double *u_, double *dudx_) = 0;
  virtual void JacF_eval(double x_, double *u_, double **dls_out_) = 0;
  virtual void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) = 0;

};


#endif
