#ifndef APV_ODE_SYS_HH
#define APV_ODE_SYS_HH

#include "apv_ode.hh"
#include <cstring>
#include <cmath>

struct known_ode: public ode_system
{
  known_ode(int eor_, int ndep_, const char name_[]): ode_system(eor_,ndep_),
  name(new char[strlen(name_)+1])
  {strcpy(name,name_);}
  ~known_ode() {delete [] name;}

  char * const name;
};
  class Riccati_ode: public known_ode
  {
    public:
      Riccati_ode(double coeff1_=2.0, double coeff2_=-1.0):
        known_ode(1,1,"Riccati"),
        coeff1(coeff1_), coeff2(coeff2_) {}
      ~Riccati_ode() {}

      void dudx_eval(double x_, double *u_, double *dudx_)
        {dudx_[0] = (coeff1*(u_[0]/x_)) + (coeff2*(x_*x_)*(u_[0]*u_[0]));}
      void JacF_eval(double x_, double *u_, double **Jac_out_)
      {
        double x_sqrd = x_*x_;
        Jac_out_[0][0] = ((coeff1*u_[0])/x_sqrd) - ((2.0*coeff2)*x_*(u_[0]*u_[0]));
        Jac_out_[0][1] = ((-1.0*coeff1)/x_) - ((2.0*coeff2)*(x_sqrd)*u_[0]);
        Jac_out_[0][2] = 1.0;
      }
      void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
        {dnp1xu_[0] = (coeff1*((u_[1]*x_)-u_[0])/(x_*x_)) + (coeff2*(2.0*x_*x_*u_[0]*u_[0])*((u_[1]/u_[0])+(1.0/x_)));}

      inline const double * get_default_IC_indep_range(int xrange_=0) {return Riccati_x_range_def[xrange_];}
      inline const double * get_default_IC_range(int icrange_=0) {return Riccati_IC_range_def[icrange_];}

    private:

      const int dparam_len = 2;
      const double  coeff1,
                    coeff2,
                    * const dparams = &coeff1;

      const double Riccati_x_range_def[1][2] =
      {
        {0.2, 2.0}
      };
      const double Riccati_IC_range_def[1][2] =
      {
        // u
        {0.02, 1.0}
      };
  };

  class Pendulum_ode: public known_ode
  {
    public:
      Pendulum_ode(double r_=2.0, double mass_=1.0, double gamma_=1.0, double Aforce_=0.0, double Fforce_=1.0, double gravity_=-9.81):
        known_ode(2,1,"Pendulum"),
        radius(r_), mass(mass_), gamma(gamma_), Aforce(Aforce_), Fforce(Fforce_), gravity(gravity_) {}
      ~Pendulum_ode() {}

      void dudx_eval(double x_, double *u_, double *dudx_)
      {
        const double  u_in = u_[0],
                      dudx_in = u_[1];
        dudx_[0] = dudx_in;
        dudx_[1] = (c1*sin(u_in)) - (c2*dudx_in) + c3*(cos(Fforce*x_));
      }
      void JacF_eval(double x_, double *u_, double **Jac_out_)
      {
        Jac_out_[0][0] = Fforce*c3*sin(Fforce*x_);
        Jac_out_[0][1] = -1.0*c1*cos(u_[0]);
        Jac_out_[0][2] = c2;
        Jac_out_[0][3] = 1.0;
      }
      void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
        {dnp1xu_[0] = (c1*cos(u_[0]))*u_[1] - (c2*u_[2]) - c3*(sin(Fforce*x_)*(Fforce));}

      inline const double * get_default_IC_indep_range(int xrange_=0) {return Pendulum_x_range_def[xrange_];}
      inline const double * get_default_IC_range(int icrange_=0) {return Pendulum_IC_range_def[icrange_][0];}

    private:

      const int dparam_len = 6;
      const double  radius,
                    mass,
                    gamma,
                    Aforce,
                    Fforce,
                    gravity,
                    * const dparams = &radius;

      const double  pi_private = 3.14159265358979323846,
                    c1 = gravity/radius,
                    c2 = gamma/(mass*radius*radius),
                    c3 = Aforce/(mass*radius*radius);

      const double Pendulum_x_range_def[1][2] =
      {
        {0.0, 4.0*pi_private}
      };
      const double Pendulum_IC_range_def[1][2][2] =
      {
        // u,              dxu
        { {-1.25, 1.25}, {-1.0, 1.0} }
      };
  };
  class Duffing_ode: public known_ode
  {
    public:
      // default parameters inducing chaotic trajectories
      Duffing_ode(double alpha_=-1.0, double beta_=1.0, double gamma_=0.5, double delta_=0.3, double omega_=1.2):
        known_ode(2,1,"Duffing"),
        alpha(alpha_), beta(beta_), gamma(gamma_), delta(delta_), omega(omega_) {}
      ~Duffing_ode() {}

      void dudx_eval(double x_, double *u_, double *dudx_)
      {
        const double  u_in = u_[0],
                      dudx_in = u_[1];
        dudx_[0] = dudx_in;
        dudx_[1] = (gamma*cos(omega*x_))-((delta*dudx_in) + u_in*(alpha + (beta*(u_in*u_in))));
      }
      void JacF_eval(double x_, double *u_, double **Jac_out_)
      {
        Jac_out_[0][0] = gamma*omega*sin(omega*x_);
        Jac_out_[0][1] = alpha + (3.0*beta*u_[0]*u_[0]);
        Jac_out_[0][2] = delta;
        Jac_out_[0][3] = 1.0;
      }
      void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
        {dnp1xu_[0] = -1.0*(delta*u_[2] + u_[1]*(alpha + 3.0*beta*u_[0]*u_[0]) + gamma*omega*sin(omega*x_));}

      inline const double * get_default_IC_indep_range(int xrange_=0) {return Duffing_x_range_def[xrange_];}
      inline const double * get_default_IC_range(int icrange_=0) {return Duffing_IC_range_def[icrange_][0];}

    private:

      const int dparam_len = 5;
      const double  alpha,
                    beta,
                    gamma,
                    delta,
                    omega,
                    * const dparams = &alpha;

      const double  pi_private = 3.14159265358979323846;
      const double Duffing_x_range_def[5][2] =
      {
        {0.0, 5.0*(2.0*pi_private)/omega},
        {0.0, 4.0*(2.0*pi_private)/omega},
        {0.0, 6.0*(2.0*pi_private)/omega},
        {0.0, 7.0*(2.0*pi_private)/omega},
        {0.0, 7.5*(2.0*pi_private)/omega}
      };
      const double Duffing_IC_range_def[1][2][2] =
      {
        // u,              dxu
        { {-1.25, 1.25}, {-1.25, 1.25} }
      };
  };


#endif
