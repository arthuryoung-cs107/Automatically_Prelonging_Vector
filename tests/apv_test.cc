#include "apv_ode_systems.hh"
#include <cstdio>

ode_solspc_meta meta(2,1);

Duffing_ode ode;
// VanDerPol_ode ode;
// Pendulum_ode ode;
// Bessel_ode ode;
// Riccati_ode ode;
// Brusselator_ode ode;

ode_solspc_meta meta0(ode.eor,ode.ndep);

int main()
{

  printf("YER.\n");

  return 0;
}
