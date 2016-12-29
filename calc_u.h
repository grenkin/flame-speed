#ifndef CALC_U_H_INCLUDED
#define CALC_U_H_INCLUDED

#include "input_data.h"

enum method_t { METHOD_TRAPEZOID, METHOD_RUNGE_KUTTA };

real_t calc_u(const calc_u_Data&, const Config&, const method_t = METHOD_TRAPEZOID);
real_t calc_z(const calc_u_Data&);

#endif // CALC_U_H_INCLUDED
