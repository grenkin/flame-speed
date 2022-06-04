#ifndef CALC_U_H_INCLUDED
#define CALC_U_H_INCLUDED

#include "input_data.h"

real_t calc_u(const calc_u_Data&, const Config&);
real_t calc_z(const calc_u_Data&);

class umax_achieved {};

#endif // CALC_U_H_INCLUDED
