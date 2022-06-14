#ifndef CALC_U_H_INCLUDED
#define CALC_U_H_INCLUDED

#include "input_data.h"

real_t calc_u(const ModelParameters&, const ModelParametersToFind&,
    const ExperimentalData&, const Config&);
real_t calc_deriv_u(int k, const ModelParameters&,
    const ModelParametersToFind&, const ExperimentalData&, const Config&);

class umax_achieved {};

#endif // CALC_U_H_INCLUDED
