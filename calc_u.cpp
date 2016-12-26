#include <cmath>
#include <fstream>
#include "calc_u.h"

enum temperature_t { DOWN, UP };
class umax_achieved {};

std::string temperature_str (const temperature_t t)
{
    if (t == DOWN)
        return "DOWN";
    else
        return "UP";
}

real_t U(const calc_u_Data &data, const real_t T)
{
    const real_t z = data.nu * (1 + 3.762 * 28. / 32.);
    return data.Q_div_cp * data.D * data.A * pow(T, data.n)
        * pow(data.phi / (data.phi + z) - (T - data.T0) / data.Q_div_cp, data.alpha)
        * pow(data.nu / (data.phi + z) - data.nu * (T - data.T0) / data.Q_div_cp, data.beta)
        * (exp(- data.E_div_R / T) - exp(- data.E_div_R / data.T0));
}

temperature_t calc_temperature(const real_t u, const calc_u_Data &data, const int N_trapezoid)
{
    const real_t z = data.nu * (1 + 3.762 * 28. / 32.);
    real_t Tb = data.T0 + data.Q_div_cp * fmin(data.phi, 1) / (data.phi + z);
    real_t h = (Tb - data.T0) / N_trapezoid;
    long double p = 0;
    for (int i = 0; i < N_trapezoid; ++i) {
        real_t Ti = data.T0 + h * i;
        real_t Ti1 = Ti + h;
        real_t D = u * u + 4 / h * (p * p / h + u * p - U(data, Ti) - U(data, Ti1));
        if (D < 0)
            return DOWN;
        p = h / 2 * (u + sqrt(D));
    }
    return UP;
}

real_t calc_u (const calc_u_Data &data, const Config &config)
{
    real_t u = config.u_init / 2;
    temperature_t temperature_res = DOWN;
    while (u <= config.max_u && temperature_res == DOWN) {
        u *= 2;
        temperature_res = calc_temperature(u, data, config.N_trapezoid);
    }
    if (u > config.max_u)
        throw umax_achieved();
    long double left = 0;
    long double right = u;
    while (right - left > config.u_eps) {
        u = (left + right) / 2;
        temperature_res = calc_temperature(u, data, config.N_trapezoid);
        if (temperature_res == UP)
            right = u;
        else if (temperature_res == DOWN)
            left = u;
    }
    return u;
}
