#include <cmath>
#include <fstream>
#include "calc_u.h"

#ifdef RUNGE_KUTTA
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#endif

enum temperature_t { DOWN, UP, L_ACHIEVED };
class umax_achieved {};

std::string temperature_str (const temperature_t t)
{
    if (t == DOWN)
        return "DOWN";
    else
        return "UP";
}

real_t calc_z(const calc_u_Data &data)
{
    return data.nu * (1 + 3.762 * 28. / 32.);
}

inline real_t calc_Tb(const calc_u_Data& data)
{
    const real_t z = calc_z(data);
    return data.T0 + data.Q_div_cp * fmin(data.phi, 1) / (data.phi + z);
}

real_t U(const calc_u_Data &data, const real_t T)
{
    const real_t z = calc_z(data);
    return data.Q_div_cp * data.D * data.A * pow(T, data.n)
        * pow(data.phi / (data.phi + z) - (T - data.T0) / data.Q_div_cp, data.alpha)
        * pow(data.nu / (data.phi + z) - data.nu * (T - data.T0) / data.Q_div_cp, data.beta)
        * (exp(- data.E_div_R / T) - exp(- data.E_div_R / data.T0));
}

temperature_t calc_temperature(const real_t u, const calc_u_Data &data, const int N_trapezoid)
{
    real_t Tb = calc_Tb(data);
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

struct rk_Data {
    calc_u_Data data;
    real_t u;
};

#ifdef RUNGE_KUTTA
int func(double t, const double y[], double dydt[], void *params)
{
    rk_Data* rk_data = (rk_Data*)params;
    real_t Tb = calc_Tb(rk_data->data);
    if (y[0] < rk_data->data.T0 || y[0] > Tb) {
        dydt[0] = y[1];
        dydt[1] = 0;
    }
    else {
        dydt[0] = y[1];
        dydt[1] = rk_data->u * y[1] - U(rk_data->data, y[0]);
    }
    return GSL_SUCCESS;
}

temperature_t calc_temperature_rk (const real_t u, const calc_u_Data& data)
{
    const double L = 10000;
    const double eps_ic = 0.001;
    const double eps_abs = 1e-4;

    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 2);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(eps_abs, 0.0);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(2);

    rk_Data rk_data = { data, u };
    gsl_odeiv2_system sys = {func, NULL, 2, &rk_data};

    double t = -L, t1 = L;
    double h = 1e-6;
    double maxh = 1e-1;
    double y[2] = { data.T0, eps_ic };

    real_t Tb = calc_Tb(data);
    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            printf ("error, return value=%d\n", status);
            break;
        }
        if (y[0] > Tb) {
            return UP;
        }
        if (y[1] < 0) {
            return DOWN;
        }
        if (h > maxh)
            h = maxh;
    }
    throw "L achieved";
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    return L_ACHIEVED;
}
#endif

real_t calc_u (const calc_u_Data &data, const Config &config, const method_t method)
{
    real_t u = config.u_init / 2;
    temperature_t temperature_res = DOWN;
    while (u <= config.max_u && temperature_res == DOWN) {
        u *= 2;
#ifdef RUNGE_KUTTA
        if (method == METHOD_TRAPEZOID)
#endif
            temperature_res = calc_temperature(u, data, config.N_trapezoid);
#ifdef RUNGE_KUTTA
        else
            temperature_res = calc_temperature_rk(u, data);
#endif
    }
    if (u > config.max_u)
        throw umax_achieved();
    long double left = 0;
    long double right = u;
    while (right - left > config.u_eps) {
        u = (left + right) / 2;
#ifdef RUNGE_KUTTA
        if (method == METHOD_TRAPEZOID)
#endif
            temperature_res = calc_temperature(u, data, config.N_trapezoid);
#ifdef RUNGE_KUTTA
        else
            temperature_res = calc_temperature_rk(u, data);
#endif
        if (temperature_res == UP)
            right = u;
        else if (temperature_res == DOWN)
            left = u;
    }
    return u;
}
