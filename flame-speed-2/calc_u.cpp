#include <cmath>
#include <fstream>
#include "calc_u.h"

enum temperature_t { DOWN, UP, L_ACHIEVED };

std::string temperature_str (const temperature_t t)
{
    if (t == DOWN)
        return "DOWN";
    else
        return "UP";
}

inline real_t calc_Tb (const ModelParameters& model_parameters,
    const ModelParametersToFind& model_parameters_to_find,
    const ExperimentalData& experimental_data)
{
    const real_t z = calc_z(model_parameters);
    return model_parameters.T0
        + experimental_data.Q_div_cp
        * fmin(experimental_data.phi, 1) / (experimental_data.phi + z);
}

real_t U (const ModelParameters& model_parameters,
    const ModelParametersToFind& data,
    const ExperimentalData& experimental_data, real_t T)
{
    const real_t z = calc_z(model_parameters);
    return experimental_data.Q_div_cp * model_parameters.D * data.A * pow(T, data.n)
        * pow(experimental_data.phi / (experimental_data.phi + z)
            - (T - model_parameters.T0) / experimental_data.Q_div_cp, data.alpha)
        * pow(model_parameters.nu / (experimental_data.phi + z)
            - model_parameters.nu * (T - model_parameters.T0)
            / experimental_data.Q_div_cp, data.beta)
        * (exp(- data.E_div_R / T) - exp(- data.E_div_R / model_parameters.T0));
}

temperature_t calc_temperature (real_t u,
    const ModelParameters& model_parameters,
    const ModelParametersToFind& model_parameters_to_find,
    const ExperimentalData& experimental_data, int N_trapezoid)
{
    real_t Tb = calc_Tb(model_parameters, model_parameters_to_find,
        experimental_data);
    real_t h = (Tb - model_parameters.T0) / N_trapezoid;
    long double p = 0;
    for (int i = 0; i < N_trapezoid; ++i) {
        real_t Ti = model_parameters.T0 + h * i;
        real_t Ti1 = Ti + h;
        real_t Discr =
            u * u + 4 / h * (p * p / h + u * p
            - U(model_parameters, model_parameters_to_find, experimental_data, Ti)
            - U(model_parameters, model_parameters_to_find, experimental_data, Ti1));
        if (Discr < 0)
            return DOWN;
        p = h / 2 * (u + sqrt(Discr));
    }
    return UP;
}

real_t calc_u (const ModelParameters& model_parameters,
    const ModelParametersToFind& model_parameters_to_find,
    const ExperimentalData& experimental_data, const Config& config)
{
    real_t u = config.u_init / 2;
    temperature_t temperature_res = DOWN;
    while (u <= config.max_u && temperature_res == DOWN) {
        u *= 2;
        temperature_res = calc_temperature(u,
            model_parameters, model_parameters_to_find, experimental_data,
            config.N_trapezoid);
    }
    if (u > config.max_u)
        throw umax_achieved();
    long double left = 0;
    long double right = u;
    while (right - left > config.u_eps) {
        u = (left + right) / 2;
        temperature_res = calc_temperature(u, model_parameters,
            model_parameters_to_find, experimental_data, config.N_trapezoid);
        if (temperature_res == UP)
            right = u;
        else if (temperature_res == DOWN)
            left = u;
    }
    return u;
}

// du/dp_k
real_t calc_deriv_u (int k, const ModelParameters& model_parameters,
    const ModelParametersToFind& model_parameters_to_find,
    const ExperimentalData& experimental_data, const Config& config)
{
    ModelParametersToFind model_parameters_to_find1 = model_parameters_to_find;
    real_t u = calc_u(model_parameters, model_parameters_to_find,
        experimental_data, config);

    real_t delta = config.delta(k);
    *model_parameters_to_find1.params[k] += delta;
    real_t u1 = calc_u(model_parameters, model_parameters_to_find1,
        experimental_data, config);

    real_t du_dpk = (u1 - u) / config.delta(k);
    return du_dpk;
}
