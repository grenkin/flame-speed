#ifndef INPUT_DATA_H_INCLUDED
#define INPUT_DATA_H_INCLUDED

#include <fstream>
#include <iostream>
#include <vector>
#include <conio.h>

typedef long double real_t;

struct ModelParameters {
    real_t T0, D, nu;
};

const int PARAMS_NUM = 5;
const std::string PARAMS_NAMES[5] = {
    "A", "E/R", "alpha", "beta", "n"};

struct ModelParametersToFind {
    real_t A, E_div_R, alpha, beta, n;

    real_t& param(int p)
    {
        return *(&A + p);
    }
};

struct BurntTemperature {
    int burnt_size;  // size of vectors
    std::vector<real_t> phi_burnt, temperature_burnt;

    BurntTemperature (const char *file_name)
    {
        std::ifstream fburnt(file_name);
        if (!fburnt) {
            std::cerr << "Can not open input file: " << file_name << "\n";
            getch();
            exit(1);
        }

        // The input file contains a sequence of pairs (phi, Tb(phi)).
        fburnt >> burnt_size;
        phi_burnt.resize(burnt_size);
        temperature_burnt.resize(burnt_size);
        for (int k = 0; k < burnt_size; ++k)
            fburnt >> phi_burnt[k] >> temperature_burnt[k];
    }

    real_t get_Tb (real_t phi) const
    {
        for (int k = 0; k < burnt_size - 1; ++k) {
            if (phi_burnt[k] <= phi && phi <= phi_burnt[k+1]) {
                // interpolate the temperature data
                real_t Tb = temperature_burnt[k]
                    + (temperature_burnt[k+1] - temperature_burnt[k])
                    * (phi - phi_burnt[k]) / (phi_burnt[k+1] - phi_burnt[k]);
                return Tb;
            }
        }

        std::cerr << "Undefined adiabatic temperature for phi = "
            << phi << "\n";
        getch();
        exit(1);
    }
};

real_t calc_z (const ModelParameters& model_parameters);

struct ExperimentalData {
    real_t phi, Q_div_cp, v;

    ExperimentalData ()
    {}

    ExperimentalData (real_t phi, real_t v,
        const BurntTemperature& burnt_temperature,
        ModelParameters model_parameters)
        : phi(phi), v(v)
    {
        real_t Tb = burnt_temperature.get_Tb(phi);
        real_t z_val = calc_z(model_parameters);
        Q_div_cp = (Tb - model_parameters.T0) * (phi + z_val) / fmin(phi, 1);
    }
};

struct ParamRange {
    real_t left, right;
    int num;  // number of intervals that [left, right] is separated into
};

struct InputParam {
    ModelParameters model_parameters;
    std::vector<ParamRange> params_ranges;
    real_t F_min;  // estimate of minimum from below of the objective function

    InputParam ()
        : params_ranges(PARAMS_NUM)
    {}
};

struct Config {
    // steps for numerical differentiation
    real_t delta_A, delta_E_div_R, delta_alpha, delta_beta, delta_n;
    int N_trapezoid;
    real_t u_eps, max_u, u_init, max_speed;

    real_t delta (int param) const
    {
        return *(&delta_A + param);
    }
};

InputParam read_input_param ();
Config read_config ();

const std::string file_input_param = "input_param.txt";
const std::string file_config = "config.txt";

#endif // INPUT_DATA_H_INCLUDED
