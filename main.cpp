#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <conio.h>
#include "input_data.h"
#include "calc_u.h"

class PairStream {
public:
    std::ostream &file;

    PairStream (std::ostream &f)
        : file(f)
    {}

    template<typename T>
    PairStream& operator<< (const T &t)
    {
        std::cout << t;
        file << t;
        return *this;
    }
};

int main (void)
{
    std::cout << "Identification of reaction rate parameters in a combustion model\n\n";
    //Read input data
    InputParam input_param = read_input_param();
    Config config = read_config();
    std::ifstream fe("flame_speed.txt");
    if (!fe) {
        std::cerr << "Can not open input file: flame_speed.txt\n";
        getch();
        exit(1);
    }
    int Dnum;
    fe >> Dnum;
    std::vector<real_t> phiD(Dnum), uD(Dnum), w(Dnum), Q_div_cp(Dnum);
    for (int j = 0; j < Dnum; ++j) {
        fe >> phiD[j] >> uD[j] >> w[j];
    }
    int Bnum;
    std::ifstream fb("burnt_temperature.txt");
    if (!fb) {
        std::cerr << "Can not open input file: burnt_temperature.txt\n";
        getch();
        exit(1);
    }
    fb >> Bnum;
    std::vector<real_t> phiB(Bnum), TB(Bnum);
    for (int k = 0; k < Bnum; ++k)
        fb >> phiB[k] >> TB[k];
    real_t z_val = calc_z(input_param.data);
    for (int j = 0; j < Dnum; ++j) {
        bool ok = false;
        for (int k = 0; k < Bnum - 1; ++k) {
            if (phiB[k] <= phiD[j] && phiD[j] <= phiB[k+1]) {
                // interpolate the temperature data
                real_t Tb = TB[k] + (TB[k+1] - TB[k]) * (phiD[j] - phiB[k]) / (phiB[k+1] - phiB[k]);
                Q_div_cp[j] = (Tb - input_param.data.T0) * (phiD[j] + z_val) / fmin(phiD[j], 1);
                ok = true;
                break;
            }
        }
        if (!ok) {
            std::cerr << "Undefined adiabatic temperature: phi = " << phiD[j] << "\n";
            getch();
            exit(1);
        }
    }

    std::ofstream flog("log.txt", std::ios_base::app);
    flog.precision(10);
    std::cout.precision(10);
    PairStream pstr(flog);

    // data0 contains unchangeable data and initial guess
    calc_u_Data data0 = input_param.data;
    real_t lambda_init = input_param.lambda_init;
    // deltas are steps for numerical differentiation
    real_t delta_A = config.delta_A;
    real_t delta_E_div_R = config.delta_E_div_R;
    real_t delta_alpha = config.delta_alpha;
    real_t delta_beta = config.delta_beta;
    real_t delta_n = config.delta_n;

    // data contains data for current flame speed calculations
    calc_u_Data data = data0;
    // ..._cur are current guesses of parameters
    real_t A_cur, E_div_R_cur, alpha_cur, beta_cur, n_cur;
    // F_cur contains the value of the objective function corresponding to parameters ..._cur
    real_t F_cur;
    // u_cur contains the values of the flame speeds corresponding to parameters ..._cur
    std::vector<real_t> u_cur(Dnum), u_new(Dnum);
    real_t lambda;

    // Define arrays of pointers for convenient access to 5 parameters
    const int PARAMS = 5;
    real_t *param_cur[PARAMS] = {&A_cur, &E_div_R_cur, &alpha_cur, &beta_cur, &n_cur};
    real_t *delta[PARAMS] = {&delta_A, &delta_E_div_R, &delta_alpha, &delta_beta, &delta_n};
    real_t *data_param[PARAMS] = {&data.A, &data.E_div_R, &data.alpha, &data.beta, &data.n};
    bool *optimize_param[PARAMS] = {&config.optimize_A, &config.optimize_E_div_R,
        &config.optimize_alpha, &config.optimize_beta, &config.optimize_n};
    std::string param_name[PARAMS] = {"A", "E/R", "alpha", "beta", "n"};

    lambda = lambda_init;
    // Define the initial guess
    A_cur = data0.A;
    E_div_R_cur = data0.E_div_R;
    alpha_cur = data0.alpha;
    beta_cur = data0.beta;
    n_cur = data0.n;
    // Calculate the flame speeds u_cur and the objective function F_cur
    F_cur = 0;
    for (int j = 0; j < Dnum; ++j) {
        data.phi = phiD[j];
        data.Q_div_cp = Q_div_cp[j];
        for (int p = 0; p < PARAMS; ++p)
            *data_param[p] = *param_cur[p];
        u_cur[j] = calc_u(data, config);
        F_cur += w[j] * pow(u_cur[j] - uD[j], 2);
    }
    for (int p = 0; p < PARAMS; ++p)
        pstr << param_name[p] << " = " << *param_cur[p] << "\n";
    pstr << "F = " << F_cur << "\nlambda = " << lambda << "\n\n";

    //lambda_not_changed is the number of subsequent iterations at which lambda was not changed
    int lambda_not_changed = 0;
    for (int iter = 1; iter <= input_param.steps; ++iter) {
        pstr << "iteration " << input_param.prev_iterations + iter << "\n";
        // F_deriv and F_deriv2 are approximations of the first and the second derivatives of F
        real_t F_deriv[PARAMS], F_deriv2[PARAMS];
        for (int p = 0; p < PARAMS; ++p) {
            F_deriv[p] = 0;
            F_deriv2[p] = 0;
        }
        for (int j = 0; j < Dnum; ++j) {
            for (int p = 0; p < PARAMS; ++p) {
                if (!*optimize_param[p])
                    break;
                for (int q = 0; q < PARAMS; ++q) {
                    if (q == p)
                        *data_param[q] = *param_cur[q] + *delta[q];
                    else
                        *data_param[q] = *param_cur[q];
                }
                data.phi = phiD[j];
                data.Q_div_cp = Q_div_cp[j];
                real_t u_val = calc_u(data, config);
                // derivative with respect to p-th parameter
                real_t u_deriv = (u_val - u_cur[j]) / *delta[p];
                F_deriv[p] += w[j] * (u_cur[j] - uD[j]) * u_deriv;
                F_deriv2[p] += w[j] * pow(u_deriv, 2);
            }
        }
        for (int p = 0; p < PARAMS; ++p)
            flog << "F_" << param_name[p] << " = " << F_deriv[p] << "   ";
        flog << "\n";

        int lambda_decr = 0; // the number of subsequent decreases of lambda
        while (1) {
            // Calculate the new guess
            real_t param_new[PARAMS];
            for (int p = 0; p < PARAMS; ++p) {
                if (*optimize_param[p])
                    param_new[p] = *param_cur[p] - lambda / F_deriv2[p] * F_deriv[p];
                else
                    param_new[p] = *param_cur[p];
            }
            // Calculate the flame speeds and the objective function for the new guess
            real_t F_new = 0;
            for (int j = 0; j < Dnum; ++j) {
                data.phi = phiD[j];
                data.Q_div_cp = Q_div_cp[j];
                for (int p = 0; p < PARAMS; ++p)
                    *data_param[p] = param_new[p];
                u_new[j] = calc_u(data, config);
                F_new += w[j] * pow(u_new[j] - uD[j], 2);
            }
            if (F_new <= F_cur) {
                // Accept the new guess
                for (int p = 0; p < PARAMS; ++p)
                    *param_cur[p] = param_new[p];
                F_cur = F_new;
                for (int j = 0; j < Dnum; ++j)
                    u_cur[j] = u_new[j];
                break;
            }
            else {
                // Reject the new guess and decrease lambda
                lambda_not_changed = -1;
                lambda /= 2;
                pstr << "lambda = " << lambda << "\n\n";
                ++lambda_decr;
                if (lambda_decr == config.lambda_decr_max) {
                    pstr << "lambda_decr_max achieved\n";
                    break;
                }
            }
        }
        if (lambda_decr == config.lambda_decr_max)
            break;
        for (int p = 0; p < PARAMS; ++p)
            pstr << param_name[p] << " = " << *param_cur[p] << "\n";
        pstr << "F = " << F_cur << "\n\n";
        ++lambda_not_changed;
        // If lambda was not decreased for many times, then increase lambda;
        // if lambda was not increased for many times, then increase lambda again
        if (lambda_not_changed == config.lambda_threshold) {
            lambda_not_changed = 0;
            lambda *= 2;
            pstr << "lambda increased = " << lambda << "\n\n";
        }

        // Print output parameters
        std::ofstream fout("input_param.txt");
        fout.precision(10);
        fout << "T0 = " << data0.T0 << "\n";
        fout << "D = " << data0.D << "\n";
        fout << "nu = " << data0.nu << "\n\n";
        for (int p = 0; p < PARAMS; ++p)
            fout << param_name[p] << " = " << *param_cur[p] << "\n";
        fout << "\nlambda = " << lambda << "\n";
        fout << "prev_iterations = " << input_param.prev_iterations + iter << "\n\n";
        fout << "steps = " << input_param.steps;

    } // for iter

    // Print data for plotting
    std::ofstream fout_u("plot.txt");
    for (int i = 0; i < Dnum; ++i)
        fout_u << phiD[i] << "   " << u_cur[i] << "   " << uD[i] << "\n";

    // Print output parameters
    std::ofstream fout("output_param.txt");
    fout.precision(10);
    fout << "T0 = " << data0.T0 << "\n";
    fout << "D = " << data0.D << "\n";
    fout << "nu = " << data0.nu << "\n\n";
    for (int p = 0; p < PARAMS; ++p)
        fout << param_name[p] << " = " << *param_cur[p] << "\n";
    fout << "\nlambda = " << lambda << "\n";
    fout << "prev_iterations = " << input_param.prev_iterations + input_param.steps << "\n\n";
    fout << "steps = " << input_param.steps;

    return 0;
}
