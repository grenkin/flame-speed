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
    PairStream(std::ostream &f): file(f) {}

    template<typename T>
    PairStream& operator<<(const T &t)
    {
        std::cout << t;
        file << t;
        return *this;
    }

    std::ostream &file;
};

int main (void)
{
    std::cout << "Identification of reaction rate parameters in a combustion model\n\n";
    InputParam input_param = read_input_param();
    Config config = read_config();
    std::ifstream fin("experiment.txt");
    if (!fin) {
        std::cerr << "Can not open input file: experiment.txt" << "\n";
        getch();
        exit(1);
    }
    int Dnum;
    fin >> Dnum;
    std::vector<real_t> phiD(Dnum), uD(Dnum), w(Dnum);
    for (int j = 0; j < Dnum; ++j) {
        fin >> phiD[j] >> uD[j] >> w[j];
    }

    std::ofstream flog("log.txt");
    flog.precision(10);
    std::cout.precision(10);
    PairStream pstr(flog);

    calc_u_Data data0 = input_param.data;
    real_t lambda_init = input_param.lambda_init;
    real_t delta_A = config.delta_A;
    real_t delta_E_div_R = config.delta_E_div_R;
    real_t delta_alpha = config.delta_alpha;
    real_t delta_beta = config.delta_beta;
    real_t delta_n = config.delta_n;

    calc_u_Data data = data0;
    real_t A_cur, E_div_R_cur, alpha_cur, beta_cur, n_cur;
    real_t F;
    std::vector<real_t> u_cur(Dnum), u_new(Dnum);
    real_t lambda;

    const int PARAMS = 5;
    real_t *param_cur[PARAMS] = {&A_cur, &E_div_R_cur, &alpha_cur, &beta_cur, &n_cur};
    real_t *delta[PARAMS] = {&delta_A, &delta_E_div_R, &delta_alpha, &delta_beta, &delta_n};
    real_t *data_param[PARAMS] = {&data.A, &data.E_div_R, &data.alpha, &data.beta, &data.n};
    bool *optimize_param[PARAMS] = {&config.optimize_A, &config.optimize_E_div_R,
        &config.optimize_alpha, &config.optimize_beta, &config.optimize_n};
    std::string param_name[PARAMS] = {"A", "E/R", "alpha", "beta", "n"};

    lambda = lambda_init;
    A_cur = data0.A;
    E_div_R_cur = data0.E_div_R;
    alpha_cur = data0.alpha;
    beta_cur = data0.beta;
    n_cur = data0.n;
    F = 0;
    for (int j = 0; j < Dnum; ++j) {
        data.phi = phiD[j];
        for (int p = 0; p < PARAMS; ++p)
            *data_param[p] = *param_cur[p];
        u_cur[j] = calc_u(data, config);
        F += w[j] * pow(u_cur[j] - uD[j], 2);
    }
    for (int p = 0; p < PARAMS; ++p)
        pstr << param_name[p] << " = " << *param_cur[p] << "\n";
    pstr << "F = " << F << "\nlambda = " << lambda << "\n\n";

    //the number of subsequent iterations at which lambda was not changed
    int lambda_not_changed = 0;
    for (int iter = 1; iter <= input_param.steps; ++iter) {
        pstr << "iteration " << input_param.prev_iterations + iter << "\n";
        real_t F_deriv[PARAMS], F_deriv2[PARAMS];
        for (int p = 0; p < PARAMS; ++p) {
            F_deriv[p] = 0;
            F_deriv2[p] = 0;
        }
        for (int j = 0; j < Dnum; ++j) {
            for (int p = 0; p < PARAMS; ++p) {
                for (int q = 0; q < PARAMS; ++q) {
                    if (q == p)
                        *data_param[q] = *param_cur[q] + *delta[q];
                    else
                        *data_param[q] = *param_cur[q];
                }
                data.phi = phiD[j];
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
            real_t param_new[PARAMS];
            for (int p = 0; p < PARAMS; ++p) {
                if (optimize_param[p])
                    param_new[p] = *param_cur[p] - lambda / F_deriv2[p] * F_deriv[p];
                else
                    param_new[p] = *param_cur[p];
            }
            real_t F_new = 0;
            for (int j = 0; j < Dnum; ++j) {
                data.phi = phiD[j];
                for (int p = 0; p < PARAMS; ++p)
                    *data_param[p] = param_new[p];
                u_new[j] = calc_u(data, config);
                F_new += w[j] * pow(u_new[j] - uD[j], 2);
            }
            if (F_new <= F) {
                for (int p = 0; p < PARAMS; ++p)
                    *param_cur[p] = param_new[p];
                F = F_new;
                for (int j = 0; j < Dnum; ++j)
                    u_cur[j] = u_new[j];
                break;
            }
            else {
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
        pstr << "F = " << F << "\n\n";
        ++lambda_not_changed;
        if (lambda_not_changed == config.lambda_threshold) {
            lambda_not_changed = 0;
            lambda *= 2;
            pstr << "lambda increased = " << lambda << "\n\n";
        }
    } // for iter

    std::ofstream fout_u("plot.txt");
    for (int i = 0; i < Dnum; ++i)
        fout_u << phiD[i] << "   " << u_cur[i] << "   " << uD[i] << "\n";

    std::ofstream fout("output_param.txt");
    fout.precision(10);
    fout << "T0 = " << data0.T0 << "\n";
    fout << "D = " << data0.D << "\n";
    fout << "Q/cp = " << data0.Q_div_cp << "\n";
    fout << "nu = " << data0.nu << "\n\n";
    for (int p = 0; p < PARAMS; ++p)
        fout << param_name[p] << " = " << *param_cur[p] << "\n";
    fout << "\nlambda = " << lambda << "\n";
    fout << "iterations = " << input_param.prev_iterations + input_param.steps << "\n\n";
    fout << "steps = " << input_param.steps;

    return 0;
}
