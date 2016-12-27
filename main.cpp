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

    calc_u_Data data0 = input_param.data;
    real_t A_init = data0.A;
    real_t E_div_R_init = data0.E_div_R;
    real_t alpha_init = data0.alpha;
    real_t beta_init = data0.beta;
    real_t n_init = data0.n;
    real_t lambda_init = input_param.lambda_init;
    real_t delta_A = config.delta_A;
    real_t delta_E_div_R = config.delta_E_div_R;
    real_t delta_alpha = config.delta_alpha;
    real_t delta_beta = config.delta_beta;
    real_t delta_n = config.delta_n;

    std::ofstream flog("log.txt");
    flog.precision(10);
    std::cout.precision(10);
    PairStream pstr(flog);

    calc_u_Data data = data0;
    real_t A_cur, E_div_R_cur, alpha_cur, beta_cur, n_cur;
    real_t F;
    std::vector<real_t> u_cur(Dnum), u_new(Dnum);
    real_t lambda;

    lambda = lambda_init;
    A_cur = A_init;
    E_div_R_cur = E_div_R_init;
    alpha_cur = alpha_init;
    beta_cur = beta_init;
    n_cur = n_init;
    F = 0;
    for (int j = 0; j < Dnum; ++j) {
        data.phi = phiD[j];
        data.A = A_cur;
        data.E_div_R = E_div_R_cur;
        data.alpha = alpha_cur;
        data.beta = beta_cur;
        data.n = n_cur;
        u_cur[j] = calc_u(data, config);
        F += w[j] * pow(u_cur[j] - uD[j], 2);
    }
    pstr << "A = " << A_cur << "\nE/R = " << E_div_R_cur << "\nalpha = " << alpha_cur << "\nbeta = " << beta_cur << "\nn = " << n_cur << "\nF = " << F << "\nlambda = " << lambda << "\n\n";

    //the number of subsequent iterations at which lambda was not changed
    int lambda_not_changed = 0;
    for (int iter = 1; iter <= input_param.steps; ++iter) {
        pstr << "iteration " << input_param.prev_iterations + iter << "\n";
        real_t F_A = 0, F_E_div_R = 0, F_alpha = 0, F_beta = 0, F_n = 0;
        real_t d_A = 0, d_E_div_R = 0, d_alpha = 0, d_beta = 0, d_n = 0;
        for (int j = 0; j < Dnum; ++j) {
            data.phi = phiD[j];

            data.A = A_cur + delta_A;
            data.E_div_R = E_div_R_cur;
            data.alpha = alpha_cur;
            data.beta = beta_cur;
            data.n = n_cur;
            real_t u_A = calc_u(data, config);

            data.A = A_cur;
            data.E_div_R = E_div_R_cur + delta_E_div_R;
            data.alpha = alpha_cur;
            data.beta = beta_cur;
            data.n = n_cur;
            real_t u_E_div_R = calc_u(data, config);

            data.A = A_cur;
            data.E_div_R = E_div_R_cur;
            data.alpha = alpha_cur + delta_alpha;
            data.beta = beta_cur;
            data.n = n_cur;
            real_t u_alpha = calc_u(data, config);

            data.A = A_cur;
            data.E_div_R = E_div_R_cur;
            data.alpha = alpha_cur;
            data.beta = beta_cur + delta_beta;
            data.n = n_cur;
            real_t u_beta = calc_u(data, config);

            data.A = A_cur;
            data.E_div_R = E_div_R_cur;
            data.alpha = alpha_cur;
            data.beta = beta_cur;
            data.n = n_cur + delta_n;
            real_t u_n = calc_u(data, config);

            real_t deriv_A = (u_A - u_cur[j]) / delta_A;
            real_t deriv_E_div_R = (u_E_div_R - u_cur[j]) / delta_E_div_R;
            real_t deriv_alpha = (u_alpha - u_cur[j]) / delta_alpha;
            real_t deriv_beta = (u_beta - u_cur[j]) / delta_beta;
            real_t deriv_n = (u_n - u_cur[j]) / delta_n;
            flog << "delta u (A) = " << u_A - u_cur[j] << "   delta u (E/R) = " << u_E_div_R - u_cur[j] <<
                "   delta u (alpha) = " << u_alpha - u_cur[j] << "   delta u (beta) = " << u_beta - u_cur[j] <<
                "   delta u (n) = " << u_n - u_cur[j] << "\n";
            F_A += w[j] * (u_cur[j] - uD[j]) * deriv_A;
            F_E_div_R += w[j] * (u_cur[j] - uD[j]) * deriv_E_div_R;
            F_alpha += w[j] * (u_cur[j] - uD[j]) * deriv_alpha;
            F_beta += w[j] * (u_cur[j] - uD[j]) * deriv_beta;
            F_n += w[j] * (u_cur[j] - uD[j]) * deriv_n;
            d_A += w[j] * pow(deriv_A, 2);
            d_E_div_R += w[j] * pow(deriv_E_div_R, 2);
            d_alpha += w[j] * pow(deriv_alpha, 2);
            d_beta += w[j] * pow(deriv_beta, 2);
            d_n += w[j] * pow(deriv_n, 2);
        }
        flog << "F_A = " << F_A << "   F_E/R = " << F_E_div_R <<
            "   F_alpha = " << F_alpha << "   F_beta = " << F_beta <<
            "   F_n = " << F_n << "\n";

        int lambda_decr = 0; // the number of subsequent decreases of lambda
        while (1) {
            real_t A_new = A_cur, E_div_R_new = E_div_R_cur,
                alpha_new = alpha_cur, beta_new = beta_cur, n_new = n_cur;
            if (config.optimize_A)
                A_new = A_cur - lambda / d_A * F_A;
            if (config.optimize_E_div_R)
                E_div_R_new = E_div_R_cur - lambda / d_E_div_R * F_E_div_R;
            if (config.optimize_alpha)
                alpha_new = alpha_cur - lambda / d_alpha * F_alpha;
            if (config.optimize_beta)
                beta_new = beta_cur - lambda / d_beta * F_beta;
            if (config.optimize_n)
                n_new = n_cur - lambda / d_n * F_n;
            real_t F_new = 0;
            for (int j = 0; j < Dnum; ++j) {
                data.phi = phiD[j];
                data.A = A_new;
                data.E_div_R = E_div_R_new;
                data.alpha = alpha_new;
                data.beta = beta_new;
                data.n = n_new;
                u_new[j] = calc_u(data, config);
                F_new += w[j] * pow(u_new[j] - uD[j], 2);
            }
            if (F_new <= F) {
                A_cur = A_new;
                E_div_R_cur = E_div_R_new;
                alpha_cur = alpha_new;
                beta_cur = beta_new;
                n_cur = n_new;
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
        pstr << "A = " << A_cur << "\nE/R = " << E_div_R_cur << "\nalpha = " << alpha_cur << "\nbeta = " << beta_cur << "\nn = " << n_cur << "\nF = " << F << "\n\n";
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

    return 0;
}
