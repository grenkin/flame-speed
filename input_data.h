#ifndef INPUT_DATA_H_INCLUDED
#define INPUT_DATA_H_INCLUDED

const std::string file_input_param = "input_param.txt";
const std::string file_config = "config.txt";

typedef long double real_t;

struct calc_u_Data {
    real_t T0, D, Q_div_cp, nu;
    real_t A, E_div_R, alpha, beta, n;
    real_t phi;
};

struct InputParam {
    calc_u_Data data;
    real_t lambda_init;
    int steps;
    int prev_iterations;
};

struct Config {
    bool optimize_A, optimize_E_div_R, optimize_alpha, optimize_beta, optimize_n;
    real_t delta_A, delta_E_div_R, delta_alpha, delta_beta, delta_n;
    real_t lambda_threshold, lambda_decr_max;
    int N_trapezoid;
    real_t u_eps, max_u, u_init;
};

InputParam read_input_param ();
Config read_config ();

#endif // INPUT_DATA_H_INCLUDED
