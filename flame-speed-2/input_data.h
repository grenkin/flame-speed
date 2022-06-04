#ifndef INPUT_DATA_H_INCLUDED
#define INPUT_DATA_H_INCLUDED

const std::string file_input_param = "input_param.txt";
const std::string file_config = "config.txt";

typedef long double real_t;

struct calc_u_Data {
    real_t T0, D, nu;
    real_t A, E_div_R, alpha, beta, n;
    real_t phi, Q_div_cp;
};
// можно ввести класс CalcFlameSpeed с методом вычисления скорости

struct ParamRange {
    real_t left, right;
    int num;  // number of intervals that [left, right] separates into
};

struct InputParam {
    calc_u_Data data;
    ParamRange A, E_div_R, alpha, beta, n;
};

struct Config {
    real_t delta_A, delta_E_div_R, delta_alpha, delta_beta, delta_n;
    int N_trapezoid;
    real_t u_eps, max_u, u_init, max_speed;
};

InputParam read_input_param ();
Config read_config ();

#endif // INPUT_DATA_H_INCLUDED
