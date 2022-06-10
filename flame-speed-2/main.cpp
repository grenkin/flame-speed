#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <cstdlib>
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

enum sign_t {PLUS, MINUS};

int sign (sign_t s)
{
    return s == PLUS ? 1 : -1;
}

struct Interval {
    real_t left, right;
    bool empty;

    Interval ()
    {}

    Interval (real_t l, real_t r)
        : left(l), right(r), empty(false)
    {}
};

Interval intersect_intervals (Interval a, Interval b)
{
    if (a.left > b.right || b.left > a.right) {
        Interval ans;
        ans.empty = true;
        return ans;
    }
    else {
        return Interval(fmax(a.left, b.left), fmin(a.right, b.right));
    }
}

bool inside_the_interval (real_t x, Interval interval)
{
    return x >= interval.left && x <= interval.right;
}

int main (void)
{
    std::cout << "Identification of reaction rate parameters in a combustion model\n\n";
    // Read input data
    InputParam input_param = read_input_param();
    Config config = read_config();
    // Read burnt temperature data
    BurntTemperature burnt_temperature("burnt_temperature.txt");
    // Read flame speed data
    std::ifstream fspeed("flame_speed.txt");
    if (!fspeed) {
        std::cerr << "Can not open input file: flame_speed.txt\n";
        getch();
        exit(1);
    }
    int data_size;
    fspeed >> data_size;
    std::vector<ExperimentalData> experimental_data(data_size);
    for (int j = 0; j < data_size; ++j) {
        real_t phi, v;
        fspeed >> phi >> v;
        experimental_data[j] = ExperimentalData(phi, v, burnt_temperature,
            input_param.model_parameters);
    }

    std::ofstream flog("log.txt");
    flog.precision(10);
    std::cout.precision(10);
    PairStream pstr(flog);

    // deltas are steps for numerical differentiation
    real_t **delta; // [ModelParametersToFind::PARAMS_NUM];
    delta = config.delta;

    // data contains parameters to optimize
    ModelParametersToFind& data = input_param.model_parameters_to_find;

    // Next, i = 1..m, j, k = 1..n where n = 5, m = data_size
    // signs of derivatives of u_i(x)
    sign_t u_deriv_sign[ModelParametersToFind::PARAMS_NUM] =
        {PLUS, MINUS, MINUS, MINUS, PLUS};
    // signs of second derivatives of u_i(x)
    sign_t u_deriv2_sign[ModelParametersToFind::PARAMS_NUM]
        [ModelParametersToFind::PARAMS_NUM] =
        {
            {MINUS, MINUS, MINUS, MINUS, PLUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {PLUS, MINUS, MINUS, MINUS, PLUS}
        };

    // steps of parameters ranges
    real_t range_step[ModelParametersToFind::PARAMS_NUM];
    for (int p = 0; p < ModelParametersToFind::PARAMS_NUM; ++p) {
        range_step[p] = (input_param.params_ranges[p].right
            - input_param.params_ranges[p].left)
            / input_param.params_ranges[p].num;
    }

    std::ofstream fout("output.txt");
    fout.precision(10);

    // TODO: Реализовать класс для перебора всех кортежей

    for (int i_A = 0; i_A < input_param.A.num; ++i_A) {
        for (int i_E_div_R = 0; i_E_div_R < input_param.E_div_R.num;
            ++i_E_div_R)
        {
            for (int i_alpha = 0; i_alpha < input_param.alpha.num; ++i_alpha) {
                for (int i_beta = 0; i_beta < input_param.beta.num; ++i_beta) {
                    for (int i_n = 0; i_n < input_param.n.num; ++i_n) {
                        real_t ll[PARAMS_NUM], rr[PARAMS_NUM];
                        for (int flag = 0; flag <= 1; ++flag) {
                       /* Interval interval_A(
                            input_param.A.left + r_A * i_A,
                            input_param.A.left + r_A * (i_A + 1));
                        Interval interval_E_div_R(
                            input_param.E_div_R.left + r_E_div_R * i_E_div_R,
                            input_param.E_div_R.left
                                + r_E_div_R * (i_E_div_R + 1)
                        );
                        Interval interval_alpha(
                            input_param.alpha.left + r_alpha * i_alpha,
                            input_param.alpha.left + r_alpha * (i_alpha + 1));
                        Interval interval_beta(
                            input_param.beta.left + r_beta * i_beta,
                            input_param.beta.left + r_beta * (i_beta + 1));
                        Interval interval_n(
                            input_param.n.left + r_n * i_n,
                            input_param.n.left + r_n * (i_n + 1));
                        */

                        Interval intervals[PARAMS_NUM] = {
                            Interval(input_param.A.left + r_A * i_A,
                                input_param.A.left + r_A * (i_A + 1)),
                            Interval(input_param.E_div_R.left + r_E_div_R * i_E_div_R,
                                input_param.E_div_R.left
                                + r_E_div_R * (i_E_div_R + 1)),
                            Interval(input_param.alpha.left + r_alpha * i_alpha,
                                input_param.alpha.left + r_alpha * (i_alpha + 1)),
                            Interval(input_param.beta.left + r_beta * i_beta,
                                input_param.beta.left + r_beta * (i_beta + 1)),
                            Interval(input_param.n.left + r_n * i_n,
                                input_param.n.left + r_n * (i_n + 1))
                        };

                        pstr << "A = " << intervals[0].left << " - "
                            << intervals[0].right << "\n";
                        pstr << "E/R = " << intervals[1].left << " - "
                            << intervals[1].right << "\n";
                        pstr << "alpha = " << intervals[2].left << " - "
                            << intervals[2].right << "\n";
                        pstr << "beta = " << intervals[3].left << " - "
                            << intervals[3].right << "\n";
                        pstr << "n = " << intervals[4].left << " - "
                            << intervals[4].right << "\n";

                        // Далее координаты точек в 5-мерном пространстве
                        //    задаются со знаком
                        //    (чтобы все u_i возрастали по всем x_j)
                        // точка с минимальными скоростями
                        real_t x_left[PARAMS_NUM];
                        // точка с максимальными скоростями
                        real_t x_right[PARAMS_NUM];
                        // x_min_dudxk[k] - точка с минимальными производными
                        //    скоростей по x_k (k - первый индекс)
                        real_t x_min_dudxk[PARAMS_NUM][PARAMS_NUM];
                        // x_max_dudxk[k] - точка с максимальными производными
                        //    скоростей по x_k (k - первый индекс)
                        real_t x_max_dudxk[PARAMS_NUM][PARAMS_NUM];
                        for (int param = 0; param < PARAMS_NUM; ++param) {
                            // TODO: чтобы не писать каждый раз индекс param,
                            // реализовать функцию, которая будет вычислять
                            // все 4 значения
                            if (u_deriv_sign[param] == PLUS) {
                                x_left[param] = intervals[param].left;
                                x_right[param] = intervals[param].right;
                            }
                            else if (u_deriv_sign[param] == MINUS) {
                                x_left[param] = -intervals[param].right;
                                x_right[param] = -intervals[param].left;
                            }
                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                // if (u_deriv2_sign[k][param] == PLUS) - это
                                //    если бы не было замены переменных x_j на -x_j
                                if (sign(u_deriv2_sign[k][param])
                                    * sign(u_deriv_sign[k])
                                    * sign(u_deriv_sign[param]) == 1)
                                {
                                    x_min_dudxk[k][param] = x_left[param];
                                    x_max_dudxk[k][param] = x_right[param];
                                }
                                else {
                                    x_min_dudxk[k][param] = x_right[param];
                                    x_max_dudxk[k][param] = x_left[param];
                                }
                            }
                        }

                        // maximum and minimum flame speed
                        real_t max_ui[data_size], min_ui[data_size];
                        // maximum and minimum derivatives of the flame speed
                        real_t max_dui_dxk[data_size][PARAMS_NUM],
                            min_dui_dxk[data_size][PARAMS_NUM];
                        bool too_large_speed = false;
                        for (int i = 0; i < data_size; ++i) {
                            data.phi = phi_data[i];
                            data.Q_div_cp = Q_div_cp_data[i];

                            for (int p = 0; p < PARAMS_NUM; ++p)
                                *data_param[p] = x_left[p] * sign(u_deriv_sign[p]);
                            try {
                                min_ui[i] = calc_u(data, config);
                            }
                            catch (umax_achieved) {
                                min_ui[i] = config.max_u;
                            }

                            for (int p = 0; p < PARAMS_NUM; ++p)
                                *data_param[p] = x_right[p] * sign(u_deriv_sign[p]);
                            try {
                                max_ui[i] = calc_u(data, config);
                            }
                            catch (umax_achieved) {
                                max_ui[i] = config.max_u;
                            }

                          /*  pstr << "\nmin u[" << i+1 << "] = " << min_ui[i]
                                << ",  max u[" << i+1 << "] = " << max_ui[i] << "\n";
*/
                            if (min_ui[i] > config.max_speed) {
                                too_large_speed = true;
                                pstr << "Too large speed\n\n";
                                break;
                            }

                            // TODO: реализовать отдельную функцию вычисления
                            //    производной

                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                for (int p = 0; p < PARAMS_NUM; ++p) {
                                    *data_param[p] = x_min_dudxk[k][p]
                                        * sign(u_deriv_sign[p]);
                                }
                                real_t u0 = calc_u(data, config);
                                *data_param[k] += *delta[k];
                                real_t u1 = calc_u(data, config);
                                *data_param[k] -= *delta[k];
                                min_dui_dxk[i][k] = (u1 - u0) / *delta[k]
                                    * sign(u_deriv_sign[k]);

                                for (int p = 0; p < PARAMS_NUM; ++p) {
                                    *data_param[p] = x_max_dudxk[k][p]
                                        * sign(u_deriv_sign[p]);
                                }
                                u0 = calc_u(data, config);
                                *data_param[k] += *delta[k];
                                u1 = calc_u(data, config);
                                *data_param[k] -= *delta[k];
                                max_dui_dxk[i][k] = (u1 - u0) / *delta[k]
                                    * sign(u_deriv_sign[k]);

                             /*   pstr << "du[" << i+1 << "]/d" << param_name[k]
                                    << " = " << min_dui_dxk[i][k] << " - "
                                    << max_dui_dxk[i][k] << "\n"; */
                                // TODO: поставить минусы перед du/dx
                            }
                        }

                        if (too_large_speed)
                            // this parameters range doesn't take interest to us
                            continue;

                     /*   fout << "A = " << intervals[0].left << " - "
                            << intervals[0].right << "\n";
                        fout << "E/R = " << intervals[1].left << " - "
                            << intervals[1].right << "\n";
                        fout << "alpha = " << intervals[2].left << " - "
                            << intervals[2].right << "\n";
                        fout << "beta = " << intervals[3].left << " - "
                            << intervals[3].right << "\n";
                        fout << "n = " << intervals[4].left << " - "
                            << intervals[4].right << "\n\n"; */

                        // найти a[i][j], b[i]; посчитать eps1, eps2, eps
                        // найти множество S_eps

                        real_t a[data_size][PARAMS_NUM], b[data_size];
                       /* for (int i = 0; i < data_size; ++i) {
                            // !!! программа плохо спроектирована,
                            //   потому что эти присваивания легко забыть
                            data.phi = phi_data[i];
                            data.Q_div_cp = Q_div_cp_data[i];

                            for (int j = 0; j < PARAMS_NUM; ++j) {
                               // a[i][j] = min_dui_dxk[i][j];

                               //a[i][j] = max_dui_dxk[i][j];

                               a[i][j] = max_dui_dxk[i][j];

                              // a[i][j] = (max_dui_dxk[i][j] + min_dui_dxk[i][j]) / 2;

                               /* real_t u0 = min_dui_dxk[i][j];
                                for (int p = 0; p < PARAMS_NUM; ++p) {
                                    *data_param[p] = x_left[p]
                                        * sign(u_deriv_sign[p]);
                                }
                                *data_param[j] = x_right[j]
                                    * sign(u_deriv_sign[j]);
                                real_t u1 = calc_u(data, config);
                                a[i][j] = (u1 - u0) / (x_right[j] - x_left[j]);
                                */
                           /* }
                            // вычисляем b[i] в точке x_left,
                            //    в которой u = min_ui[i]
                            b[i] = min_ui[i];
                            for (int j = 0; j < PARAMS_NUM; ++j)
                                b[i] -= a[i][j] * x_left[j];
                        }*/

                        real_t eps0 = 0.0;
                        for (int i = 0; i < data_size; ++i) {
                            for (int j = 0; j < PARAMS_NUM; ++j)
                                //a[i][j] = (max_ui[i] - min_ui[i]) / (x_right[j] - x_left[j]);
                              if (!flag)
                                //a[i][j] = min_dui_dxk[i][j];  // так лучше?
                                a[i][j] = (max_dui_dxk[i][j] + min_dui_dxk[i][j]) / 2;
                              else
                                 a[i][j] = max_dui_dxk[i][j];
                              //  a[i][j] = min_dui_dxk[i][j];
                                // a[i][j] = (min_dui_dxk[i][j] + max_dui_dxk[i][j]) / 2;
                            // вычисляем b[i] в точке x_left,
                            //    в которой u = min_ui[i]
                            b[i] = min_ui[i] + eps0;
                            for (int j = 0; j < PARAMS_NUM; ++j)
                                b[i] -= a[i][j] * x_left[j];
                        }

                        // calculate maximum of the linear approximation of u
//                        for (int p = 0; p < PARAMS_NUM; ++p)
                            //*data_param[p] = x_right[p] * sign(u_deriv_sign[p]);
                        for (int i = 0; i < data_size; ++i) {
                            real_t sum = b[i];
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                sum += a[i][j] * x_right[j];
                            }
                            //pstr << "max_lin_u[i] = " << calc_u(data, config) << "\n";
                            pstr << "max_lin_u[i] = " << sum << "\n";
                        }

                        real_t sum_aik_i[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            sum_aik_i[k] = 0;
                            for (int i = 0; i < data_size; ++i)
                                sum_aik_i[k] += a[i][k];
                        }

                        real_t eps1[data_size], eps2[data_size][PARAMS_NUM];
                        for (int i = 0; i < data_size; ++i) {
                            real_t sum = 0;
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                sum += (max_dui_dxk[i][j] - min_dui_dxk[i][j])
                                    * (x_right[j] - x_left[j]);
                            }
                            eps1[i] = sum;

                            real_t max_u_minus_v =
                                fmax(fabs(max_ui[i] - v_data[i]),
                                    fabs(min_ui[i] - v_data[i]));
                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                eps2[i][k] = max_u_minus_v
                                    * (max_dui_dxk[i][k] - min_dui_dxk[i][k]) / a[i][k];
                            }
                        }

                      /*  for (int i = 0; i < data_size; ++i) {
                            pstr << "eps1[" << i + 1 << "] = " << eps1[i] << "\n";
                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                pstr << "eps2[" << i + 1 << ", " << param_name[k]
                                    << "] = " << eps2[i][k] << "\n";
                            }
                        } */

                        real_t eps1_max = 0.0, eps2_max = 0.0;
                        for (int i = 0; i < data_size; ++i) {
                            eps1_max = fmax(eps1_max, eps1[i]);
                            for (int k = 0; k < PARAMS_NUM; ++k)
                                eps2_max = fmax(eps2_max, eps2[i][k]);
                        }

                        real_t eps = eps1_max + eps2_max;
                        pstr << "eps = " << eps << "\n";
                      /*  for (int k = 0; k < PARAMS_NUM; ++k) {
                            pstr << "sum_aik_i[" << param_name[k] << "] = "
                                << sum_aik_i[k] << "\n";
                        } */

                        // TODO: реализовать отдельную функцию для решения СЛАУ

                        real_t l[PARAMS_NUM], r[PARAMS_NUM];  // S_eps range
                        gsl_matrix *mat = gsl_matrix_alloc(PARAMS_NUM, PARAMS_NUM);
                        gsl_vector *vec = gsl_vector_alloc(PARAMS_NUM);
                        gsl_vector *x_vec = gsl_vector_alloc(PARAMS_NUM);
                        // fill matrix of a linear system
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                real_t sum = 0;
                                for (int i = 0; i < data_size; ++i)
                                    sum += a[i][j] * a[i][k];
                                gsl_matrix_set(mat, k, j, sum);  // mat[k, j] = sum
                               // pstr << "mat[" << k << ", " << j << "] = " << sum << " ";
                            }
                            pstr << "\n";
                        }

                        // find Cholesky decomposition of the matrix
                        gsl_linalg_cholesky_decomp(mat);

                    /*    // fill right-hand side of the linear system
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sum = eps * sum_aik_i[k];
                            for (int i = 0; i < data_size; ++i) {
                                sum -= a[i][k] * (b[i] - v_data[i]);
                            }
                            gsl_vector_set(vec, k, sum);
                        }

                        gsl_linalg_cholesky_solve(mat, vec, x_vec);
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            r[k] = gsl_vector_get(x_vec, k);

                        // fill right-hand side of a new linear system
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sum = - eps * sum_aik_i[k];
                            for (int i = 0; i < data_size; ++i) {
                                sum -= a[i][k] * (b[i] - v_data[i]);
                            }
                            gsl_vector_set(vec, k, sum);
                        }

                      //  gsl_linalg_cholesky_decomp(mat);
                        gsl_linalg_cholesky_solve(mat, vec, x_vec);
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            l[k] = gsl_vector_get(x_vec, k);
*/

                        // fill right-hand side of a new linear system
                        real_t lin_sol[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sum = 0;
                            for (int i = 0; i < data_size; ++i) {
                                sum -= a[i][k] * (b[i] - v_data[i]);
                            }
                            gsl_vector_set(vec, k, sum);
                        }


                       // gsl_linalg_cholesky_decomp(mat);
                        gsl_linalg_cholesky_solve(mat, vec, x_vec);
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            lin_sol[k] = gsl_vector_get(x_vec, k);

                        // вычисляем b[i] в точке x_right,
                        //    в которой u = max_ui[i]
                        for (int i = 0; i < data_size; ++i) {
                            b[i] = max_ui[i] - eps0;
                            for (int j = 0; j < PARAMS_NUM; ++j)
                                b[i] -= a[i][j] * x_right[j];
                        }

                        real_t lin_sol2[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sum = 0;
                            for (int i = 0; i < data_size; ++i) {
                                sum -= a[i][k] * (b[i] - v_data[i]);
                            }
                            gsl_vector_set(vec, k, sum);
                        }
                        gsl_linalg_cholesky_solve(mat, vec, x_vec);
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            lin_sol2[k] = gsl_vector_get(x_vec, k);

                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            l[k] = lin_sol[k];
                            r[k] = lin_sol2[k];

                            if (!flag) {
                                ll[k] = l[k];
                                rr[k] = r[k];
                            }
                            else {
                                ll[k] = fmin(fmin(ll[k], l[k]), r[k]);
                                rr[k] = fmax(fmax(rr[k], r[k]), l[k]);
                            }
                        }


                        /*
                        // deviations_rhs[k]
                        real_t deviations_rhs[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            deviations_rhs[k] = eps * sum_aik_i[k];
                        }
                        // basis_sol[j][k], k=1..n содержит решение СЛАУ,
                        // где правая часть - базисный вектор e[j]
                        real_t basis_sol[PARAMS_NUM][PARAMS_NUM];
                        for (int j = 0; j < PARAMS_NUM; ++j) {
                            // берем правую часть с базисным вектором e[j]

                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                if (k == j)
                                   // gsl_vector_set(vec, k, 1.0);
                                    gsl_vector_set(vec, k, deviations_rhs[k]);
                                else
                                    gsl_vector_set(vec, k, 0.0);
                            }
                            gsl_linalg_cholesky_solve(mat, vec, x_vec);
                            for (int k = 0; k < PARAMS_NUM; ++k)
                                basis_sol[j][k] = gsl_vector_get(x_vec, k);
                            std::cout << "basis_sol[" << j << "]:\n";
                            for (int k = 0; k < PARAMS_NUM; ++k)
                                std::cout << basis_sol[j][k] << " ";
                            std::cout << "\n";
                        }

                        // max deviations from lin_sol[k]
                        real_t deviations_sol[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sum = 0.0;
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                sum += fabs(basis_sol[j][k]); //* deviations_rhs[k]);
                            }
                            deviations_sol[k] = sum;
                            l[k] = lin_sol[k] - deviations_sol[k];
                            r[k] = lin_sol[k] + deviations_sol[k];
                        }
                        */

                     /*   // calculate linear objective function
                        real_t lin_obj_func = 0.0;
                        for (int i = 0; i < data_size; ++i) {
                            real_t term = 0;
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                term += a[i][j] * lin_sol[j];
                            }
                            term += b[i] - v_data[i];
                            lin_obj_func += 0.5 * pow(term, 2);
                        }
                       */ // calculate nonlinear objective function
                        real_t nonlin_obj_func = 0.0, nonlin_obj_func2 = 0.0;
                        for (int i = 0; i < data_size; ++i) {
                            data.phi = phi_data[i];
                            data.Q_div_cp = Q_div_cp_data[i];
                            for (int p = 0; p < PARAMS_NUM; ++p) {
                                *data_param[p] = lin_sol[p]
                                    * sign(u_deriv_sign[p]);
                            }
                            real_t u;
                            try {
                                u = calc_u(data, config);
                            }
                            catch (umax_achieved) {
                                u = config.max_u;
                            }
                            real_t term = u - v_data[i];
                            nonlin_obj_func += 0.5 * pow(term, 2);

                            for (int p = 0; p < PARAMS_NUM; ++p) {
                                *data_param[p] = lin_sol2[p]
                                    * sign(u_deriv_sign[p]);
                            }
                            try {
                                u = calc_u(data, config);
                            }
                            catch (umax_achieved) {
                                u = config.max_u;
                            }
                            term = u - v_data[i];
                            nonlin_obj_func2 += 0.5 * pow(term, 2);
                        }

                        gsl_vector_free(x_vec);
                        gsl_vector_free(vec);
                        gsl_matrix_free(mat);

                       /* for (int k = 0; k < PARAMS_NUM; ++k) {
                            real_t sp1 = 0, sp2 = 0;
                            for (int j = 0; j < PARAMS_NUM; ++j) {
                                // вычисляем элемент матрицы A[k, j]
                                real_t sum = 0;
                                for (int i = 0; i < data_size; ++i) {
                                    sum += a[i][j] * a[i][k];
                                }
                                sp1 += sum * l[j];
                                sp2 += sum * r[j];
                            }
                            std::cout << sp1 << " " << sp2 << " " << eps * sum_aik_i[k] << "\n";
                        }
                        std::cout << "\n"; */

                        // if l[k] > r[k] then swap them
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            if (l[k] > r[k]) {
                                real_t tmp = l[k];
                                l[k] = r[k];
                                r[k] = tmp;
                             //   pstr << "swap\n";
                            }
                        }

                        pstr << "l: ";
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            pstr << l[k] << " ";
                        pstr << "\n r: ";
                        for (int k = 0; k < PARAMS_NUM; ++k)
                            pstr << r[k] << " ";
                        pstr << "\n\n";



                       /* for (int j = 0; j < PARAMS_NUM; ++j) {
                            pstr << "lin_sol[" << param_name[j] << "] = "
                                << lin_sol[j] * sign(u_deriv_sign[j])
                                << "   lin_sol2[" << param_name[j] << "] = "
                                << lin_sol2[j] * sign(u_deriv_sign[j]) << "\n";
                        }
                        pstr << "\n"; */

                        // TODO: вывести пересечение отдельно по каждой координате -
                        //    одна координата в одну строчку

                        // проверяем, что решение линейной задачи
                        //    принадлежит рассматриваемому диапазону
                        bool inside_the_range = true;
                        for (int j = 0; j < PARAMS_NUM; ++j) {
                            inside_the_range = inside_the_range
                                && inside_the_interval(
                                    lin_sol[j] * sign(u_deriv_sign[j]),
                                        intervals[j]);
                        }
                     /*   if (!inside_the_range)
                            pstr << "The solution of linear problem is not "
                                "inside the range.\n\n";
                        }
                        else {
                            pstr << "Empty intersection\n\n";
                        } */

                    /*    // calculate objective function
                        real_t sum = 0;
                        for (int i = 0; i < data_size; ++i) {
                            real_t term = 0;
                            for (int j = 0; j < PARAMS_NUM; ++j)
                                term += (max_dui_dxk[i][j] - min_dui_dxk[i][j])
                                    * (x_right[j] - x_left[j]);
                            term *= fmax(fabs(max_ui[i] - v_data[i]),
                                fabs(min_ui[i] - v_data[i]));
                            sum += term;
                        }
                        pstr << "Difference between nonlinear and linear "
                            "objective function is less than " << sum << "\n";*/
                       // pstr << "Linear objective function is "
                         pstr << "\nNonlinear objective function is "
                            << nonlin_obj_func << "\n";
                         pstr << "\nNonlinear objective function 2 is "
                            << nonlin_obj_func2 << "\n";

                    } // flag

                        // create intervals [l, r]
                        Interval S_eps[PARAMS_NUM];
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            if (u_deriv_sign[k] == PLUS) {
                                S_eps[k].left = ll[k];
                                S_eps[k].right = rr[k];
                            }
                            else {
                                S_eps[k].right = -ll[k];
                                S_eps[k].left = -rr[k];
                            }
                        }
                        // intersect [l, r] with [x_left, x_right]
                      /*  Interval sol[PARAMS_NUM];
                        bool empty_sol = false;
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            sol[k] = intersect_intervals(
                                S_eps[k], intervals[k]);
                            empty_sol = empty_sol || sol[k].empty;
                        } */
                        pstr << "S_eps :\n";
                        for (int k = 0; k < PARAMS_NUM; ++k) {
                            pstr << param_name[k] << "_left = "
                                << S_eps[k].left << "\n"
                                << param_name[k] << "_right = "
                                << S_eps[k].right << "\n";
                        }
                        pstr << "\n";

                        // output intersection
                     /*   if (!empty_sol) {
                            pstr << "Intersection:\n";
                            if (!inside_the_range)
                                fout << "not inside the range\n";
                            for (int k = 0; k < PARAMS_NUM; ++k) {
                                //fout
                                pstr << param_name[k] << "_left = "
                                    << sol[k].left << "\n "
                                    << param_name[k] << "_right = "
                                    << sol[k].right << "\n";
                            }
                            fout << "\n";
                            pstr << "\n";
                        }*/
                    }
                }
            }
        }
    }

    return 0;
}
