#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <cstdlib>
#include "input_data.h"
#include "calc_u.h"
#include "enumeration.h"

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
    {
        if (l > r)
            empty = true;
    }
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

// Solve A^T (Ay + b) = 0 with interval A and b
// A is m*n-matrix, b is m-vector
// The result is n-vector
std::vector<Interval> solve_interval_linear_system (
    int m, int n, const std::vector<std::vector<Interval>>& a,
    const std::vector<Interval>& b)
{
    std::vector<Interval> sol(n, Interval(0.0, 0.0));
    Enumeration e_a(n * (n + 1) / 2);
    // real_t A[n][n];
    gsl_matrix *mat = gsl_matrix_alloc(n, n);
    gsl_vector *vec = gsl_vector_alloc(n);
    gsl_vector *y_vec = gsl_vector_alloc(n);
    while (!e_a.end) {
        int counter = 0;
        // fill matrix of a linear system
        for (int j = 0; j < n; j++) {
            for (int k = j; k < n; k++) {
                real_t sum = 0;
                if (e_a.v[counter] == 0) {
                    for (int i = 0; i < m; i++)
                        sum += a[i][j].left * a[i][k].left;
                }
                else if (e_a.v[counter] == 1) {
                    for (int i = 0; i < m; i++)
                        sum += a[i][j].right * a[i][k].right;
                }
                counter++;
                // A[j][k] = A[k][j] = sum;
                gsl_matrix_set(mat, j, k, sum);
                gsl_matrix_set(mat, k, j, sum);
            }
        }
        // find Cholesky decomposition of the matrix
        gsl_linalg_cholesky_decomp(mat);

        Enumeration e_b(n);
        while (!e_b.end) {
            // fill right-hand side of the linear system
            for (int k = 0; k < n; k++) {
                real_t sum = 0;
                for (int i = 0; i < m; i++) {
                    real_t b_val;
                    if (e_b.v[k] == 0)
                        b_val = b[i].left;
                    else if (e_b.v[k] == 1)
                        b_val = b[i].right;
                    sum -= (a[i][k].left + a[i][k].right) / 2 * b_val;
                }
                gsl_vector_set(vec, k, sum);
            }

            // solve the linear system
            gsl_linalg_cholesky_solve(mat, vec, y_vec);
            for (int k = 0; k < n; k++) {
                real_t y_k = gsl_vector_get(y_vec, k);
                sol[k].left = fmin(sol[k].left, y_k);
                sol[k].right = fmax(sol[k].right, y_k);
                std::cout << y_k << "   ";
            }
            std::cout << "\n";

            e_b.next();
        }
        e_a.next();
    }

    /*for (int j = 0; j < n; ++j)
        sol[j] = Interval(0, 0);*/

    return sol;
}

std::vector<Interval> process_params_intervals (
    const Interval intervals[PARAMS_NUM],
    const ModelParameters& model_parameters,
    const std::vector<ExperimentalData>& experimental_data,
    const Config& config,
    const sign_t u_deriv_sign[PARAMS_NUM],
    const sign_t u_deriv2_sign[PARAMS_NUM][PARAMS_NUM], PairStream& out
)
{
    std::vector<Interval> new_intervals(PARAMS_NUM);
    for (int p = 0; p < PARAMS_NUM; ++p) {
        out << PARAMS_NAMES[p] << " = "
            << intervals[p].left << " .. " << intervals[p].right << "\n";
    }
    out << "\n";

    // data contains parameters to optimize
    ModelParametersToFind data;

    // x[p] = params[p] * deriv_sign[p] so that flame speed increases
    //    as any parameter x[p] increases
    // x_left and x_right are points with minimal and maximal speeds resp.
    real_t x_left[PARAMS_NUM], x_right[PARAMS_NUM];
    // x_min_dudxk[k] and x_max_dudxk[k] are points with
    //    minimal and maximal du_i/dx_k
    real_t x_min_dudxk[PARAMS_NUM][PARAMS_NUM],
        x_max_dudxk[PARAMS_NUM][PARAMS_NUM];
    // calculate the variables declared above
    for (int p = 0; p < PARAMS_NUM; ++p) {
        if (u_deriv_sign[p] == PLUS) {
            x_left[p] = intervals[p].left;
            x_right[p] = intervals[p].right;
        }
        else if (u_deriv_sign[p] == MINUS) {
            x_left[p] = -intervals[p].right;
            x_right[p] = -intervals[p].left;
        }
    }
    for (int p = 0; p < PARAMS_NUM; ++p) {
        for (int k = 0; k < PARAMS_NUM; ++k) {
            // we might write
            // if (u_deriv2_sign[k][param] == PLUS)
            // if it were not change of variables x_j to -x_j
            if (sign(u_deriv2_sign[k][p])
                * sign(u_deriv_sign[k])
                * sign(u_deriv_sign[p]) == 1)
            {
                x_min_dudxk[k][p] = x_left[p];
                x_max_dudxk[k][p] = x_right[p];
            }
            else {
                x_min_dudxk[k][p] = x_right[p];
                x_max_dudxk[k][p] = x_left[p];
            }
        }
    }

    int data_size = experimental_data.size();
    // maximal and minimal flame speed
    real_t max_ui[data_size], min_ui[data_size];
    // maximal and minimal derivatives of the flame speed
    real_t max_dui_dxk[data_size][PARAMS_NUM],
        min_dui_dxk[data_size][PARAMS_NUM];
    for (int i = 0; i < data_size; ++i) {
        // calculate minimal and maximal flame speed
        for (int p = 0; p < PARAMS_NUM; ++p)
            *data.params[p] = x_left[p] * sign(u_deriv_sign[p]);
        try {
            min_ui[i] = calc_u(
                model_parameters, data, experimental_data[i], config);
        }
        catch (umax_achieved) {
            min_ui[i] = config.max_u;
        }

        for (int p = 0; p < PARAMS_NUM; ++p)
            *data.params[p] = x_right[p] * sign(u_deriv_sign[p]);
        try {
            max_ui[i] = calc_u(
                model_parameters, data, experimental_data[i], config);
        }
        catch (umax_achieved) {
            max_ui[i] = config.max_u;
        }

        out << "i = " << i << "\n";
        out << "u = " << min_ui[i] << " .. " << max_ui[i] << "\n";

        // calculate minimal and maximal derivatives of the flame speed
        for (int k = 0; k < PARAMS_NUM; ++k) {
            for (int p = 0; p < PARAMS_NUM; ++p)
                *data.params[p] = x_min_dudxk[k][p] * sign(u_deriv_sign[p]);
            try {
                min_dui_dxk[i][k] = calc_deriv_u(
                    k, model_parameters, data, experimental_data[i], config)
                    * sign(u_deriv_sign[k]);
            }
            catch (umax_achieved) {
                min_dui_dxk[i][k] = 0;
            }

            for (int p = 0; p < PARAMS_NUM; ++p)
                *data.params[p] = x_max_dudxk[k][p] * sign(u_deriv_sign[p]);
            try {
                max_dui_dxk[i][k] = calc_deriv_u(
                    k, model_parameters, data, experimental_data[i], config)
                    * sign(u_deriv_sign[k]);
            }
            catch (umax_achieved) {
                max_dui_dxk[i][k] = 0;
            }

            out << (u_deriv_sign[k] == PLUS ? "" : "-")
                << "du/d" << PARAMS_NAMES[k] << " = "
                << min_dui_dxk[i][k] << " .. " << max_dui_dxk[i][k] << "\n";
        }

        // check condition that max_dui_dxk - min_dui_dxk < min_dui_dxk
        real_t sum = 0;
        for (int k = 0; k < PARAMS_NUM; ++k) {
            sum += (max_dui_dxk[i][k] - min_dui_dxk[i][k])
                * (x_right[k] - x_left[k]);
        }
        out << "grad u_i * (x_right - x_left) = " << sum << " ; "
            << "u_max - u_min = " << max_ui[i] - min_ui[i];
        if (sum < max_ui[i] - min_ui[i])
            out << "  condition fulfilled\n";
        else
            out << "  condition not fulfilled\n";
    }
    out << "\n";

    // an interval matrix and an interval vector
    std::vector<std::vector<Interval>>
        a(data_size, std::vector<Interval>(PARAMS_NUM));
    std::vector<Interval> b(data_size);
    // calculate a and b
    for (int i = 0; i < data_size; ++i) {
        for (int j = 0; j < PARAMS_NUM; ++j) {
            a[i][j] = Interval(min_dui_dxk[i][j], max_dui_dxk[i][j]);
            b[i] = Interval(min_ui[i] - experimental_data[i].v,
                max_ui[i] - experimental_data[i].v);
        }
    }

    // Solve A^T (Ay + b) = 0
    std::vector<Interval> y = solve_interval_linear_system(
        data_size, PARAMS_NUM, a, b);
    // calculate [x_right + y_min, x_left + y_max]
    for (int k = 0; k < PARAMS_NUM; ++k) {
        real_t l = x_right[k] + y[k].left, r = x_left[k] + y[k].right;
        if (l > r)
            new_intervals[k].empty = true;
        else {
            if (u_deriv_sign[k] == PLUS)
                new_intervals[k] = Interval(l, r);
            else
                new_intervals[k] = Interval(-r, -l);
        }
    }

    return new_intervals;
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

    // Next, i = 1..m, j, k = 1..n where n = 5, m = data_size
    // signs of derivatives of u_i(x)
    sign_t u_deriv_sign[PARAMS_NUM] =
        {PLUS, MINUS, MINUS, MINUS, PLUS};
    // signs of second derivatives of u_i(x)
    sign_t u_deriv2_sign[PARAMS_NUM][PARAMS_NUM] =
        {
            {MINUS, MINUS, MINUS, MINUS, PLUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {MINUS, PLUS, PLUS, PLUS, MINUS},
            {PLUS, MINUS, MINUS, MINUS, PLUS}
        };

    // steps of parameters ranges
    real_t range_step[PARAMS_NUM];
    for (int p = 0; p < PARAMS_NUM; ++p) {
        range_step[p] = (input_param.params_ranges[p].right
            - input_param.params_ranges[p].left)
            / input_param.params_ranges[p].num;
    }

    std::ofstream fout("output.txt");
    fout.precision(10);

    // set bounds for enumeration
    std::vector<int> bounds(PARAMS_NUM);
    for (int p = 0; p < PARAMS_NUM; ++p)
        bounds[p] = input_param.params_ranges[p].num;

    Enumeration enumeration(bounds);
    while (!enumeration.end) {
        // enumerate (i0, i1, i2, i3, i4): 0 <= ip < params_ranges[p].num
        Interval intervals[PARAMS_NUM];  // intervals of parameters values
        for (int p = 0; p < PARAMS_NUM; ++p) {
            int i_p = enumeration.v[p];
            intervals[p] = Interval(input_param.params_ranges[p].left
                    + range_step[p] * i_p,
                input_param.params_ranges[p].left
                    + range_step[p] * (i_p + 1)
            );
        }

        std::vector<Interval> new_intervals = process_params_intervals(
            intervals, input_param.model_parameters, experimental_data, config,
            u_deriv_sign, u_deriv2_sign, pstr);

        pstr << "New intervals:\n";
        for (int p = 0; p < PARAMS_NUM; ++p) {
            pstr << PARAMS_NAMES[p] << " = ";
            if (new_intervals[p].empty)
                pstr << "empty\n";
            else
                pstr << new_intervals[p].left << " .. "
                    << new_intervals[p].right << "\n";
        }
        pstr << "\n";

        enumeration.next();
    }

    return 0;
}
