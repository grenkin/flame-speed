#include <cmath>
#include <cstdlib>
#include "input_data.h"
#include "calc_u.h"
#include "enumeration.h"

using namespace std;

class PairStream {
public:
    ostream &file;

    PairStream (ostream &f)
        : file(f)
    {}

    template<typename T>
    PairStream& operator<< (const T &t)
    {
        cout << t;
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
    if (a.empty || b.empty || a.left > b.right - 1e-12
        || b.left > a.right - 1e-12)
    {
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

real_t calc_func (
    const ModelParametersToFind& data, const ModelParameters& model_parameters,
    const vector<ExperimentalData>& experimental_data,
    const Config& config, vector<real_t>& u)
{
    int m = experimental_data.size();
    real_t F = 0.0;
    for (int i = 0; i < m; i++) {
        real_t u_i = calc_u(model_parameters, data, experimental_data[i], config);
        F += pow(u_i - experimental_data[i].v, 2);
        u[i] = u_i;
    }
    return F;
}

void calc_func_gradient (
    const ModelParametersToFind& data, const ModelParameters& model_parameters,
    const vector<ExperimentalData>& experimental_data,
    const Config& config, bool u_filled, const vector<real_t>& u,
    real_t& F, vector<real_t>& grad_F, vector<real_t>& F_deriv2)
{
    int m = experimental_data.size();
    F = 0.0;
    for (int j = 0; j < PARAMS_NUM; j++) {
        grad_F[j] = 0.0;
        F_deriv2[j] = 0.0;
    }

    for (int i = 0; i < m; i++) {
        real_t u_i;
        if (u_filled)
            u_i = u[i];
        else
            u_i = calc_u(model_parameters, data, experimental_data[i], config);
        F += pow(u_i - experimental_data[i].v, 2);
        for (int j = 0; j < PARAMS_NUM; j++) {
            real_t dui_dxj = calc_deriv_u(j, model_parameters, data,
                experimental_data[i], config, u_i);
            grad_F[j] += 2 * (u_i - experimental_data[i].v) * dui_dxj;
            F_deriv2[j] += 2 * pow(dui_dxj, 2);
        }
    }
    // F and grad_F were computed
}

template<typename T>
void print_parameters (ModelParametersToFind data, T& out)
{
    for (int p = 0; p < PARAMS_NUM; ++p)
        out << PARAMS_NAMES[p] << " = " << data.param(p) << "\n";
    out << "\n";
}

template<typename T>
void print_intervals (const vector<Interval>& intervals, T& out)
{
    for (int p = 0; p < PARAMS_NUM; ++p) {
        out << PARAMS_NAMES[p] << " = "
            << intervals[p].left << " .. " << intervals[p].right << "\n";
    }
    out << "\n";
}

void start_gradient_descent (
    const ModelParametersToFind& starting_point,
    const ModelParameters& model_parameters,
    const vector<ExperimentalData>& experimental_data,
    const Config& config,
    ModelParametersToFind& final_point, real_t& final_F,
    PairStream& out)
{
    out << "Start gradient descent\n";
    print_parameters(starting_point, out);

    ModelParametersToFind data = starting_point;

    real_t F;  // current value of the objective function
    vector<real_t> u(experimental_data.size());  // current flame speed
    // grad_F and F_deriv2 are approximations of the first
    //   and the second derivatives of F
    vector<real_t> grad_F(PARAMS_NUM), F_deriv2(PARAMS_NUM);
    bool u_filled = false;

    real_t lambda = config.gradient_descent_step_size;

    // lambda_not_changed is the number of subsequent iterations
    //   at which lambda was not changed
    int lambda_not_changed = 0;
    for (int iter = 1; iter <= config.gradient_descent_steps; ++iter) {
        out << "iteration " << iter << "\n";
        calc_func_gradient(data, model_parameters, experimental_data, config,
            u_filled, u, F, grad_F, F_deriv2);
        for (int p = 0; p < PARAMS_NUM; ++p)
            out << "F_" << PARAMS_NAMES[p] << " = " << grad_F[p] << "   ";
        out << "\n";
        // invariant: F (but not grad_F and F_deriv2, see (1)) corresponds
        //   to the results of computations for data,
        //   and if u_filled == true then u also corresponds to data

        int lambda_decr = 0;  // the number of subsequent decreases of lambda
        while (1) {
            // Calculate the new guess
            ModelParametersToFind new_data;
            for (int p = 0; p < PARAMS_NUM; ++p)
                new_data.param(p) = data.param(p)
                    - lambda / F_deriv2[p] * grad_F[p];
            // Calculate the flame speeds and the objective function for the new guess
            real_t F_new;
            try {
                F_new = calc_func(new_data, model_parameters,
                    experimental_data, config, u);
            }
            catch (umax_achieved) {
                out << "umax_achieved:\n";
                print_parameters(new_data, out);
                throw umax_achieved();
            }
            u_filled = true;
            if (F_new <= F) {
                // Accept the new guess
                data = new_data;
                F = F_new;  // (1)
                break;
            }
            else {
                // Reject the new guess and decrease lambda
                lambda_not_changed = -1;
                lambda /= 2;
                out << "lambda = " << lambda << "\n\n";
                ++lambda_decr;
                if (lambda_decr == config.lambda_decr_max) {
                    out << "lambda_decr_max achieved\n";
                    break;
                }
                if (lambda < config.lambda_min) {
                    out << "lambda_min achieved\n";
                    break;
                }
            }
        }
        if (lambda_decr == config.lambda_decr_max)
            break;
        if (lambda < config.lambda_min)
            break;
        print_parameters(data, out);
        out << "F = " << F << "\n\n";
        ++lambda_not_changed;
        // If lambda was not decreased for many times, then increase lambda;
        // if lambda was not increased for many times, then increase lambda again
        if (lambda_not_changed == config.lambda_threshold) {
            lambda_not_changed = 0;
            lambda *= 2;
            out << "lambda increased = " << lambda << "\n\n";
        }
    }

    final_point = data;
    final_F = F;
}

bool accept_params_intervals (
    const vector<Interval>& intervals,
    const ModelParameters& model_parameters,
    const vector<ExperimentalData>& experimental_data,
    const Config& config,
    const sign_t u_deriv_sign[PARAMS_NUM],
    const sign_t u_deriv2_sign[PARAMS_NUM][PARAMS_NUM], PairStream& out,
    real_t& F_min, real_t& F_max)
{
    vector<Interval> new_intervals(PARAMS_NUM);
    print_intervals(intervals, out);

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
            data.param(p) = x_left[p] * sign(u_deriv_sign[p]);
        try {
            min_ui[i] = calc_u(
                model_parameters, data, experimental_data[i], config);
        }
        catch (umax_achieved) {
            min_ui[i] = config.max_u;
            out << "umax_achieved\n";
        }

        for (int p = 0; p < PARAMS_NUM; ++p)
            data.param(p) = x_right[p] * sign(u_deriv_sign[p]);
        try {
            max_ui[i] = calc_u(
                model_parameters, data, experimental_data[i], config);
        }
        catch (umax_achieved) {
            max_ui[i] = config.max_u;
            out << "umax_achieved\n";
        }

        out << "i = " << i << "\n";
        out << "u - v = " << min_ui[i] - experimental_data[i].v
            << " .. " << max_ui[i] - experimental_data[i].v << "\n";

        // calculate minimal and maximal derivatives of the flame speed
        for (int k = 0; k < PARAMS_NUM; ++k) {
            bool use_x_left = true;
            for (int p = 0; p < PARAMS_NUM; ++p) {
                data.param(p) = x_min_dudxk[k][p] * sign(u_deriv_sign[p]);
                use_x_left = use_x_left && (x_min_dudxk[k][p] == x_left[p]);
            }
            // use_x_left = forall p : x_min_dudxk[k][p] == x_left[p]
            try {
                real_t ui = 0.0;
                if (use_x_left)
                    ui = min_ui[i];
                min_dui_dxk[i][k] = calc_deriv_u(
                    k, model_parameters, data, experimental_data[i], config, ui)
                    * sign(u_deriv_sign[k]);
            }
            catch (umax_achieved) {
                min_dui_dxk[i][k] = config.max_u;
                out << "umax_achieved\n";
            }

            bool use_x_right = true;
            for (int p = 0; p < PARAMS_NUM; ++p) {
                data.param(p) = x_max_dudxk[k][p] * sign(u_deriv_sign[p]);
                use_x_right = use_x_right && (x_max_dudxk[k][p] == x_right[p]);
            }
            // use_x_right = forall p : x_max_dudxk[k][p] == x_right[p]
            try {
                real_t ui = 0.0;
                if (use_x_right)
                    ui = max_ui[i];
                max_dui_dxk[i][k] = calc_deriv_u(
                    k, model_parameters, data, experimental_data[i], config, ui)
                    * sign(u_deriv_sign[k]);
            }
            catch (umax_achieved) {
                max_dui_dxk[i][k] = config.max_u;
                out << "umax_achieved\n";
            }

            out << (u_deriv_sign[k] == PLUS ? "" : "-")
                << "du/d" << PARAMS_NAMES[k] << " = "
                << min_dui_dxk[i][k] << " .. " << max_dui_dxk[i][k] << "\n";
        }
    }
    out << "\n";

    // calculate estimates of minimum from below
    //   and maximum of the objective function
    F_min = 0.0;
    F_max = 0.0;
    for (int i = 0; i < data_size; ++i) {
        if ((min_ui[i] - experimental_data[i].v)
            * (max_ui[i] - experimental_data[i].v) > 0)
        {
            F_min += pow(fmin(fabs(min_ui[i] - experimental_data[i].v),
                fabs(max_ui[i] - experimental_data[i].v)), 2);
        }
        F_max += pow(fmax(fabs(min_ui[i] - experimental_data[i].v),
            fabs(max_ui[i] - experimental_data[i].v)), 2);
    }

    // calculate necessary condition of minimum
    vector<real_t> min_df_dxk(PARAMS_NUM), max_df_dxk(PARAMS_NUM);
    for (int k = 0; k < PARAMS_NUM; ++k) {
        real_t sum_min = 0, sum_max = 0;
        for (int i = 0; i < data_size; ++i) {
            sum_min += (min_ui[i] - experimental_data[i].v)
                * (min_ui[i] - experimental_data[i].v > 0 ? min_dui_dxk[i][k]
                   : max_dui_dxk[i][k]);
            sum_max += (max_ui[i] - experimental_data[i].v)
                * (max_ui[i] - experimental_data[i].v > 0 ? max_dui_dxk[i][k]
                   : min_dui_dxk[i][k]);
        }
        min_df_dxk[k] = sum_min;
        max_df_dxk[k] = sum_max;
        out << (u_deriv_sign[k] == PLUS ? "" : "-")
            << "df/d" << PARAMS_NAMES[k] << " = "
            << min_df_dxk[k] << " .. " << max_df_dxk[k] << "\n";
    }
    out << "\n";

    bool accept = true;
    for (int k = 0; k < PARAMS_NUM; ++k)
        accept = accept && (min_df_dxk[k] <= 0 && 0 <= max_df_dxk[k]);

    return accept;
}

real_t calc_sigma(const ModelParameters& model_parameters,
    const ModelParametersToFind& data,
    const vector<ExperimentalData>& experimental_data, const Config& config)
{
    int m = experimental_data.size();
    real_t sum = 0.0;
    for (int i = 0; i < m; i++) {
        sum += pow(
            calc_u(model_parameters, data, experimental_data[i], config)
            - experimental_data[i].v, 2);
    }
    real_t sigma = sqrt(sum / m);
    return sigma;
}

// try to union a pair of ranges, if done return true
bool union_pair_of_ranges (vector<vector<Interval>>& ranges)
{
    for (int i = 0; i < ranges.size(); i++) {
        for (int j = 0; j < ranges.size(); j++) {
            vector<int> diff_intervals;
            for (int p = 0; p < PARAMS_NUM; p++) {
                // if (ranges[i][p] != ranges[j][p])
                if (ranges[i][p].left != ranges[j][p].left ||
                    ranges[i][p].right != ranges[j][p].right)
                {
                    diff_intervals.push_back(p);
                }
            }
            if (diff_intervals.size() == 1) {
                int p = diff_intervals[0];
                if (ranges[i][p].left == ranges[j][p].right) {
                    ranges[j][p].right = ranges[i][p].right;
                    ranges.erase(ranges.begin() + i);
                    return true;
                }
                else if (ranges[i][p].right == ranges[j][p].left) {
                    ranges[i][p].right = ranges[j][p].right;
                    ranges.erase(ranges.begin() + j);
                    return true;
                }
            }
        }
    }
    return false;
}

void union_ranges (vector<vector<Interval>>& ranges)
{
    while (union_pair_of_ranges(ranges)) {}
}

int main (void)
{
    cout << "Identification of reaction rate parameters in a combustion model\n\n";

    // Read input data
    InputParam input_param = read_input_param();
    Config config = read_config();

    // Read burnt temperature data
    BurntTemperature burnt_temperature("burnt_temperature.txt");

    // Read flame speed data
    ifstream fspeed("flame_speed.txt");
    if (!fspeed) {
        cerr << "Can not open input file: flame_speed.txt\n";
        getch();
        exit(1);
    }
    int data_size;
    fspeed >> data_size;
    vector<ExperimentalData> experimental_data(data_size);
    for (int j = 0; j < data_size; ++j) {
        real_t phi, v;
        fspeed >> phi >> v;
        experimental_data[j] = ExperimentalData(phi, v, burnt_temperature,
            input_param.model_parameters);
    }

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

    ofstream flog("log.txt");
    flog.precision(10);
    cout.precision(10);
    PairStream pstr(flog);
    ofstream fout("output.txt");
    fout.precision(10);

    vector<Interval> intervals(PARAMS_NUM);  // intervals of parameters values
    for (int p = 0; p < PARAMS_NUM; ++p) {
        intervals[p] = Interval(input_param.params_ranges[p].left,
            input_param.params_ranges[p].right);
    }

    fout << "Input range:\n";
    print_intervals(intervals, fout);
    fout.close();

    // Enumerate all corners of the range and start the gradient descent
    //   from each corner
    Enumeration enumeration(PARAMS_NUM);
    while (!enumeration.end) {
        ModelParametersToFind data, new_data;
        for (int p = 0; p < PARAMS_NUM; p++) {
            data.param(p) = (enumeration.v[p] == 0 ?
                intervals[p].left : intervals[p].right);
        }
        real_t F;
        ofstream fout("output.txt", std::ios_base::app);
        fout << "\n-------------\nStarting point:\n";
        print_parameters(data, fout);
        try {
            start_gradient_descent(data, input_param.model_parameters,
                experimental_data, config, new_data, F, pstr);
            fout << "Final point:\n";
            print_parameters(new_data, fout);
            fout << "F = " << F << "\n";
            fout << "sigma = " << calc_sigma(input_param.model_parameters,
                new_data, experimental_data, config) << "\n";
/*
            // Calculate projections of the final point to the hyperplanes
            for (int p = 0; p < PARAMS_NUM; p++) {
                ModelParametersToFind new_data_projection = new_data;
                new_data_projection.param(p) = data.param(p);
                fout.close();
                fout.open("output.txt", std::ios_base::app);
                fout << "\n-------------\nProjection:\n";
                print_parameters(new_data_projection, fout);
                start_gradient_descent(new_data_projection,
                    input_param.model_parameters,
                    experimental_data, config, new_data, F, pstr);
                fout << "Final point:\n";
                print_parameters(new_data, fout);
                fout << "F = " << F << "\n";
                fout << "sigma = " << calc_sigma(input_param.model_parameters,
                    new_data, experimental_data, config) << "\n";
            }
*/
        }
        catch (umax_achieved) {
            fout << "umax_achieved\n";
        }
        enumeration.next();
    }

    return 0;
}
