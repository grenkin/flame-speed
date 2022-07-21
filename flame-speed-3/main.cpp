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

bool accept_params_intervals (
    const std::vector<Interval>& intervals,
    const ModelParameters& model_parameters,
    const std::vector<ExperimentalData>& experimental_data,
    const Config& config,
    const sign_t u_deriv_sign[PARAMS_NUM],
    const sign_t u_deriv2_sign[PARAMS_NUM][PARAMS_NUM], PairStream& out,
    real_t& F_min, real_t& F_max)
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
        out << "u = " << min_ui[i] << " .. " << max_ui[i] << "\n";

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
    std::vector<real_t> min_df_dxk(PARAMS_NUM), max_df_dxk(PARAMS_NUM);
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
    const std::vector<ExperimentalData>& experimental_data, const Config& config)
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
bool union_pair_of_ranges (std::vector<std::vector<Interval>>& ranges)
{
    for (int i = 0; i < ranges.size(); i++) {
        for (int j = 0; j < ranges.size(); j++) {
            std::vector<int> diff_intervals;
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

void union_ranges (std::vector<std::vector<Interval>>& ranges)
{
    while (union_pair_of_ranges(ranges)) {}
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
    //std::ofstream fout1("output_rej.txt");
    //fout1.precision(10);

    // set bounds for enumeration
    std::vector<int> bounds(PARAMS_NUM);
    for (int p = 0; p < PARAMS_NUM; ++p)
        bounds[p] = input_param.params_ranges[p].num;

    Enumeration enumeration(bounds);
    std::vector<std::vector<Interval>> accepted_ranges, rejected_ranges;
    while (!enumeration.end) {
        // enumerate (i0, i1, i2, i3, i4): 0 <= ip < params_ranges[p].num
        std::vector<Interval> intervals(PARAMS_NUM);  // intervals of parameters values
        for (int p = 0; p < PARAMS_NUM; ++p) {
            int i_p = enumeration.v[p];
            intervals[p] = Interval(input_param.params_ranges[p].left
                    + range_step[p] * i_p,
                input_param.params_ranges[p].left
                    + range_step[p] * (i_p + 1)
            );
        }

        real_t F_min, F_max;
        bool accept_intervals = accept_params_intervals(
            intervals, input_param.model_parameters, experimental_data, config,
            u_deriv_sign, u_deriv2_sign, pstr, F_min, F_max);
        pstr << "F_min >= " << F_min << "   F_max <= " << F_max << "\n";

        if (accept_intervals && F_min < input_param.F_min) {
            for (int p = 0; p < PARAMS_NUM; ++p) {
                pstr << PARAMS_NAMES[p] << " = ";
               // fout << PARAMS_NAMES[p] << " = ";
                pstr << intervals[p].left << " .. "
                    << intervals[p].right << "\n";
               // fout << intervals[p].left << " .. "
               //     << intervals[p].right << "\n";
            }
            pstr << "\n";
            //fout << "\n";
            accepted_ranges.push_back(intervals);
        }
        else {
            pstr << "Intervals not accepted\n\n";
            /*for (int p = 0; p < PARAMS_NUM; ++p) {
                fout1 << PARAMS_NAMES[p] << " = ";
                fout1 << intervals[p].left << " .. "
                    << intervals[p].right << "\n";
            }
            fout1 << "\n";*/
            rejected_ranges.push_back(intervals);
        }

        enumeration.next();
    }

    union_ranges(accepted_ranges);
    union_ranges(rejected_ranges);

    fout << "Accepted intervals:\n\n";
    for (int i = 0; i < accepted_ranges.size(); ++i) {
        std::vector<Interval> intervals = accepted_ranges[i];
        for (int p = 0; p < PARAMS_NUM; ++p) {
            fout << PARAMS_NAMES[p] << " = ";
            fout << intervals[p].left << " .. "
                << intervals[p].right << "\n";
        }
        fout << "\n";
    }
    fout << "\n\nRejected intervals:\n\n";
    for (int i = 0; i < rejected_ranges.size(); ++i) {
        std::vector<Interval> intervals = rejected_ranges[i];
        for (int p = 0; p < PARAMS_NUM; ++p) {
            fout << PARAMS_NAMES[p] << " = ";
            fout << intervals[p].left << " .. "
                << intervals[p].right << "\n";
        }
        fout << "\n";
    }

    return 0;
}
