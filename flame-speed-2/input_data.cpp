#include <iostream>
#include <fstream>
#include <conio.h>
#include <boost/program_options.hpp>
#include "input_data.h"

namespace po = boost::program_options;

struct ParamIsNotSet {
    std::string param;
    ParamIsNotSet (std::string param)
        : param(param) {}
};

struct ParamIsNotPositive {
    std::string param;
    ParamIsNotPositive (std::string param)
        : param(param) {}
};

void get_double_param (const po::variables_map &vm, const char *param,
    real_t &val)
{
    if (vm.count(param))
        val = vm[param].as<double>();
    else
        throw ParamIsNotSet(param);
}

void get_positive_double_param (const po::variables_map &vm,
    const char *param, real_t &val)
{
    if (vm.count(param)) {
        val = vm[param].as<double>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_int_param (const po::variables_map &vm, const char *param, int &val)
{
    if (vm.count(param)) {
        val = vm[param].as<int>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_range (const po::variables_map &vm, const char *param,
    ParamRange &range)
{
    std::string param_name = std::string(param) + "_left";
    get_double_param(vm, param_name.c_str(), range.left);
    param_name = std::string(param) + "_right";
    get_double_param(vm, param_name.c_str(), range.right);
    param_name = std::string(param) + "_num";
    get_int_param(vm, param_name.c_str(), range.num);
}

void get_positive_range (const po::variables_map &vm, const char *param,
    ParamRange &range)
{
    std::string param_name = std::string(param) + "_left";
    get_positive_double_param(vm, param_name.c_str(), range.left);
    param_name = std::string(param) + "_right";
    get_positive_double_param(vm, param_name.c_str(), range.right);
    param_name = std::string(param) + "_num";
    get_int_param(vm, param_name.c_str(), range.num);
}

void print_error (const std::string msg)
{
    std::cerr << msg << "\n";
    getch();
    exit(1);
}

real_t calc_z (const ModelParameters& model_parameters)
{
    return model_parameters.nu * (1 + 3.762 * 28. / 32.);
}

InputParam read_input_param ()
{
    try {
        po::options_description desc("Input parameters");
        desc.add_options()
            ("T0", po::value<double>(), "T0")
            ("D", po::value<double>(), "D")
            ("nu", po::value<double>(), "nu")
            ("A_left", po::value<double>(), "A_left")
            ("A_right", po::value<double>(), "A_right")
            ("A_num", po::value<int>(), "A_num")
            ("E/R_left", po::value<double>(), "E/R_left")
            ("E/R_right", po::value<double>(), "E/R_right")
            ("E/R_num", po::value<int>(), "E/R_num")
            ("alpha_left", po::value<double>(), "alpha_left")
            ("alpha_right", po::value<double>(), "alpha_right")
            ("alpha_num", po::value<int>(), "alpha_num")
            ("beta_left", po::value<double>(), "beta_left")
            ("beta_right", po::value<double>(), "beta_right")
            ("beta_num", po::value<int>(), "beta_num")
            ("n_left", po::value<double>(), "n_left")
            ("n_right", po::value<double>(), "n_right")
            ("n_num", po::value<int>(), "n_num")
        ;

        po::variables_map vm;
        std::ifstream ifs(file_input_param.c_str());
        if (!ifs)
            print_error("Can not open input file: " + file_input_param);
        else {
            po::store(parse_config_file(ifs, desc), vm);
            po::notify(vm);
        }

        InputParam input_param;
        get_positive_double_param(vm, "T0", input_param.model_parameters.T0);
        get_positive_double_param(vm, "D", input_param.model_parameters.D);
        get_positive_double_param(vm, "nu", input_param.model_parameters.nu);
        get_positive_range(vm, "A", input_param.params_ranges[0]);
        get_positive_range(vm, "E/R", input_param.params_ranges[1]);
        get_range(vm, "alpha", input_param.params_ranges[2]);
        get_range(vm, "beta", input_param.params_ranges[3]);
        get_range(vm, "n", input_param.params_ranges[4]);
        return input_param;
    }
    catch (ParamIsNotSet e) {
        print_error("Parameter \"" + e.param + "\" is not specified");
    }
    catch (ParamIsNotPositive e) {
        print_error("Parameter \"" + e.param + "\" is not positive");
    }
    catch(std::exception& e) {
        print_error(std::string("error: ") + e.what());
    }
    catch(...) {
        print_error("Exception of unknown type!");
    }
}

Config read_config ()
{
    try {
        po::options_description desc("Configuration parameters");
        desc.add_options()
            ("delta_A", po::value<double>(), "delta_A")
            ("delta_E/R", po::value<double>(), "delta_E/R")
            ("delta_alpha", po::value<double>(), "delta_alpha")
            ("delta_beta", po::value<double>(), "delta_beta")
            ("delta_n", po::value<double>(), "delta_n")
            ("N_trapezoid", po::value<int>(), "N_trapezoid")
            ("u_eps", po::value<double>(), "u_eps")
            ("u_init", po::value<double>(), "u_init")
            ("max_u", po::value<double>(), "max_u")
            ("max_speed", po::value<double>(), "max_speed")
        ;

        po::variables_map vm;
        std::ifstream ifs(file_config.c_str());
        if (!ifs)
            print_error("Can not open input file: " + file_config);
        else {
            po::store(parse_config_file(ifs, desc), vm);
            po::notify(vm);
        }

        Config config;
        get_positive_double_param(vm, "delta_A", config.delta_A);
        get_positive_double_param(vm, "delta_E/R", config.delta_E_div_R);
        get_positive_double_param(vm, "delta_alpha", config.delta_alpha);
        get_positive_double_param(vm, "delta_beta", config.delta_beta);
        get_positive_double_param(vm, "delta_n", config.delta_n);
        get_int_param(vm, "N_trapezoid", config.N_trapezoid);
        get_positive_double_param(vm, "u_eps", config.u_eps);
        get_positive_double_param(vm, "u_init", config.u_init);
        get_positive_double_param(vm, "max_u", config.max_u);
        get_positive_double_param(vm, "max_speed", config.max_speed);
        return config;
    }
    catch (ParamIsNotSet e) {
        print_error("Parameter \"" + e.param + "\" is not specified");
    }
    catch (ParamIsNotPositive e) {
        print_error("Parameter \"" + e.param + "\" is not positive");
    }
    catch(std::exception& e) {
        print_error(std::string("error: ") + e.what());
    }
    catch(...) {
        print_error("Exception of unknown type!");
    }
}
