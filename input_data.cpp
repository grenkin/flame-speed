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

struct ParamWrongBool {
    std::string param;
    ParamWrongBool (std::string param)
        : param(param) {}
};

void get_double_param (const po::variables_map &vm, const char *param, real_t &val)
{
    if (vm.count(param))
        val = vm[param].as<double>();
    else
        throw ParamIsNotSet(param);
}

void get_pos_double_param (const po::variables_map &vm, const char *param, real_t &val)
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

void get_bool_param (const po::variables_map &vm, const char *param, bool &val)
{
    if (vm.count(param)) {
        int val1 = vm[param].as<int>();
        if (val != 0 && val != 1)
            throw ParamWrongBool(param);
        val = (bool)val1;
    }
    else
        throw ParamIsNotSet(param);
}

void print_error (const std::string msg)
{
    std::cerr << msg << "\n";
    getch();
    exit(1);
}

InputParam read_input_param ()
{
    try {
        po::options_description desc("Input parameters");
        desc.add_options()
            ("T0", po::value<double>(), "T0")
            ("D", po::value<double>(), "D")
            ("nu", po::value<double>(), "nu")
            ("A", po::value<double>(), "A")
            ("E/R", po::value<double>(), "E/R")
            ("alpha", po::value<double>(), "alpha")
            ("beta", po::value<double>(), "beta")
            ("n", po::value<double>(), "n")
            ("lambda", po::value<double>(), "lambda")
            ("prev_iterations", po::value<int>(), "prev_iterations")
            ("steps", po::value<int>(), "steps")
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
        get_pos_double_param(vm, "T0", input_param.data.T0);
        get_pos_double_param(vm, "D", input_param.data.D);
        get_pos_double_param(vm, "nu", input_param.data.nu);
        get_pos_double_param(vm, "A", input_param.data.A);
        get_pos_double_param(vm, "E/R", input_param.data.E_div_R);
        get_double_param(vm, "alpha", input_param.data.alpha);
        get_double_param(vm, "beta", input_param.data.beta);
        get_double_param(vm, "n", input_param.data.n);
        get_pos_double_param(vm, "lambda", input_param.lambda_init);
        get_int_param(vm, "prev_iterations", input_param.prev_iterations);
        get_int_param(vm, "steps", input_param.steps);
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
            ("optimize_A", po::value<int>(), "optimize_A")
            ("optimize_E/R", po::value<int>(), "optimize_E/R")
            ("optimize_alpha", po::value<int>(), "optimize_alpha")
            ("optimize_beta", po::value<int>(), "optimize_beta")
            ("optimize_n", po::value<int>(), "optimize_n")
            ("delta_A", po::value<double>(), "delta_A")
            ("delta_E/R", po::value<double>(), "delta_E/R")
            ("delta_alpha", po::value<double>(), "delta_alpha")
            ("delta_beta", po::value<double>(), "delta_beta")
            ("delta_n", po::value<double>(), "delta_n")
            ("lambda_threshold", po::value<int>(), "lambda_threshold")
            ("lambda_decr_max", po::value<int>(), "lambda_decr_max")
            ("N_trapezoid", po::value<int>(), "N_trapezoid")
            ("u_eps", po::value<double>(), "u_eps")
            ("u_init", po::value<double>(), "u_init")
            ("max_u", po::value<double>(), "max_u")
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
        get_bool_param(vm, "optimize_A", config.optimize_A);
        get_bool_param(vm, "optimize_E/R", config.optimize_E_div_R);
        get_bool_param(vm, "optimize_alpha", config.optimize_alpha);
        get_bool_param(vm, "optimize_beta", config.optimize_beta);
        get_bool_param(vm, "optimize_n", config.optimize_n);
        get_pos_double_param(vm, "delta_A", config.delta_A);
        get_pos_double_param(vm, "delta_E/R", config.delta_E_div_R);
        get_pos_double_param(vm, "delta_alpha", config.delta_alpha);
        get_pos_double_param(vm, "delta_beta", config.delta_beta);
        get_pos_double_param(vm, "delta_n", config.delta_n);
        get_int_param(vm, "lambda_threshold", config.lambda_threshold);
        get_int_param(vm, "lambda_decr_max", config.lambda_decr_max);
        get_int_param(vm, "N_trapezoid", config.N_trapezoid);
        get_pos_double_param(vm, "u_eps", config.u_eps);
        get_pos_double_param(vm, "u_init", config.u_init);
        get_pos_double_param(vm, "max_u", config.max_u);
        return config;
    }
    catch (ParamIsNotSet e) {
        print_error("Parameter \"" + e.param + "\" is not specified");
    }
    catch (ParamIsNotPositive e) {
        print_error("Parameter \"" + e.param + "\" is not positive");
    }
    catch (ParamWrongBool e) {
        print_error("Wrong boolean parameter \"" + e.param + "\"");
    }
    catch(std::exception& e) {
        print_error(std::string("error: ") + e.what());
    }
    catch(...) {
        print_error("Exception of unknown type!");
    }
}
