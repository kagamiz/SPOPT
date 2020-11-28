#include "SPOPT.hpp"

#include <yaml-cpp/yaml.h>
#include <boost/program_options.hpp>

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description opt("Allowed options");

    std::string problemDataConfigFileName, solverConfigFileName;

    opt.add_options()
        ("help",                                                                             					                       "print this help message")
        ("problem_data_config_file,p",	po::value<std::string>(&problemDataConfigFileName)->default_value("examples/problem_data.yml"), "path of the problem data configuration file")
        ("solver_config_file,s",	    po::value<std::string>(&solverConfigFileName)->default_value("examples/solver.yml"),	           "path of the solver configuration file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opt), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << opt << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    SPOPT::ProblemData problemData(problemDataConfigFileName);

    YAML::Node solverConfig = YAML::LoadFile(solverConfigFileName);
    std::map<std::string, SPOPT::Solver *> solvers;
    solvers["DualLagrangianSolver"] = new SPOPT::DualLagrangianSolver();
    auto solverData = std::move(solvers[solverConfig["solverType"].as<std::string>()]);
    solverData->LoadConfig(solverConfigFileName);

    problemData.ConstructSDP();
    solverData->Solve(problemData);

    return 0;
}