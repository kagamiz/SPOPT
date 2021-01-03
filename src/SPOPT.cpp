#include "SPOPT.hpp"

#include <yaml-cpp/yaml.h>
#include <boost/program_options.hpp>

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description opt("Allowed options");

    std::string problemDataConfigFileName = "", solverConfigFileName = "";
    std::string fileName = "";

    opt.add_options()
        ("help",                                                                            "print this help message")
        ("problem_data_config_file,p",	po::value<std::string>(&problemDataConfigFileName), "path of the problem data configuration file")
        ("solver_config_file,s",	    po::value<std::string>(&solverConfigFileName),	    "path of the solver configuration file")
        ("tssos,t",                     po::value<std::string>(&fileName),                  "converts problem data to a tssos format")
        ("mat,m",                       po::value<std::string>(&fileName),                  "converts problem data to a matlab format")
        ("extract,e",                   po::value<std::string>(&fileName),                  "extract solution(s) from a given truncated moment sequence")
        ("view_problem,v",                                                                  "shows the converted problem")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opt), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << opt << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    SPOPT::ProblemData problemData(problemDataConfigFileName);
    problemData.ConstructSDP();

    if (vm.count("solver_config_file")) {
        YAML::Node solverConfig = YAML::LoadFile(solverConfigFileName);
        std::map<std::string, SPOPT::Solver *> solvers;
        solvers["DualLagrangianSolver"] = new SPOPT::DualLagrangianSolver();
        solvers["HSDESolver"] = new SPOPT::HSDESolver();

        auto solverData = std::move(solvers[solverConfig["solverType"].as<std::string>()]);
        solverData->LoadConfig(solverConfigFileName);

        solverData->Solve(problemData);

        delete solvers["DualLagrangianSolver"];
        delete solvers["HSDESolver"];
    }

    if (vm.count("mat")) {
        problemData.OutputMatFile(fileName);
    }

    if (vm.count("tssos")) {
        problemData.OutputJuliaFile(fileName);
    }

    if (vm.count("extract")) {
        std::ifstream fin;
        fin.open(fileName);
        if (fin.fail()) {
            std::cerr << "opening file " << fileName << " failed !" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::vector<double> tms;
        double moment;
        while (fin >> moment) {
            tms.emplace_back(moment);
        }
        fin.close();

        std::vector<Eigen::VectorXd> sols = problemData.ExtractSolutionsFrom(tms);

        std::cout << "Extracted " << sols.size() << " solution(s)." << std::endl;
        for (auto sol : sols) {
            std::cout << "\n";
            std::cout << sol << std::endl;
        }
    }

    if (vm.count("view_problem")) {
        problemData.ShowPOP();
    }

    return 0;
}