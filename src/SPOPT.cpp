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
        ("all_real_eigenpairs,a",       po::value<std::string>(&fileName),                  "show all real eigenpair in ascending order of the eigenvalue")
        ("get_all_real_eigenpairs,g",                                                       "enumerate all real eigenpairs")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opt), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << opt << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    SPOPT::ProblemData problemData(problemDataConfigFileName);

    if (!vm.count("get_all_real_eigenpairs")) {
        problemData.ConstructSDP();
    }
    
    if (vm.count("solver_config_file")) {
        YAML::Node solverConfig = YAML::LoadFile(solverConfigFileName);
        std::map<std::string, SPOPT::Solver *> solvers;
        solvers["DualLagrangianSolver"] = new SPOPT::DualLagrangianSolver();
        solvers["HSDESolver"] = new SPOPT::HSDESolver();
        solvers["MOSEKSolver"] = new SPOPT::MOSEKSolver();
        solvers["SCSSolver"] = new SPOPT::SCSSolver();

        auto solverData = solvers[solverConfig["solverType"].as<std::string>("MOSEKSolver")];
        solverData->LoadConfig(solverConfigFileName);

        if (!vm.count("get_all_real_eigenpairs")) {
            auto [infimum, tms] = solverData->Solve(problemData);
            std::ofstream ofs;
            ofs.open("tms.txt");
            for (auto &elem : tms) {
                ofs << elem << std::endl;
            }
            ofs.close();

            std::vector<std::vector<double>> sols = problemData.ExtractSolutionsFrom(tms, true, infimum);

            std::cout << "Extracted " << sols.size() << " solution(s)." << std::endl;
            for (auto sol : sols) {
                problemData.Analyze(sol);
            }

            for (auto it = solvers.begin(); it != solvers.end(); it++) {
                delete it->second;
            }
        }
        else {
            solverData->EnumerateStationaryPoints(problemDataConfigFileName);
        }
    }

    #ifdef __BUILD_WITH_MATLAB__
    if (vm.count("mat")) {
        problemData.OutputMatFile(fileName);
    }
    #endif

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

        std::vector<std::vector<double>> sols = problemData.ExtractSolutionsFrom(tms);

        std::cout << "Extracted " << sols.size() << " solution(s)." << std::endl;
        for (auto sol : sols) {
            problemData.Analyze(sol);
        }
    }

    if (vm.count("view_problem")) {
        problemData.ShowPOP();
    }

    if (vm.count("all_real_eigenpairs")) {
        std::ifstream fin;
        fin.open(fileName);
        if (fin.fail()) {
            std::cerr << "opening file " << fileName << " failed !" << std::endl;
            exit(EXIT_FAILURE);
        }
        int n, m;
        fin >> n >> m;
        std::vector<std::pair<double, std::vector<double>>> realEigenVectors(m);
        for (int i = 0; i < m; i++) {
            std::vector<double> eigenvector(n);
            for (auto &elem : eigenvector) {
                fin >> elem;
            }
            realEigenVectors.emplace_back(problemData.Apply(eigenvector), eigenvector);
        }
        fin.close();
        std::sort(realEigenVectors.begin(), realEigenVectors.end());

        std::cout << m << std::endl;
        for (int i = 0; i < std::min(m, 10); i++) {
            printf("%lf : ", realEigenVectors[i].first);
            for (int j = 0; j < n; j++) {
                printf("%lf, ", realEigenVectors[i].second[j]);
            }
            printf("\n");
        }
    }

    return 0;
}
