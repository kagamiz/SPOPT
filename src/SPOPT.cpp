#include "SPOPT.hpp"

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "Usage : SPOPT config_file_name" << std::endl;
        exit(1);
    }

    SPOPT::ProblemData problemData(argv[1]);

    problemData.ConstructSDP();

    return 0;
}