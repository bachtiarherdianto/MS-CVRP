/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#include "Solver.h"

int main(int argc, char **argv) {
    if(argc == 1) {
        std::cout << "Missing input instance.\n\n";
        display_help();
        exit(1);
    }

    if(argc == 2 && std::string(argv[1]) == "--help") {
        display_help();
        exit(0);
    }

    auto parameters = Parameters(argv[1]);
    for(int n = 2; n < argc; n+=2) {
        auto token = std::string(argv[n]);

        if(n + 1 >= argc) {
            std::cout << "Missing value for '" << token << "'.\n\n";
            display_help();
            exit(0);
        }

        parameters.set(token, std::string(argv[n+1]));
    }

    std::cout << "File name: " << parameters.get_file_input() << std::endl;

    Solver solve = Solver(parameters);

    std::cout << "Time budget: " << solve.runTimeMax << " seconds\n";
    std::cout << "Seed: " << parameters.get_seed() << std::endl;

    srand(parameters.get_seed());
    solve.run();

    std::cout << "\nEnd of optimization"  << std::endl;
    std::cout << "Total Cost: " << solve.globalBest.SolObjcost << std::endl;

    return 0;
}