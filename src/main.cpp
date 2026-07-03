#include "solver/AdvectionDiffusionSolver1D.hpp"
#include "solver/AdvectionDiffusionSolver2D.hpp"
#include "solver/BaseSolver.hpp"
#include "solver/Burgers2D.hpp"
#include "solver/IncNS2D.hpp"
#include <iostream>
#include <cstring>
#include <cmath>


int main(int argc, char* argv[]) {

    // Catch the config file
    bool dryRunFlag = false;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n"
        << "Config file missing\n";
        return 1;
    }
    if (argc >= 3 && std::strcmp(argv[2], "dry") == 0) dryRunFlag = true;


    ConfigParser cfg(argv[1]);
    std::string solver_type = cfg.getString("SOLVER");

    std::cout << "========================================\n";
    std::cout << "              SOLVER MVP_CFD            \n";
    std::cout << "========================================\n\n";

    std::cout << "SOLVER: " << solver_type << "\n";


    try{
        std::unique_ptr<BaseSolver> solver;
        if (solver_type == "INCOMPRESSIBLE_NS_2D") {
            solver = std::make_unique<IncNS2D>(cfg);
        } 
        else if (solver_type == "BURGERS_2D") {
            solver = std::make_unique<Burgers2D>(cfg);
        } 
        else if (solver_type == "ADVECTION_DIFFUSION_2D") {
            solver = std::make_unique<AdvectionDiffusionSolver2D>(cfg);
        }
        else if (solver_type == "ADVECTION_DIFFUSION_1D") {
            solver = std::make_unique<AdvectionDiffusionSolver1D>(cfg);
        }
        else {
            std::cerr << "Error: solver not recognized\n";
            return 1;
        }

        // Setting initial conditions
        solver->setInitialCondition();

        // Setting boundary conditions
        solver->setBoundaryConditions();
        
        //Solve
        if (!dryRunFlag){
            solver->solve();

            std::cout << "\n========================================\n";
            std::cout << "          Simulation complete           \n";
            std::cout << "========================================\n";
        }
        else{
            std::cout << "\n========================================\n";
            std::cout << "            Dry run complete            \n";
            std::cout << "========================================\n";
        }

    }
    catch (const GeneralRuntimeError& e){
        std::cerr << "\n[FATAL ERROR] " << e.what() << '\n';
        return EXIT_FAILURE; 
    }
    return EXIT_SUCCESS;
}
