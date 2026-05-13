#include "solver/AdvectionDiffusionSolver1D.hpp"
#include "solver/AdvectionDiffusionSolver2D.hpp"
#include "solver/BaseSolver.hpp"
#include "solver/Burgers2D.hpp"
#include <iostream>
#include <cmath>


int main(int argc, char* argv[]) {

    // Catch the config file

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n"
        << "Config file missing\n";
        return 1;
    }

    ConfigParser cfg(argv[1]);
    std::string solver_type = cfg.getString("SOLVER");

    std::cout << "========================================" << std::endl;
    std::cout << "              SOLVER MVP_CFD            " << std::endl;
    std::cout << "========================================" << std::endl << std::endl;

    std::cout << "SOLVER: " << solver_type << std::endl;


    try{
        std::unique_ptr<BaseSolver> solver;
        if (solver_type == "BURGERS_2D") {
            solver = std::make_unique<Burgers2D>(cfg);
        } 
        else if (solver_type == "ADVECTION_DIFFUSION_2D") {
            solver = std::make_unique<AdvectionDiffusionSolver2D>(cfg);
        }
        else if (solver_type == "ADVECTION_DIFFUSION_1D") {
            solver = std::make_unique<AdvectionDiffusionSolver1D>(cfg);
        }
        else {
            std::cerr << "Error: solver not recognized" << std::endl;
            return 1;
        }

        // Setting initial conditions
        solver->setInitialCondition();

        // Setting boundary conditions
        solver->setBoundaryConditions();
        
        //Solve
        solver->solve();

        std::cout << "\n========================================" << std::endl;
        std::cout << "          Simulation complete           " << std::endl;
        std::cout << "========================================" << std::endl;
    }
    catch (const GeneralRuntimeError& e){
        std::cerr << "\n[FATAL ERROR] " << e.what() << '\n';
        return EXIT_FAILURE; 
    }
    return EXIT_SUCCESS;
}
