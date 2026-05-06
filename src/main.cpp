#include "solver/AdvectionDiffusionSolver1D.h"
#include "solver/AdvectionDiffusionSolver2D.h"
#include "solver/BaseSolver.h"
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
        if (solver_type == "ADVECTION_DIFFUSION_2D") {
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

        //Initialising the output file;
        solver->createOutputFile();
        
        //Solve
        solver->solve();
    }
    catch (const GeneralRuntimeError& e){
        std::cerr << "\n[FATAL ERROR] " << e.what() << '\n';
        return EXIT_FAILURE; 
    }

    return EXIT_SUCCESS;

    // // Condizione iniziale: gaussiana
    // solver.setInitialCondition(gaussiana);
    
    // // Condizioni al contorno periodiche
    // solver.setBoundaryConditions(
    //     AdvectionDiffusionSolver::BoundaryCondition::PERIODIC
    // );

    // // Salva condizione iniziale
    // solver->saveToFile("results/prova_t0.dat");
    
    // // Salva risultato
    // solver->saveToFile("results/prova_t05.dat");
    // solver->printStats();

    
    /*
    // ==========================================
    // ESEMPIO 1: Avvezione pura (D=0)
    // ==========================================
    std::cout << "\n\n=== ESEMPIO 1: Avvezione pura ===" << std::endl;
    {
        int nx = 201;
        double L = 1.0;
        double c = 1.0;   // Velocità avvezione
        double D = 0.0;   // Nessuna diffusione
        
        AdvectionDiffusionSolver solver(
            nx, L, c, D,
            AdvectionDiffusionSolver::SpatialScheme::UPWIND,  // Upwind per avvezione
            AdvectionDiffusionSolver::TimeScheme::RK4
        );
        
        // Condizione iniziale: gaussiana
        solver.setInitialCondition(gaussiana);
        
        // Condizioni al contorno periodiche
        solver.setBoundaryConditions(
            AdvectionDiffusionSolver::BoundaryCondition::PERIODIC
        );
        
        // Salva condizione iniziale
        solver.saveToFile("avvezione_t0.dat");
        
        // Risolvi
        double dt = 0.001;  // CFL ~ 0.5
        double tEnd = 0.5;
        solver.solve(dt, tEnd);
        
        // Salva risultato
        solver.saveToFile("avvezione_t05.dat");
        solver.printStats();
    }
    
    // ==========================================
    // ESEMPIO 2: Diffusione pura (c=0)
    // ==========================================
    std::cout << "\n\n=== ESEMPIO 2: Diffusione pura ===" << std::endl;
    {
        int nx = 201;
        double L = 1.0;
        double c = 0.0;   // Nessuna avvezione
        double D = 0.01;  // Coefficiente diffusione
        
        AdvectionDiffusionSolver solver(
            nx, L, c, D,
            AdvectionDiffusionSolver::SpatialScheme::CENTRAL,
            AdvectionDiffusionSolver::TimeScheme::EULER_EXPLICIT
        );
        
        // Condizione iniziale: onda quadra
        solver.setInitialCondition(ondaQuadra);
        
        // Condizioni al contorno Dirichlet
        solver.setBoundaryConditions(
            AdvectionDiffusionSolver::BoundaryCondition::DIRICHLET,
            0.0,  // u(0) = 0
            0.0   // u(L) = 0
        );
        
        // Salva condizione iniziale
        solver.saveToFile("diffusione_t0.dat");
        
        // Risolvi
        double dt = 0.0001;  // Numero diffusione ~ 0.2
        double tEnd = 0.1;
        solver.solve(dt, tEnd);
        
        // Salva risultato
        solver.saveToFile("diffusione_t01.dat");
        solver.printStats();
    }
    
    // ==========================================
    // ESEMPIO 3: Avvezione-Diffusione
    // ==========================================
    std::cout << "\n\n=== ESEMPIO 3: Avvezione-Diffusione ===" << std::endl;
    {
        int nx = 201;
        double L = 1.0;
        double c = 0.5;    // Velocità avvezione
        double D = 0.005;  // Coefficiente diffusione
        
        AdvectionDiffusionSolver solver(
            nx, L, c, D,
            AdvectionDiffusionSolver::SpatialScheme::QUICK,  // Schema alto ordine
            AdvectionDiffusionSolver::TimeScheme::RK4
        );
        
        // Condizione iniziale: sinusoidale
        solver.setInitialCondition(sinusoidale);
        
        // Condizioni al contorno periodiche
        solver.setBoundaryConditions(
            AdvectionDiffusionSolver::BoundaryCondition::PERIODIC
        );
        
        // Salva condizione iniziale
        solver.saveToFile("adv_diff_t0.dat");
        
        // Risolvi per diversi tempi
        double dt = 0.0005;
        
        solver.solve(dt, 0.2);
        solver.saveToFile("adv_diff_t02.dat");
        
        solver.solve(dt, 0.5);
        solver.saveToFile("adv_diff_t05.dat");
        
        solver.solve(dt, 1.0);
        solver.saveToFile("adv_diff_t10.dat");
        
        solver.printStats();
    }
    
    // ==========================================
    // ESEMPIO 4: Numero di Peclet alto (avvezione dominante)
    // ==========================================
    std::cout << "\n\n=== ESEMPIO 4: Alto numero di Peclet ===" << std::endl;
    {
        int nx = 401;
        double L = 1.0;
        double c = 1.0;     // Velocità avvezione alta
        double D = 0.001;   // Diffusione bassa -> Pe = cL/D = 1000
        
        std::cout << "Numero di Peclet Pe = " << (c * L / D) << std::endl;
        
        AdvectionDiffusionSolver solver(
            nx, L, c, D,
            AdvectionDiffusionSolver::SpatialScheme::UPWIND,
            AdvectionDiffusionSolver::TimeScheme::RK4
        );
        
        // Condizione iniziale: gaussiana
        solver.setInitialCondition(gaussiana);
        
        // Condizioni al contorno
        solver.setBoundaryConditions(
            AdvectionDiffusionSolver::BoundaryCondition::DIRICHLET,
            0.0, 0.0
        );
        
        solver.saveToFile("peclet_t0.dat");
        
        double dt = 0.0002;
        solver.solve(dt, 0.3);
        
        solver.saveToFile("peclet_t03.dat");
        solver.printStats();
    }
    */
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "          Simulation complete           " << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
