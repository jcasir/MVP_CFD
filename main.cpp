#include "AdvectionDiffusionSolver.h"
#include <iostream>
#include <cmath>


int main(int argc, char* argv[]) {

    // Catch the config file

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>\n"
        << "Config file missing\n";
        return 1;
    }

    std::string config_file {argv[1]};

    std::cout << "========================================" << std::endl;
    std::cout << "  SOLVER CFD 1D: Avvezione-Diffusione  " << std::endl;
    std::cout << "========================================" << std::endl;

    AdvectionDiffusionSolver solver(config_file);
    
    // // Condizione iniziale: gaussiana
    // solver.setInitialCondition(gaussiana);
    
    // // Condizioni al contorno periodiche
    // solver.setBoundaryConditions(
    //     AdvectionDiffusionSolver::BoundaryCondition::PERIODIC
    // );
    
    // Salva condizione iniziale
    solver.saveToFile("avvezione_t0.dat");
    
    // Risolvi
    solver.solve();
    
    // Salva risultato
    solver.saveToFile("avvezione_t05.dat");
    solver.printStats();

    
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
    std::cout << "  Tutti gli esempi completati!         " << std::endl;
    std::cout << "  I file .dat possono essere           " << std::endl;
    std::cout << "  visualizzati con gnuplot:            " << std::endl;
    std::cout << "  $ gnuplot -p -e \"plot 'file.dat'\"  " << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
