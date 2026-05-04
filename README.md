# MVP_CFD — Advection-Diffusion Solver

A numerical solver in C++ for the advection-diffusion equation, currently supporting 1D and 2D formulations. Designed as a modular, extensible framework with future support for RANS equations planned.

## Equation Solved

```
∂u/∂t + c·∂u/∂x = D·∂²u/∂x²
```

where:
- `u(x,t)` — unknown variable (concentration, temperature, etc.)
- `c` — advection speed [m/s]
- `D` — diffusion coefficient [m²/s]

## Features

### Spatial Schemes
| Scheme | Order | Notes |
|--------|-------|-------|
| Central Difference | 2nd | May produce spurious oscillations at high Pe |
| Upwind | 1st | Stable and monotone, recommended for advection-dominated flows |
| QUICK | 3rd | Higher accuracy, requires finer grids |

### Time Schemes
| Scheme | Order | Notes |
|--------|-------|-------|
| Explicit Euler | 1st | Simple, subject to CFL and diffusion stability limits |
| Runge-Kutta 4 | 4th | High accuracy, recommended for precision simulations |

### Boundary Conditions
- **Dirichlet** — fixed value at boundaries
- **Neumann** — fixed derivative at boundaries
- **Periodic** — u(0) = u(L)

### Initial Conditions
- **Gaussian** — smooth bell-shaped profile
- **Square wave** — sharp step profile
- **Sinusoidal** — smooth periodic oscillation

## Dimensionless Numbers

### Courant-Friedrichs-Lewy (CFL)
```
CFL = c·Δt/Δx ≤ 1   (explicit schemes)
```

### Diffusion Number
```
Diff = D·Δt/Δx² ≤ 0.5   (explicit schemes)
```

### Peclet Number
```
Pe = c·L/D
```
- Pe << 1 — diffusion dominated
- Pe >> 1 — advection dominated
- Pe ~ 1 — mixed regime

## Build

The project uses a Makefile with three build configurations:

```bash
make                  # Release build (optimized, -O3)
make BUILD=debug      # Debug build (-g, no optimization)
make BUILD=asan       # AddressSanitizer build (memory/UB checks)
```

Executables produced:
- `solver1d` — release
- `solver1d_debug` — debug
- `solver1d_asan` — asan

```bash
make run              # Build (release) and run with config.cfg
make clean            # Remove objects and executables
make cleanall         # Also remove output files
make help             # Show all available targets
```

## Usage

The quickest way to run is:
```bash
make run
```

This builds the project (if needed) and runs it automatically with config.cfg. To use a different config file, run the executable directly:

```bash
./solver1d <config_file>
# example:
./solver1d config.cfg
```

## Configuration File

The solver is fully configured via a `.cfg` text file. Lines starting with `#` are comments.

```ini
# Solver type
SOLVER = ADVECTION_DIFFUSION_1D      # or ADVECTION_DIFFUSION_2D

# Grid
GRID_POINTS    = 201
DOMAIN_LENGHT  = 1.0

# Physics
ADVECTION_SPEED = 1.0
DIFFUSION_COEFF = 0.01

# Time settings
TIME_STEP    = 0.001
INITIAL_TIME = 0.0
END_TIME     = 1.5

# Numerical schemes
SPATIAL_SCHEME = UPWIND              # CENTRAL | UPWIND | QUICK
TIME_SCHEME    = RK4                 # EULER_EXPLICIT | RK4

# Initial condition
INITIAL_CONDITIONS = GAUSSIAN        # GAUSSIAN | SQUARE_WAVE | SINUSOIDAL

# Boundary conditions
BOUNDARY_CONDITIONS = PERIODIC       # DIRICHLET | NEUMANN | PERIODIC
BC_LEFT  = 0.0
BC_RIGHT = 0.0

# Output
VERBOSE          = true
OUTPUT_FILE      = results/output.csv
MESH_FILE        = results/mesh.csv
OUTPUT_FREQUENCY = 10                # Save every N time steps
```

## Results and Visualization

Results are saved as `.csv` files in the `results/` directory:
- `mesh.csv` — grid coordinates
- `output.csv` — solution `u` at each saved time step (one row per step)

Visualization is handled by a Python script. Details to be added.

## Project Structure

```
MVP_CFD/
├── include/
│   └── solver/
│       ├── BaseSolver.h
│       ├── AdvectionDiffusionSolver1D.h
│       └── AdvectionDiffusionSolver2D.h
├── src/
│   ├── main.cpp
│   └── solver/
│       ├── BaseSolver.cpp
│       ├── AdvectionDiffusionSolver1D.cpp
│       └── AdvectionDiffusionSolver2D.cpp
├── config.cfg
├── Makefile
└── README.md
```

## Roadmap

- [x] 1D Advection-Diffusion solver
- [ ] 2D Advection-Diffusion solver *(in progress)*
- [ ] Implicit schemes (Crank-Nicolson)
- [ ] Non-uniform grids
- [ ] Variable coefficients c(x), D(x)
- [ ] Source terms
- [ ] RANS equations
- [ ] OpenMP/MPI parallelisation
- [ ] VTK output for Paraview

## References

1. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*
2. Versteeg & Malalasekera (2007). *An Introduction to Computational Fluid Dynamics*
3. Ferziger & Perić (2002). *Computational Methods for Fluid Dynamics*

## License

Free for educational and research use.