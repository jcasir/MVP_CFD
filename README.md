# MVP_CFD — Advection-Diffusion Solver

A numerical solver in C++ for the advection-diffusion equation, currently supporting 1D and 2D formulations. Designed as a modular, extensible framework with future support for RANS equations planned.

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
- **Gaussian** — smooth bell-shaped profile 	(1D only)
- **Square wave** — sharp step profile			(1D & 2D)
- **Sinusoidal** — smooth periodic oscillations	(1D only)

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
SOLVER = ADVECTION_DIFFUSION_1D      #  ADVECTION_DIFFUSION_1D | ADVECTION_DIFFUSION_2D | BURGERS_2D

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

1D case:
Results are saved as `.csv` files in the `results/` directory:
- `mesh.csv` — grid coordinates
- `output.csv` — solution `u` at each saved time step (one row per step)

Visualization is handled by a Python script inside the `results/` directory: visualize.py

To visualize the results simply run:

```bash
pyhton3 visualize.py
```

2D case:
Results are saved as `.vtu` files in the `results/` directory. The `.pvd` file is used to visualize all the results saved in the `.vtu` files.

To visualize the results open the `.pvd` file using Paraview.


## Roadmap

- [x] 1D Advection-Diffusion solver
- [x] 2D Advection-Diffusion solver
- [x] 2D Burgers solver
- [ ] Implicit schemes (Crank-Nicolson)
- [ ] Non-uniform grids
- [ ] Source terms
- [ ] RANS equations
- [ ] OpenMP/MPI parallelisation
- [x] VTK output for Paraview

## References

1. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*
2. Lorena A. Barba *CFD Python: 12 steps to Navier-Stokes*
3. Politecnico di Milano. *Lecture notes from university courses.*
4. Tom-Robin Teschner. *cfd.university* 

## License

Free for educational and research use.