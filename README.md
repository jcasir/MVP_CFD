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
- `solver` — release
- `solver_debug` — debug
- `solver_asan` — asan

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
./solver <config_file>
# example:
./solver AdvDiff1D.cfg
```

## Configuration File

The solver is fully configured via a `.cfg` text file. Lines starting with `#` are comments.

```ini
# MVP CFD solver configuration file

# Solver type
SOLVER = BURGERS_2D 	 			# ADVECTION_DIFFUSION_1D | ADVECTION_DIFFUSION_2D | BURGERS_2D

# Grid options
GRID_POINTS_X = 201
GRID_POINTS_Y = 201
DOMAIN_LENGHT_X = 2.0
DOMAIN_LENGHT_Y = 2.0

# Simulation settings
DIFFUSION_COEFF = 0.0
TIME_STEP = 0.001
INITIAL_TIME = 0.0
END_TIME = 3.0

# Spatial and Time schemes
SPATIAL_SCHEME = QUICK   			# UPWIND | CENTRAL | QUICK
TIME_SCHEME = RK4					# EULER_EXPLICIT | RK4

# Initial Conditions Settings
INITIAL_CONDITIONS = SQUARE_WAVE

# Possible value for INITIAL_CONDITIONS:
#   GAUSSIAN    (1D only)  → requires: AMPLITUDE_IC, X0_IC, SIGMA_IC
#   SQUARE_WAVE (1D & 2D)  → requires: AMPLITUDE_IC, BASELINE_VALUE_IC, RANGE_START_IC, RANGE_END_IC (1D)
#                                       AMPLITUDE_IC, BASELINE_VALUE_IC, RANGE_START_X_IC, RANGE_END_X_IC,
#                                                     RANGE_START_Y_IC, RANGE_END_Y_IC (2D)
#   SINUSOIDAL  (1D only)  → requires: AMPLITUDE_IC
#   CONSTANT    (1D & 2D)  → requires: AMPLITUDE_IC

============================ Boundary Condition settings ==============================

# Advection-Diffusion solver (scalar field: e.g. temperature, concentration)
# Supports 1D and 2D configurations depending on grid setup
#
# Boundary condition is global (same type applied to all boundaries)
# Allowed values: DIRICHLET | NEUMANN | PERIODIC
#
# BC_* values define scalar boundary values (used only for DIRICHLET and NEUMANN)

BOUNDARY_CONDITIONS = PERIODIC  # DIRICHLET | NEUMANN | PERIODIC

BC_LEFT = 0.0
BC_RIGHT = 0.0
BC_TOP = 0.0
BC_BOTTOM = 0.0

=======================================================================================

# Burgers equation solver (2D vector field: U, V)
# Boundary conditions are defined per face (bottom, top, left, right)
#
# Supports Dirichlet (D) and Neumann (N) boundary conditions
#
# For periodic boundary conditions, set:
#   BOUNDARY_CONDITION = PERIODIC
# In this case, periodicity overrides individual boundary type settings.
#
# Example: lid-driven configuration (top lid moving in x-direction)

BOUNDARY_CONDITIONS = DIRICHLET  # DIRICHLET | NEUMANN | PERIODIC

BOUNDARY_CONDITIONS_VALUES_U = {0.0, 1.0, 0.0, 0.0}
BOUNDARY_CONDITIONS_VALUES_V = {0.0, 0.0, 0.0, 0.0}

=======================================================================================

# Boundary Conditions Settings for Incompressible Navier Stokes solver
# Layout: { bottom, top, left, right }
#
# BOUNDARY_CONDITIONS_VALUES_X defines the prescribed values for each boundary face.
# Example: BOUNDARY_CONDITIONS_VALUES_U = {0.0, 1.0, 0.0, 0.0}
#
# BOUNDARY_CONDITIONS_TYPES_X defines the type of boundary condition per face.
# Allowed values: Dirichlet (D), Neumann (N)
# Example: BOUNDARY_CONDITIONS_TYPES_U = {D, N, N, N}
#
# At least one Dirichlet boundary condition is required to ensure a well-posed problem
#
# For periodic boundary conditions, set:
#   BOUNDARY_CONDITION = PERIODIC
# In this case, periodicity overrides individual boundary type settings.
#
# Otherwise set:
#   BOUNDARY_CONDITION = NON_PERIODIC
# and use BOUNDARY_CONDITIONS_TYPES_X to define D/N conditions.
#
# Example: lid-driven cavity
# Top wall moves right (u = 1.0), all other walls are no-slip (u = 0.0)

BOUNDARY_CONDITIONS_VALUES_U = {0.0, 1.0, 0.0, 0.0}
BOUNDARY_CONDITIONS_VALUES_V = {0.0, 0.0, 0.0, 0.0}
BOUNDARY_CONDITIONS_VALUES_P = {0.0, 0.0, 0.0, 0.0}

BOUNDARY_CONDITIONS_TYPES_U = {D,D,D,D}
BOUNDARY_CONDITIONS_TYPES_V = {D,D,D,D}
BOUNDARY_CONDITIONS_TYPES_P = {N,D,N,N}

BOUNDARY_CONDITION = NON_PERIODIC

=======================================================================================

#Output flag
VERBOSE = true

#Results
OUTPUT_DIR = results/
OUTPUT_FILE = output
OUTPUT_FREQUENCY = 20

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
- [ ] RANS equations
- [ ] Implicit schemes
- [ ] Non-uniform grids
- [ ] OpenMP/MPI parallelisation
- [x] VTK output for Paraview

## References

1. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*
2. Lorena A. Barba *CFD Python: 12 steps to Navier-Stokes*
3. Politecnico di Milano. *Lecture notes from university courses.*
4. Tom-Robin Teschner. *cfd.university* 

## License

Free for educational and research use.