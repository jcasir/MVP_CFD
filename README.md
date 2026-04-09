# Risolutore CFD 1D: Equazione di Avvezione-Diffusione

Un solver numerico in C++ per l'equazione di avvezione-diffusione monodimensionale.

## Equazione Risolta

```
∂u/∂t + c·∂u/∂x = D·∂²u/∂x²
```

dove:
- `u(x,t)` è la variabile incognita (concentrazione, temperatura, ecc.)
- `c` è la velocità di avvezione [m/s]
- `D` è il coefficiente di diffusione [m²/s]
- `x` è la coordinata spaziale
- `t` è il tempo

## Caratteristiche

### Schemi Spaziali
1. **Central Difference** (2° ordine)
   - Accurato ma può presentare oscillazioni spurie
   - Ottimo per diffusione dominante (Pe < 2)

2. **Upwind** (1° ordine)
   - Stabile e monotono
   - Più diffusivo numericamente
   - Consigliato per avvezione dominante

3. **QUICK** (3° ordine)
   - Quadratic Upstream Interpolation
   - Maggiore accuratezza
   - Necessita di griglia più fine

### Schemi Temporali
1. **Eulero Esplicito** (1° ordine)
   - Semplice e veloce
   - Limitazioni di stabilità: CFL < 1, numero diffusione < 0.5

2. **Runge-Kutta 4** (4° ordine)
   - Alta accuratezza temporale
   - Maggiore stabilità
   - Consigliato per simulazioni di precisione

### Condizioni al Contorno
1. **Dirichlet**: Valore fisso ai bordi
2. **Neumann**: Derivata fissa ai bordi
3. **Periodiche**: u(0) = u(L)

## Numeri Adimensionali

### Numero di Courant-Friedrichs-Lewy (CFL)
```
CFL = c·Δt/Δx
```
**Condizione di stabilità** (Eulero esplicito): CFL ≤ 1

### Numero di Diffusione
```
Diff = D·Δt/Δx²
```
**Condizione di stabilità** (Eulero esplicito): Diff ≤ 0.5

### Numero di Peclet
```
Pe = c·L/D
```
- Pe << 1: Diffusione dominante
- Pe >> 1: Avvezione dominante
- Pe ~ 1: Regime misto

## Compilazione

```bash
# Compilazione standard (ottimizzata)
make

# Compilazione debug
make debug

# Pulizia
make clean
```

## Utilizzo

### Esempio Base

```cpp
#include "AdvectionDiffusionSolver.h"

// Parametri del problema
int nx = 201;           // Numero di punti griglia
double L = 1.0;         // Lunghezza dominio [m]
double c = 1.0;         // Velocità avvezione [m/s]
double D = 0.01;        // Coefficiente diffusione [m²/s]

// Crea solver
AdvectionDiffusionSolver solver(
    nx, L, c, D,
    AdvectionDiffusionSolver::SpatialScheme::CENTRAL,
    AdvectionDiffusionSolver::TimeScheme::RK4
);

// Imposta condizione iniziale
solver.setInitialCondition([](double x) {
    return exp(-pow(x - 0.5, 2) / 0.01);
});

// Imposta condizioni al contorno
solver.setBoundaryConditions(
    AdvectionDiffusionSolver::BoundaryCondition::DIRICHLET,
    0.0,  // Valore a x=0
    0.0   // Valore a x=L
);

// Risolvi
double dt = 0.001;      // Passo temporale [s]
double tEnd = 1.0;      // Tempo finale [s]
solver.solve(dt, tEnd);

// Salva risultati
solver.saveToFile("output.dat");
solver.printStats();
```

### Scelta dei Parametri Numerici

1. **Passo spaziale (Δx)**:
   - Risolvi le strutture più piccole: Δx < λ/10
   - Numero di Peclet di griglia: Pe_grid = c·Δx/D < 2 (per schema centrale)

2. **Passo temporale (Δt)**:
   - CFL ≤ 1 per stabilità
   - Diff ≤ 0.5 per stabilità
   - Δt < min(Δx/c, Δx²/2D)

3. **Schema numerico**:
   - Pe > 2: usa Upwind o QUICK
   - Pe < 2: usa Central Difference
   - Alta precisione: usa RK4

## Esempi Inclusi

Il file `main.cpp` contiene 4 esempi dimostrativi:

1. **Avvezione pura** (D=0): Trasporto di un profilo gaussiano
2. **Diffusione pura** (c=0): Smussamento di un'onda quadra
3. **Avvezione-Diffusione**: Evoluzione di un profilo sinusoidale
4. **Alto Peclet**: Regime dominato dall'avvezione

Esegui con:
```bash
make run
```

## Visualizzazione Risultati

I risultati sono salvati in formato ASCII (file `.dat`) con colonne `x u(x,t)`.

### Gnuplot
```bash
gnuplot -p -e "plot 'avvezione_t0.dat' w l, 'avvezione_t05.dat' w l"
```

### Python
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('output.dat')
plt.plot(data[:,0], data[:,1])
plt.xlabel('x')
plt.ylabel('u')
plt.show()
```

## Applicazioni

- Trasporto di inquinanti in fiume
- Diffusione termica con convezione
- Trasporto di specie chimiche
- Propagazione di segnali
- Modelli di traffico

## Sviluppi Futuri

- [ ] Schemi impliciti (Crank-Nicolson)
- [ ] Griglia non uniforme
- [ ] Coefficienti variabili c(x), D(x)
- [ ] Termini sorgente
- [ ] Equazioni 2D/3D
- [ ] Parallelizzazione OpenMP/MPI
- [ ] Output VTK per Paraview

## Riferimenti

1. LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*
2. Versteeg & Malalasekera (2007). *An Introduction to Computational Fluid Dynamics*
3. Ferziger & Perić (2002). *Computational Methods for Fluid Dynamics*

## Licenza

Codice libero per uso educativo e di ricerca.

## Autore

Progetto didattico per apprendimento CFD
