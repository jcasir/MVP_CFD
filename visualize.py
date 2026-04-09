#!/usr/bin/env python3
"""
Script per visualizzare i risultati del solver CFD 1D
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

def plot_single_file(filename, **kwargs):
    """Plotta un singolo file di dati"""
    if not os.path.exists(filename):
        print(f"Errore: file {filename} non trovato")
        return None
    
    data = np.loadtxt(filename, comments='#')
    x = data[:, 0]
    u = data[:, 1]
    
    label = kwargs.pop('label', filename)
    plt.plot(x, u, linewidth=2, label=label, **kwargs)
    return x, u

def plot_comparison(files, title="Confronto Soluzioni"):
    """Plotta più file per confronto"""
    plt.figure(figsize=(10, 6))
    
    for filename in files:
        plot_single_file(filename)
    
    plt.xlabel('x', fontsize=12)
    plt.ylabel('u(x,t)', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_evolution(prefix, times, title="Evoluzione Temporale"):
    """Plotta l'evoluzione temporale di una soluzione"""
    plt.figure(figsize=(12, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    
    for i, t in enumerate(times):
        filename = f"{prefix}_t{t:02d}.dat"
        if os.path.exists(filename):
            plot_single_file(filename, color=colors[i], 
                           label=f't = {t/10:.1f}')
    
    plt.xlabel('x', fontsize=12)
    plt.ylabel('u(x,t)', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

def plot_all_examples():
    """Plotta tutti gli esempi del main"""
    
    # Esempio 1: Avvezione
    if os.path.exists('avvezione_t0.dat') and os.path.exists('avvezione_t05.dat'):
        plt.figure(figsize=(10, 6))
        plot_single_file('avvezione_t0.dat', label='t = 0.0', linestyle='--')
        plot_single_file('avvezione_t05.dat', label='t = 0.5')
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x,t)', fontsize=12)
        plt.title('Esempio 1: Avvezione Pura', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
    
    # Esempio 2: Diffusione
    if os.path.exists('diffusione_t0.dat') and os.path.exists('diffusione_t01.dat'):
        plt.figure(figsize=(10, 6))
        plot_single_file('diffusione_t0.dat', label='t = 0.0', linestyle='--')
        plot_single_file('diffusione_t01.dat', label='t = 0.1')
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x,t)', fontsize=12)
        plt.title('Esempio 2: Diffusione Pura', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
    
    # Esempio 3: Avvezione-Diffusione
    files = ['adv_diff_t0.dat', 'adv_diff_t02.dat', 
             'adv_diff_t05.dat', 'adv_diff_t10.dat']
    if all(os.path.exists(f) for f in files):
        plt.figure(figsize=(12, 7))
        times = [0.0, 0.2, 0.5, 1.0]
        colors = plt.cm.viridis(np.linspace(0, 1, len(files)))
        for i, (f, t) in enumerate(zip(files, times)):
            plot_single_file(f, color=colors[i], label=f't = {t:.1f}')
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x,t)', fontsize=12)
        plt.title('Esempio 3: Avvezione-Diffusione', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
    
    # Esempio 4: Alto Peclet
    if os.path.exists('peclet_t0.dat') and os.path.exists('peclet_t03.dat'):
        plt.figure(figsize=(10, 6))
        plot_single_file('peclet_t0.dat', label='t = 0.0', linestyle='--')
        plot_single_file('peclet_t03.dat', label='t = 0.3')
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x,t)', fontsize=12)
        plt.title('Esempio 4: Alto Numero di Peclet', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
    
    plt.show()

def main():
    """Funzione principale"""
    if len(sys.argv) == 1:
        print("Visualizzazione di tutti gli esempi...")
        plot_all_examples()
    elif len(sys.argv) == 2:
        print(f"Visualizzazione di {sys.argv[1]}...")
        plt.figure(figsize=(10, 6))
        plot_single_file(sys.argv[1])
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x,t)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()
    else:
        print(f"Confronto di {len(sys.argv)-1} file...")
        plot_comparison(sys.argv[1:])

if __name__ == "__main__":
    print("""
    =========================================
      Visualizzatore Risultati Solver 1D
    =========================================
    
    Uso:
      python3 visualize.py                  # Mostra tutti gli esempi
      python3 visualize.py file.dat         # Mostra un file
      python3 visualize.py file1.dat file2.dat ...  # Confronta più file
    """)
    
    main()
