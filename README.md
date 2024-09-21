# Molecular Dynamics Simulation in FORTRAN90

## Overview
This project implements a molecular dynamics simulation using FORTRAN90. The program models the behavior of particles in a three-dimensional space, leveraging fundamental physics principles to simulate interactions and movements over time.

## Features
- **Particle Simulation:** Simulates a system of particles using Lennard-Jones potential.
- **Memory Management:** Efficient allocation of memory for particle properties.
- **Periodic Boundary Conditions:** Implements periodic boundaries to simulate an infinite system.
- **Data Output:** Generates trajectory data for visualization and analysis.

## Modules and Subroutines
- **variables:** Contains all necessary variables, including particle positions, velocities, and forces.
- **md:** Main program that orchestrates the simulation.
- **Memoria:** Allocates memory for particle properties.
- **LeeDatos:** Reads simulation parameters from an input file.
- **Posiciones:** Initializes particle positions in a cubic lattice.
- **Velocidades:** Assigns random velocities to particles.
- **Fuerza:** Calculates forces between particles based on the Lennard-Jones potential.
- **Repetir:** Updates particle positions and velocities, writing data to output files.
- **Pelicula:** Handles output for trajectory data and visualizations.
- **Histograma:** Generates histograms of velocity distributions.

## Requirements
- FORTRAN90 compiler (e.g., gfortran)

## Usage
1. Prepare the input file `run.dat` with the following parameters:
   - Number of particles (`nat`)
   - Density (`rho`)
   - Box dimensions (`lz`)
   - Time step (`dt`)
   - Interaction parameters (`sigma`, `eps`)
   - Number of steps (`npasos`)
   - Cut-off radius (`rcut`)

2. Compile the program:
   ```bash
   gfortran -o md md.f90
