# MolecularDynamicsSimulation
My masters project 2021

Hard sphere simulation
Simulates an ensemble of particles each with the same radius in a square box. The particles act as hard spheres and deflect off eachother and the the walls conserving both momentum and kinetic energy. The simulation is event driven meaning it works out the time of the next collision using the current positions and velocities of all the particles then fast forwards to this time.
An interpolation function works out the positions and velocities between these collisions. There are other functions for analysing results such as heatmaps, velocity distributions, diffusion graphs, particle paths, pair correlation and energy.
The simulation is unitless. To convert to SI units multiply distances by the r_0 and times by sqrt(mr_0^2/K_B*T) where r_0 is the radius of the particles, m the mass of each particle, K_B the Boltzmann constant and T the temperature.

HSM animation

Long Range potentials

LRP animation