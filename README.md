# MolecularDynamicsSimulation
My masters project 2021

Hard sphere simulation
Simulates an ensemble of particles each with the same radius in a square box. The particles act as hard spheres and deflect off eachother and the the walls conserving both momentum and kinetic energy. The simulation is event driven meaning it works out the time of the next collision using the current positions and velocities of all the particles then fast forwards to this time.
An interpolation function works out the positions and velocities between these collisions. There are other functions for analysing results such as heatmaps, velocity distributions, diffusion graphs, particle paths, pair correlation and energy.
The simulation is unitless. To convert to SI units multiply distances by the r_0 and times by sqrt(m*r_0^2/K_B*T) where r_0 is the radius of the particles, m the mass of each particle, K_B the Boltzmann constant and T the temperature.

HSM animation
Simulates and then animates the hard sphere model simulation. Set to 2 dimensions as this is the one it will look best in.

Long Range potentials
Simulates an ensemble of particles using the Lennard-Jones potential U(r) = epsilon*((r_0/r)^12 - 2 (r_0/r)^6) where r is the distance between 2 particles and r_0 is the minimum potential energy. The simulation uses Verlet integration to get the next timestep.
Like the hard sphere simulation, there are functions for analysis such as heatmaps, velocity distributions, diffusion graphs, particle paths, pair correlation and energy.
The simulation is again unitless. To convert to SI units multiply distances by r_0 and times by sqrt(m*r_0^2/epsilon). For Argon gas r_0=3.8*10^-10 m, m=6.7*10^-26 kg and epsilon=1.6*10^-21 J. To work out the unitless temperature use K_B*T/epsilon.

LRP animation
Simulates and then animates the long range potential simulation. Set to 2 dimensions as this is the one it will look best in.