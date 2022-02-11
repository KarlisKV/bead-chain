# magnetic-bead-chain-dipole-dipole-interaction-energy
Magnetic bead chain simulation in a gravitational field. Shows a graph and calculates the
potential energy of magnetic beads using an optimization algorithm (Scipy's Nelder-Mead)
that finds the positions of beads given initial directions for each bead's magnetic moment. 
Takes in to account gravity, which can be turned off by setting it to 0, in that case will 
return a linear chain. 

Runtime:
2 beads = 10s
10 beads = 5min
20 beads = 40min
In other words, time complexity of O(N^2)

Runs on Python 3
