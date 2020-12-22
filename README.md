# CSC417-Fluid-Simulation
CSC417 Final Project, CPU implementation of SPH and FLIP method using Eigen and libigl.

## Install
Follow the <a href="https://github.com/dilevin/CSC417-a1-mass-spring-1d">installation</a> of the assignments.

## Execution
Once built, you can execute the project from inside the `build/`. Issue the following for simulation using SPH method:
```
./fluid SPH
```
Issue the following for simulation using FLIP method:
```
./fluid FLIP
```

For SPH method, press "a/A" to move one timestep forward, and press "r/R" to reset positions and velocities for all particles.
For FLIP method, press "a/A" to move one timestep forward for water simulation, press "b/B" to move one timestep forward for smoke simulation, and press "r/R" to reset positions and velocities for all particles.

## Video
https://youtu.be/tieIpwKc4KM

## Report 
See <a href="https://github.com/lihd1003/CSC417-Fluid-Simulation/blob/master/paper.pdf">paper.pdf</a>
