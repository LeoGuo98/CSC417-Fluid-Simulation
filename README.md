# CSC417-Fluid-Simulation
CSC417 Final Project, CPU implementation of SPH and FLIP method using Eigen and libigl.

## Install
Download project.zip and unzip it, you would see a folder named as CSC417-Fluid-Simulation. Move into that folder, and issue:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

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
Coming soon

## Report 
Coming soon
