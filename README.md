# mpc_motion_planner

A joint space motion planner based on Model Predictive Control (MPC) to find the minimum time trajectory between the current state and target state (position and velocity), while respecting box constraints on joint position, velocity, acceleartion and torque.

The MPC library used is [polyMPC](https://gitlab.epfl.ch/listov/polympc) and we use [Ruckig](https://github.com/pantor/ruckig) as an initial guess for the solver.

## Installation

### Requires:

- [pinocchio](https://github.com/stack-of-tasks/pinocchio)

- [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (Version > 3.3)

### Installation steps

Clone the repository with the submodules, then compile the examples:

```
git clone --recurse-submodules git@github.com:AlbericDeLajarte/mpc_motion_planner.git 

cd mpc_motion_planner && mkdir build && cd build
cmake ..
make
```