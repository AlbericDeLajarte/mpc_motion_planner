# mpc_motion_planner

A joint space motion planner based on Model Predictive Control (MPC) to find the minimum time trajectory between the current state and target state (position and velocity), while respecting box constraints on joint position, velocity, acceleartion and torque.

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

## Comments:
- Turn on release mode for ~50 time increase in speed