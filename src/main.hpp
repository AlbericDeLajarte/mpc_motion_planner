#include <iostream>
#include <fstream>
#include "robot_ocp.hpp"
#include "polympc_redef.hpp" 

#include "pandaWrapper.hpp"

#include <ruckig/ruckig.hpp>

class MotionPlanner{

    public:

        Eigen::Matrix<double, NDOF, 1> init_position{0.0, 0.0, 0.0, -1.5708, 0.0, 1.8675, 0.0};
        Eigen::Matrix<double, NDOF, 1> init_velocity{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
};