#include <iostream>
#include <fstream>
#include "robot_ocp.hpp"
#include "polympc_redef.hpp" 

#include "pandaWrapper.hpp"

#include <ruckig/ruckig.hpp>

class MotionPlanner{

    public:

        std::array<double, NDOF> init_position;
        std::array<double, NDOF> init_velocity;
    
        MotionPlanner(){

            init_position = {0.0, 0.0, 0.0, -1.5708, 0.0, 1.8675, 0.0};
            init_velocity = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

            // target_position = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
            // target_velocity = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};


        }

};