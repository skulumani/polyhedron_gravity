#include "cube.hpp"

#include <Eigen/Dense>

#include <iostream>

int main() {

    Eigen::Vector3d axes(2, 2, 2);
    Eigen::Vector3d state(1000, 0, 0);

    Cube cube_analytical;
    
    cube_analytical.set_grav_constant(1);
    cube_analytical.set_sigma(1);
    cube_analytical.set_axes(axes);
    
    cube_analytical.potential(state);

    std::cout << "G: " << cube_analytical.get_grav_constant() << std::endl;
    std::cout << "Sigma: " << cube_analytical.get_sigma() << std::endl;
    std::cout << "Axes: " << cube_analytical.get_axes().transpose() << std::endl;
    std::cout << "Potential: " << cube_analytical.get_potential() << std::endl;
    std::cout << "True Potential: " << 8 / state[0] << std::endl;
    std::cout << "Acceleration: " << cube_analytical.get_acceleration().transpose() << std::endl;

    return 0;
}
