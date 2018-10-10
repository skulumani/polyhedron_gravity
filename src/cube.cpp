#include "cube.hpp"

Cube::Cube( void ) {
    this->sigma = 1;
    this->mU = 0;
    this->mU_grad = Eigen::Vector3d::Zeros(3);

    this->axes << 1, 1, 1; 
}

Cube::Cube(const double& a_in, const double& b_in, const double& c_in) {
    this->sigma = 1;
    this->mU = 0;
    this->mU_grad = Eigen::Vector3d::Zeros(3);
    this->axes << a_in, b_in, c_in;
}


