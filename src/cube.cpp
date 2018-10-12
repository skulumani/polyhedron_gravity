#include "cube.hpp"

#include <cmath>

Cube::Cube( void ) {
    this->msigma = 1;
    this->mU = 0;
    this->mU_grad = Eigen::Vector3d::Zero(3);

    this->maxes << 1, 1, 1; 
}

Cube::Cube(const double& a_in, const double& b_in, const double& c_in) {
    this->msigma = 1;
    this->mU = 0;
    this->mU_grad = Eigen::Vector3d::Zero(3);
    this->maxes << a_in, b_in, c_in;
}

void Cube::potential(const Eigen::Ref<const Eigen::Vector3d>& state) {
    double xmax = maxes[0] - state[0];
    double xmin = - maxes[0] - state[0];

    double ymax = maxes[1] - state[1];
    double ymin = - maxes[1] - state[1];

    double zmax = maxes[2] - state[2];
    double zmin = - maxes[2] - state[2];
    
    double r = sqrt(pow(state[0], 2) + pow(state[1], 2) + pow(state[2], 2));

    double U_max = 
          ymax * zmax * log(xmax + r) - pow(xmax, 2) / 2 * atan(ymax * zmax / (xmax * r)) 
        + xmax * zmax * log(ymax + r) - pow(ymax, 2) / 2 * atan(xmax * zmax / (ymax * r)) 
        + xmax * ymax * log(zmax + r) - pow(zmax, 2) / 2 * atan(xmax * ymax / (zmax * r));
    double U_min = 
          ymin * zmin * log(xmin + r) - pow(xmin, 2) / 2 * atan(ymin * zmin / (xmin * r)) 
        + xmin * zmin * log(ymin + r) - pow(ymin, 2) / 2 * atan(xmin * zmin / (ymin * r)) 
        + xmin * ymin * log(zmin + r) - pow(zmin, 2) / 2 * atan(xmin * ymin / (zmin * r));
    mU = - mG * msigma * (U_max - U_min);
}
