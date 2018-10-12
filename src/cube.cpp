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
    
    double rmax = sqrt(pow(zmax,2.0) + pow(ymax,2.0) + pow(xmax, 2.0));
    double rmin = sqrt(pow(zmin, 2.0) + pow(ymin, 2.0) + pow(xmin, 2.0));

    double U_max = 
          ymax * zmax * log(xmax + rmax) - pow(xmax, 2.0) / 2 * atan(ymax * zmax / (xmax * rmax)) 
        + xmax * zmax * log(ymax + rmax) - pow(ymax, 2.0) / 2 * atan(xmax * zmax / (ymax * rmax)) 
        + xmax * ymax * log(zmax + rmax) - pow(zmax, 2.0) / 2 * atan(xmax * ymax / (zmax * rmax));
    double U_min = 
          ymin * zmin * log(xmin + rmin) - pow(xmin, 2.0) / 2 * atan(ymin * zmin / (xmin * rmin)) 
        + xmin * zmin * log(ymin + rmin) - pow(ymin, 2.0) / 2 * atan(xmin * zmin / (ymin * rmin)) 
        + xmin * ymin * log(zmin + rmin) - pow(zmin, 2.0) / 2 * atan(xmin * ymin / (zmin * rmin));
    mU = - mG * msigma * (U_max - U_min);
}
