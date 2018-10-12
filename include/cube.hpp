/**
    Analytical potential and acceleration of a cube/sphere shapes

    @author Shankar Kulumani
    @version 9 October 2018
*/
#ifndef CUBE_H
#define CUBE_H

#include <Eigen/Dense>

class Cube {
    private:
        // member variables to hold the potential
        double mG = 6.673e-20; /**< Gravitational constant - km^3/kg/sec^2 */
        double msigma = 1 ; /**< Density - kg/km^3 */

        double mU;
        Eigen::Vector3d mU_grad;
        Eigen::Vector3d maxes;
        
    public:

        // constructor
        Cube( void );
        virtual ~Cube( void ) {};

        // other constructors
        Cube(const double& a_in, const double& b_in, const double& c_in);
        
        void potential(const Eigen::Ref<const Eigen::Vector3d>& state);
        void acceleration(const Eigen::Ref<const Eigen::Vector3d>& state);
        
        double get_potential( void ) { return mU; }
        Eigen::Vector3d get_acceleration( void ) { return mU_grad; }

        void set_grav_constant( const double& G_in ) { this->mG = G_in; }
        void set_sigma( const double& sigma_in ) { this->msigma = sigma_in; }
        void set_axes(const Eigen::Ref<const Eigen::Vector3d>& axes_in) { this->maxes = axes_in; }

        double get_grav_constant( void ) const { return mG; }
        double get_sigma( void ) const { return msigma; } 
        Eigen::Vector3d get_axes( void) const { return maxes; }
};
#endif
