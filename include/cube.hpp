/**
    Analytical potential and acceleration of a cube/sphere shapes

    @author Shankar Kulumani
    @version 9 October 2018
*/

#include <Eigen/Dense>

class Cube {
    private:
        // member variables to hold the potential
        double G = 6.673e-20; /**< Gravitational constant - km^3/kg/sec^2 */
        double sigma; /**< Density - kg/km^3 */

        double mU;
        Eigen::Vector3d mU_grad;

    public:

        // constructor
        Cube( void ) {};
        virtual ~Cube( void ) {};

        // other constructors
        Cube(const double& a_in, const double& b_in, const double& c_in);
        
        void potential(const Eigen::Ref<const Eigen::Vector3d>& state);
        void acceleration(const Eigen::Ref<const Eigen::Vector3d>& state);
        
        double get_potential( void ) { return mU; }
        Eigen::Vector3d get_acceleration( void ) { return mU_grad; }
        double get_grav_constant( void ) const { return G; }
        double get_sigma( void ) const { return sigma; } 
};
