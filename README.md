[![Build Status](https://travis-ci.org/skulumani/polyhedron_gravity.svg?branch=master)](https://travis-ci.org/skulumani/polyhedron_gravity)

## Polyhedron Gravity

This project defines the gravitational potential of a polyhedron asteroid. 
It implements the method derived by Werner and Scheeres.

* Exterior gravitation of a polyhedron derived and compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia 
* The gravitational potential of a homogeneous polyhedron or don't cut corners


## Installation

~~~
brew install cgal cmake libomp
~~~

## Requirements

1. CMake > 3.12

2. Eigen > 3.3.5

3. CGAL > 4.12-1

4. Boost > 1.68

~~~
git submodule update --init --recursive 
~~~

CMake doesn't search for matlab in libigl

Might need to edit 

~~~
extern/libigl/shared/cmake/libigl.cmake
~~~

## TODO

* Remove dependency on libigl
* Bundle CGAL requirements with app (similar to how igl works)
* Remove boost dependency (probably easy)
* Analytical verification of potential https://arxiv.org/pdf/1206.3857.pdf

~~~
/ The attracting shape is a cube of dimensions 1 x 1 x 1 m of density rho = 1e6 kg/m^3
	// The analytic potential and acceleration at (1,2,3) (m) in the shape model's barycentric frame is computed
	// Assumes that G = 6.67408e-11 m^3 / (kg * s ^2)

		// Queried point for pgm validation
		arma::vec p4 = {1, 2, 3};

		arma::vec acc_true = {
			-1.273782722739791e-06,
			-2.548008881415967e-06,
			-3.823026510474731e-06
		};
		double pot_true = 0.26726619638669064 * arma::datum::G * density;
~~~
