#include "potential.hpp"
#include "loader.hpp"
#include "mesh.hpp"

#include "gtest/gtest.h"

#include <iostream>

TEST(TestAsteroid, CubeGravity) {
    std::shared_ptr<MeshData> mesh_data = Loader::load("./integration/cube.obj");
    Eigen::Matrix<double, 1, 3> state;
    state << 1, 2, 3;

    Asteroid ast("cube", mesh_data);
    ast.polyhedron_potential(state);

    const double U_true = 1.7834673284883827e-08;
    Eigen::Matrix<double, 3, 1> U_grad_true;
    U_grad_true <<-1.27357660e-09, -2.54759656e-09, -3.82240787e-09;
    Eigen::Matrix<double, 3, 3> U_grad_mat_true;
    U_grad_mat_true << -1.00116623e-09,  5.45412364e-10,  8.18746722e-10,
 5.45412364e-10, -1.82972253e-10,  1.63835424e-09,
 8.18746722e-10,  1.63835424e-09,  1.18413848e-09;
    const double Ulaplace_true = 9.260647804154587e-25;

    ASSERT_NEAR(ast.get_potential(), U_true, 1e-7);
    ASSERT_TRUE(ast.get_acceleration().isApprox(U_grad_true, 1e-7));
    ASSERT_TRUE(ast.get_gradient_mat().isApprox(U_grad_mat_true, 1e-7));
    ASSERT_NEAR(ast.get_laplace(), Ulaplace_true, 1e-7);
}

TEST(TestAsteroid, TetrahedronGravity) {
    std::shared_ptr<MeshData> mesh_data = Loader::load("./integration/tetrahedron.obj");
    Eigen::Matrix<double, 1, 3> state;
    state << 1, 2, 3;

    Asteroid ast("cube", mesh_data);
    ast.polyhedron_potential(state);

    const double U_true = -9.15229252142212e-09;
    Eigen::Matrix<double, 3, 1> U_grad_true;
    U_grad_true << 6.55438157e-10, 1.31045978e-09, 1.95837246e-09;
    Eigen::Matrix<double, 3, 3> U_grad_mat_true;
    U_grad_mat_true <<  5.12578025e-10, -2.82609296e-10, -4.21082987e-10,
        -2.82609296e-10,  8.87352037e-11, -8.41551959e-10,
        -4.21082987e-10, -8.41551959e-10, -6.01313229e-10;
 
    const double Ulaplace_true = 1.1575809755193234e-25;

    ASSERT_NEAR(ast.get_potential(), U_true, 1e-7);
    ASSERT_TRUE(ast.get_acceleration().isApprox(U_grad_true, 1e-7));
    ASSERT_TRUE(ast.get_gradient_mat().isApprox(U_grad_mat_true, 1e-7));
    ASSERT_NEAR(ast.get_laplace(), Ulaplace_true, 1e-7);
    
}

TEST(TestAsteroid, AnalyticalCube) {
    std::shared_ptr<MeshData> mesh = Loader::load("./data/shape_model/CUBE/cube.obj");
    Asteroid ast("cube", mesh);
    ast.set_grav_constant(6.67408e-11 * (1.0 / 1000.0) * (1.0 / 1000.0) * (1.0 / 1000.0));
    ast.set_sigma(1e6 * (1000.0 ) * (1000.0) * (1000.0));
    Eigen::Matrix<double, 1, 3> state;
    state << 1, 2, 3;
    
    Eigen::Vector3d U_grad_true;
    U_grad_true << -1.273782722739791e-06, -2.548008881415967e-06,
                -3.823026510474731e-06;

    ast.polyhedron_potential(state);
    
    ASSERT_NEAR(ast.get_potential(),0.26726619638669064 * ast.get_grav_constant() * ast.get_sigma(), 1e-10);
    ASSERT_TRUE(ast.get_acceleration().isApprox(U_grad_true, 1e-10));
}

TEST(TestAsteroid, CastaliaGravity) {
    std::shared_ptr<MeshData> mesh_data = Loader::load("./data/shape_model/CASTALIA/castalia.obj");
    Eigen::Matrix<double, 1, 3> state;
    state << 1, 2, 3;

    Asteroid ast("castalia", mesh_data);
    ast.polyhedron_potential(state);

    const double U_true = 1.7834673284883827e-08;
    Eigen::Matrix<double, 3, 1> U_grad_true;
    U_grad_true << -1.721526e-09, -3.5357462e-09, -5.31929459e-09;
    Eigen::Matrix<double, 3, 3> U_grad_mat_true;
    U_grad_mat_true << -1.38277686e-09, 7.19253451e-10, 1.08545592e-09,
                    7.19253451e-10, -2.63370870e-10, 2.26872578e-09,
                    1.08545592e-09, 2.26872578e-09, 1.64614773e-09;
    const double Ulaplace_true = -1.44639742e-23;

    ASSERT_NEAR(ast.get_potential(), U_true, 1e-7);
    ASSERT_TRUE(ast.get_acceleration().isApprox(U_grad_true, 1e-7));
    ASSERT_TRUE(ast.get_gradient_mat().isApprox(U_grad_mat_true, 1e-7));
    ASSERT_NEAR(ast.get_laplace(), Ulaplace_true, 1e-7);

}

TEST(TestAsteroid, CastaliaSurfaceSlope) {
    std::shared_ptr<MeshData> mesh_data = Loader::load("./data/shape_model/CASTALIA/castalia.obj");
    Asteroid ast("castalia", mesh_data);

    Eigen::VectorXd face_slope = ast.surface_slope( );

    EXPECT_NEAR(face_slope.maxCoeff(), 0.669671, 1e-3);
    EXPECT_NEAR(face_slope.minCoeff(), 0.00451144, 1e-3);
}
