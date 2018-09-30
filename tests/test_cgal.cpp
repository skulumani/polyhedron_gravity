#include "cgal.hpp"
#include "mesh.hpp"
#include "loader.hpp"

#include <gtest/gtest.h>

#include <iostream>

// The fixture for testing class Foo.
class TestRayCaster: public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  TestRayCaster() {
    // You can do set-up work for each test here.
    // initialize the mesh
    Ve_true << -0.5,   -0.5, -0.5,
            -0.5, -0.5, 0.5,
            -0.5, 0.5,  -0.5,
            -0.5, 0.5,  0.5,
            0.5,  -0.5, -0.5,
            0.5,  -0.5, 0.5,
            0.5,  0.5,  -0.5,
            0.5,  0.5,  0.5;

    Fe_true << 1, 7, 5,
            1, 3, 7,
            1, 4, 3,
            1, 2, 4,
            3, 8, 7,
            3, 4, 8,
            5, 7, 8,
            5, 8, 6,
            1, 5, 6,
            1, 6, 2,
            2, 6, 8,
            2, 8, 4;
    Fe_true = Fe_true.array() - 1;
  }

  virtual ~TestRayCaster() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
    const std::string input_file = "./integration/cube.obj";
    const std::string itokawa_file = "./data/shape_model/ITOKAWA/itokawa_low.obj";
    Eigen::Matrix<double, 8, 3> Ve_true;
    Eigen::Matrix<int, 12, 3> Fe_true;
};

TEST_F(TestRayCaster, NoIntersectionCube) {
    std::shared_ptr<MeshData> mesh;

    mesh = Loader::load(input_file);
    RayCaster caster(mesh);

    Eigen::Vector3d psource(3), ptarget(3);
    psource << 2, 0, 0;
    ptarget << 5, 0, 0;
    
    Eigen::RowVector3d intersection = caster.castray(psource, ptarget);

    ASSERT_TRUE(intersection.isApprox(Eigen::RowVector3d::Zero()));
}

TEST_F(TestRayCaster, ItokawaIntersection) {
    std::shared_ptr<MeshData> mesh;

    mesh = Loader::load(itokawa_file);
    RayCaster caster(mesh);

    Eigen::RowVector3d psource(3), ptarget(3);
    psource << 5, 0, 0;
    ptarget << 0, 0, 0;

    Eigen::RowVector3d intersection = caster.castray(psource, ptarget);
    Eigen::RowVector3d intersection_true(3);
    intersection_true << 0.29, 0, 0;

    ASSERT_NEAR(intersection(0), intersection_true(0.29), 1e-2);

}

TEST_F(TestRayCaster, InitMeshIntersection) {
    RayCaster caster;
    std::shared_ptr<MeshData> mesh = Loader::load("./data/shape_model/ITOKAWA/itokawa_low.obj");

    caster.init_mesh(mesh);
    Eigen::RowVector3d psource(3), ptarget(3);
    psource << 5, 0, 0;
    ptarget << -5, 0, 0;

    ASSERT_TRUE(caster.intersection(psource, ptarget));  
}

TEST_F(TestRayCaster, InitMeshNoIntersection) {
    RayCaster caster;
    std::shared_ptr<MeshData> mesh = Loader::load("./data/shape_model/ITOKAWA/itokawa_low.obj");

    caster.init_mesh(mesh);
    Eigen::RowVector3d psource(3), ptarget(3);
    psource << 5, 0, 0;
    ptarget << 10, 0, 0;

    ASSERT_FALSE(caster.intersection(psource, ptarget));  
}
