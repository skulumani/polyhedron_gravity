#include "loader.hpp"
#include "mesh.hpp"

#include "gtest/gtest.h"

TEST(LoaderTest, OBJString) {
    std::string input_filename("./integration/cube.obj");
    std::shared_ptr<MeshData> mesh;
    mesh = Loader::load(input_filename);
    ASSERT_EQ(mesh->get_faces().rows(), 12);
    ASSERT_EQ(mesh->get_verts().rows(), 8);
}


