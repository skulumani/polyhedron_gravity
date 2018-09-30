#include <igl/dot_row.h>
#include <igl/readOBJ.h>
#include <igl/sort.h>
#include <igl/unique_rows.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>

#include "gtest/gtest.h"
#include <Eigen/Dense>

#include <iostream>

TEST(TestLibigl, RowDotProduct) {
    // define to matrices
    Eigen::MatrixXd a(2, 3), b(2, 3);
    a << 1, 1, 1,
        1, 0, 0;
    b << 1, 1, 1,
        0, 1, 0;

    Eigen::MatrixXd dot_prod = igl::dot_row(a, b);
    
    ASSERT_EQ(dot_prod(0), 3);
    ASSERT_EQ(dot_prod(1), 0);

}

TEST(TestLibigl, ReadOBJ) {

    Eigen::MatrixXd V, F;
    igl::readOBJ("./integration/cube.obj", V, F);
    
    ASSERT_EQ(V.rows(), 8);
    ASSERT_EQ(F.rows(), 12);
}

TEST(TestLibigl, SortRows) {
    Eigen::MatrixXi e1_sorted_true(12, 2);
    e1_sorted_true << 0, 6, 0, 2, 0, 3, 0, 1, 2, 7, 2, 3, 4, 6, 4, 7, 0, 4, 0, 5, 1, 5, 1, 7;

    // Define n x 2 array
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ("./integration/cube.obj", V, F);
    
    Eigen::MatrixXi Fa(F.col(0)), Fb(F.col(1)), Fc(F.col(2));
    Eigen::MatrixXi e1_vertex_map(F.rows(), 2);
    e1_vertex_map << Fb, Fa;

    Eigen::MatrixXi e1_sorted, e1_index;
    igl::sort(e1_vertex_map, 2, true, e1_sorted, e1_index);

    ASSERT_TRUE(e1_sorted.isApprox(e1_sorted_true));
}

TEST(TestLibigl, UniqueRow) {
    Eigen::MatrixXi edge_unique_true(2, 2);
    edge_unique_true << 1, 0, 1, 1;

    Eigen::MatrixXi edges(3, 2);
    edges << 1, 0, 1, 1, 1, 0;

    // now find unique values using igl
    Eigen::MatrixXi edge_unique, IA, IC;
    igl::unique_rows(edges, edge_unique, IA, IC);

    ASSERT_TRUE(edge_unique_true.isApprox(edge_unique));
}
