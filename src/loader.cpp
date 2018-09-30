// Mesh Loader factory methods to create the mesh class from some input data

#include "loader.hpp"
#include "mesh.hpp"

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <memory>

// forward declare my functions
template<typename VectorType, typename IndexType>
int read_to_eigen(const std::string input_filename, Eigen::PlainObjectBase<VectorType> &V, Eigen::PlainObjectBase<IndexType> &F);

int read(const std::string input_filename, std::vector<std::vector<double>> &V, std::vector<std::vector<int>> &F);

int read(std::istream &input, std::vector<std::vector<double>> &V, std::vector<std::vector<int>> &F);

template<typename VectorType> 
void read_row(std::istringstream &ss, std::vector<VectorType> &vec); 

template<typename VectorType, typename Derived> 
int vector_array_to_eigen(std::vector<std::vector<VectorType> > &vector,
        Eigen::PlainObjectBase<Derived> &matrix);

// LOADER MEMBER FUNCTIONS
std::shared_ptr<MeshData> Loader::load(const std::string &input_filename) {
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
     
    // read from filename and create eigen arrays
    read_to_eigen(input_filename, vertices, faces);

    // instantiate the mesh and direct the pointer to it
    std::shared_ptr<MeshData> ptr = std::make_shared<MeshData>(vertices, faces);

    return ptr;
}

template <typename VectorType, typename IndexType> 
int read_to_eigen(const std::string input_filename, Eigen::PlainObjectBase<VectorType> &V,
        Eigen::PlainObjectBase<IndexType> &F) {
    // just call the stl vector version
    std::vector<std::vector<double> > V_vector;
    std::vector<std::vector<int> > F_vector;
    int read_flag = read(input_filename, V_vector, F_vector);
    V.resize(V_vector.size(), 3);
    F.resize(F_vector.size(), 3);

    if (read_flag == 0) {
        vector_array_to_eigen(V_vector,  V);
        vector_array_to_eigen(F_vector, F);
        return 0;
    } else {
        V = Eigen::MatrixXd::Zero(V_vector.size(), 3);
        F = Eigen::MatrixXi::Zero(F_vector.size(), 3);
        return 1;
    }

}

int read(const std::string input_filename, std::vector<std::vector<double>> &V, std::vector<std::vector<int>> &F) {
    std::ifstream input_stream;
    input_stream.open(input_filename);

    // check to make sure the file is opened properly
    if (!input_stream.fail()) {
        int read_flag = read(input_stream, V, F);
        return read_flag;
    } else {
        std::cout << "Error opening file filename" << std::endl;
        return 1;
    }
}

int read(std::istream& input, std::vector<std::vector<double>> &V, std::vector<std::vector<int>> &F) {

    if (input.fail()) {
        std::cout << "Error opening the file stream" << std::endl;
        return 1;
    }
    // store some strings for parsing the obj file
    std::string v("v"); // vertices
    std::string f("f"); // faces
    std::string octothorp("#"); // comments

    std::string line;

    while (std::getline(input, line)) {
        std::string row_type;
        std::istringstream row(line);

        row >> row_type;
        if (row_type == v) {
            std::vector<double> vertices;
            read_row(row, vertices);
            V.push_back(vertices);
            assert(vertices.size() == 3);
        } else if (row_type == f) {
            std::vector<int> indices;
            int v;
            while (row >> v) {
                indices.push_back(v - 1);
            }
            F.push_back(indices);
            assert(indices.size() == 3);
        }
    }
    /* input.clear(); */
    return 0;
}

template<typename VectorType> 
void read_row(std::istringstream &ss, std::vector<VectorType> &vector) {
    VectorType v;
    while (ss >> v) {
        vector.push_back(v);
    }
}

template<typename VectorType, typename Derived> 
int vector_array_to_eigen(std::vector<std::vector<VectorType> > &vector,
        Eigen::PlainObjectBase<Derived> &matrix) {
    // initialize a matrix to hold everything (assumes all are the same size
    int rows = vector.size();
    int cols = vector[0].size();
    matrix.resize(rows, cols);
    for (int ii = 0; ii < rows; ii++) {
        Eigen::Matrix<typename Derived::Scalar, 1, 3> v(vector[ii].data());
        matrix.row(ii) = v;
    }
    return 0;


}
// Explicit initialization
template int read_to_eigen<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);

template int vector_array_to_eigen<double, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);

template int vector_array_to_eigen<int, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
