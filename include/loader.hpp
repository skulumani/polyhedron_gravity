/**
    Load the data into the MeshData class

    @author Shankar Kulumani
    @version Date
*/
#ifndef LOADER_H
#define LOADER_H

#include <Eigen/Dense>
#include <memory>

// forward declaration
class MeshData;

class Loader {

    public:
        // factory methods to create a mesh
        static std::shared_ptr<MeshData> load(const std::string &input_filename);         
        static std::shared_ptr<MeshData> load(const std::istream &input_stream);
        static std::shared_ptr<MeshData> load(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
};
#endif
