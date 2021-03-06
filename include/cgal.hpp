#ifndef CGAL_H
#define CGAL_H

#include "cgal_types.hpp"
#include "mesh.hpp"

#include <Eigen/Dense>

#include <memory>
#include <cmath>

/** @class MeshDistance

    @brief Find the distance from a point to a MeshData
    
    This function uses the KD Tree to find the closest vertex to a given point.

    @author Shankar Kulumani
    @version 20181022
*/
class MeshDistance {
    public:
        MeshDistance(std::shared_ptr<MeshData> mesh_in);
        
        /** @fn void update_mesh(std::shared_ptr<MeshData> mesh_in)
                
            Updates the internal member attributes with a new mesh

            @param MeshData input mesh to update object with
            @returns void

            @author Shankar Kulumani
            @version 20181023

        */
        void update_mesh(std::shared_ptr<MeshData> mesh_in);
        /** @fn int k_nearest_neighbor(const Eigen::Ref<const Eigen::Vector3d>& pt, const int &K)
                
            Compute the minimum distance and primitive to the mesh

            @param Eigen::Vector3d Test point in mesh body frame
            @returns int Success flag
            @return int K - index to closest primitive

            @author Shankar Kulumani
            @version 20181023
        */
        int k_nearest_neighbor(const Eigen::Ref<const Eigen::Vector3d>& pt, const int &K);

    private:
        std::shared_ptr<MeshData> mesh;
        Vertex_point_pmap vppmap;
};

/** @class RayCaster

    @brief Computes intersections of rays with a mesh
    
    @author Shankar Kulumani
    @version 20181024
*/
class RayCaster {
    public:
        RayCaster( void );
        virtual ~RayCaster( void ) {};

        RayCaster(std::shared_ptr<const MeshData> mesh_in);
        RayCaster(const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3> >& V_in,
                  const Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 3> >& F_in);
        
        void init_mesh(std::shared_ptr<const MeshData> mesh_in);
        void init_mesh(const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3> >& V_in,
                       const Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 3> >&  F_in);
        
        void accelerate( void );
        
        bool intersection(const Eigen::Ref<const Eigen::Vector3d>& psource,
                          const Eigen::Ref<const Eigen::Vector3d>& ptarget);

        // cast ray function
        Eigen::Matrix<double, 1, 3> castray(const Eigen::Ref<const Eigen::Vector3d>& psource,
                const Eigen::Ref<const Eigen::Vector3d>& ptarget);

        // cast many rays function
        Eigen::Matrix<double, Eigen::Dynamic, 3> castarray(const Eigen::Ref<const Eigen::Vector3d> &psource,
                                                           const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3> > &targets);
    
        // update the raycaster with a new mesh ptr
        void update_mesh(std::shared_ptr<const MeshData> mesh_in);
        void update_mesh(const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3> >& V_in,
                         const Eigen::Ref<const Eigen::Matrix<int, Eigen::Dynamic, 3> >& F_in);

        /**
            Compute the minimum distance to the mesh

            This uses the AABB tree to find the distance from a point to the 
            polyhedron.

            @param pt Eigen::Vector3d point defining the test point
            @returns Double distance type
        */
        double minimum_distance(const Eigen::Ref<const Eigen::Vector3d> &pt);

        void minimum_primitive(const Eigen::Ref<const Eigen::Vector3d> &pt);
    private:
        // needs the mesh to operate on
        std::shared_ptr<const MeshData> mesh;
        AABB_Tree tree; // holds the AABB tree for CGAL distance computations
};

#endif
