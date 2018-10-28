#ifndef MESH_H
#define MESH_H

#include "cgal_types.hpp"

#include <Eigen/Dense>

#include <vector>

/** @class MeshData

    @brief Class to hold a polyhedron mesh.
    
    This is a wrapper of the Surface_mesh from CGAL and is the foundation used for the asteroid polyhedron potential model.

    @author Shankar Kulumani
    @version 28 October 2018
*/
class MeshData {
    public:
        // constructor
        MeshData( void ) {}
        virtual ~MeshData( void ) {}

        // other constructors
        MeshData(const Eigen::Ref<const Eigen::MatrixXd> &V,
                const Eigen::Ref<const Eigen::MatrixXi> &F);

        /** @fn void update_mesh(const Eigen::Ref<const Eigen::MatrxiXd>& V, const Eigen::Ref<const Eigen::MatrixXi>& F)
                
            Completely update the mesh with new vertices and faces. This
            will also completely regenerate the mesh properties

            @param V Eigen matrix with the vertices
            @param F Eigen matrix with teh faces
            @returns void

            @author Shankar Kulumani
            @version 9 June 2018
        */ 
        void update_mesh(const Eigen::Ref<const Eigen::MatrixXd> &V, 
                const Eigen::Ref<const Eigen::MatrixXi> &F);

    
        // TODO Move these two functions to the Asteroid class
        /** @fn build_edge_factor(const Eigen::Ref<const Eigen::Vector3d>& pos)
                
            Build the per edge factor L for the polyhedron potential model

            @param pos Position of test point
            @returns bool Success if true

            @author Shankar Kulumani
            @version 28 October 2018
        */
        bool build_edge_factor(const Eigen::Ref<const Eigen::Vector3d>& pos);

        /** @fn build_face_factor(const Eigen::Ref<const Eigen::Vctor3d>& pos)
                
            Build the per face factor w for the polyhedron potential model

            @param pos Position of test point
            @returns bool success if true

            @author Shankar Kulumani
            @version 28 October 2018
        */
        bool build_face_factor(const Eigen::Ref<const Eigen::Vector3d>& pos);

        /** @fn get_verts(void)
                
            Return the vertices of the mesh as an Eigen array

            @returns vertices Eigen arrays for the mesh nx3

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Eigen::Matrix<double, Eigen::Dynamic, 3> get_verts( void ) const;

        /** @fn Eigen::Matrix<int, Eigen::Dynamic, 3> get_faces( void) const
                
            Return the faces of the mesh as an eigen array

            @returns faces Eigen array for the faces of the mesh nx3

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Eigen::Matrix<int, Eigen::Dynamic, 3> get_faces( void ) const;
        
        /** @fn std::size_t number_of_vertices( void ) const
                
            Return the number of vertices of the mesh

            @returns num_vert Return the number of vertices

            @author Shankar Kulumani
            @version 28 October 2018
        */
        std::size_t number_of_vertices( void ) const;

        /** @fn std::size_t number_of_edges( void ) const
                
            Return number of edges of the mesh

            @returns num_edge Number of edges

            @author Shankar Kulumani
            @version 28 October 2018
        */
        std::size_t number_of_edges( void ) const;

        /** @fn std::size_t number_of_faces( void ) const
                
            Return number of faces of the mesh

            @returns num_face Number of faces

            @author Shankar Kulumani
            @version 28 October 2018
        */
        std::size_t number_of_faces( void ) const;
        
        /** @fn std::size_t number_of_halfedges( void ) const
                
            Return number of halfedges of the mesh

            @returns num_half Number of halfedges

            @author Shankar Kulumani
            @version 28 October 2018
        */
        std::size_t number_of_halfedges( void ) const;
        
        // Range types for iteration
        /** @fn Mesh::Vertex_range vertices( void ) const
                
            Return vertex range for iteration

            @returns Vertex_range Range of vertices of the mesh

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Mesh::Vertex_range vertices( void ) const;

        /** @fn Mesh::Face_range faces( void ) const
                
            Return face range for iteration

            @returns Face_range Range of faces of the mesh

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Mesh::Face_range faces( void ) const;

        /** @fn Mesh::Edge_range edges( void ) const
                
            Return edge range for iteration

            @returns Edge_range Range of edges of the mesh

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Mesh::Edge_range edges( void ) const;

        /** @fn Mesh::Halfedge_range halfedges( void ) const
                
            Return halfedge range for iteration

            @returns Halfedge_range Range of halfedges of the mesh

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Mesh::Halfedge_range halfedges( void ) const;
        
        /** @fn bool refine_faces(const std::vector<Face_index>& face_vec,
         *              std::vector<Face_index>& new_faces,
         *              std::vector<Vertex_index>& new_vertices,
         *              const int& density = 4.0)
                
            Given a set of faces, this will refine them by adding new faces/
            vertices within them. The new faces and vertices are returned.
            The surface mesh member variable is updated.

            This tends to make a very craggly shape. Better to use remesh
            instead
            
            @param face_vec Vector of faces to refine
            @param density integer defining the density control factor
            @returns new_faces Vector of new faces
            @returns new_vertices Vector of new vertices

            @author Shankar Kulumani
            @version 10 June 2018
        */
        bool refine_faces(const std::vector<Face_index>& face_vec, 
                std::vector<Face_index>& new_faces,
                std::vector<Vertex_index>& new_vertices,
                const int& density = 4.0);
        
        /** @fn remesh_faces(const std::vector<Face_index>& face_vec, const double& target_edge_length, const int& number_of_iterations)
                
            Perform isotropic remeshing of the desired faces of the mesh.

            @param face_vec Vector of faces to remesh
            @param target_edge_length Target length of the new edges
            @param number_of_iterations Iterations for isotropic_remeshing
            @returns bool True if success

            @author Shankar Kulumani
            @version 28 October 2018
        */
        bool remesh_faces(const std::vector<Face_index>& face_vec,
                const double& target_edge_length,
                const int& number_of_iterations=3);

        /** @fn std::vector<Face_index> faces_in_fov(
         *      const Eigen::Ref<const Eigen::Vector3d>& pos,
         *      const double& max_fov=0.52)
                
            Find the faces that are within a FOV of the current position

            @param pos Position of spacecraft in the asteroid fixed frame
            @returns face_vec Vector of face indices

            @author Shankar Kulumani
            @version 11 June 2018
        */
        std::vector<Face_index> faces_in_fov(
                const Eigen::Ref<const Eigen::Vector3d>& pos,
                const double& max_fov=0.52);

        /** @fn std::vector<Face_index> vertices_in_fov(
         *      const Eigen::Ref<const Eigen::Vector3d>& pos,
         *      const double& max_fov=0.52)
                
            Find the vertices that are within a FOV of the current position

            @param pos Position of spacecraft in the asteroid fixed frame
            @returns vertex_vec Vector of vertex indices

            @author Shankar Kulumani
            @version 28 October 2018
        */
        std::vector<Vertex_index> vertices_in_fov(
                const Eigen::Ref<const Eigen::Vector3d>& pos,
                const double& max_fov=0.52);
        /** @fn Eigen::Matrix<double, Eigen::Dynamic, 3> refine_faces_in_view( 
         *          const Eigen::Ref<const Eigen::Vector3d>& pos,
         *          const double& max_fov=0.52)
                
            Use the refinement function for the faces within a certain angle of the position

            @param pos Eigen vector for position
            @param max_fov Maximum cosine(angle) for faces to refine
            @returns Face_centers get the center of the new faces in view

            @author Shankar Kulumani
            @version 28 October 2018
        */
        Eigen::Matrix<double, Eigen::Dynamic, 3> refine_faces_in_view(
                const Eigen::Ref<const Eigen::Vector3d>& pos,
                const double& max_fov=0.52);

        /** @fn bool remesh_faces_in_view(const Eigen::Ref<const Eigen::Vector3d>& pos,
         *      const double& max_fov=0.52, const double& edge_length=0.01)
                
            Use isotropic remeshing to remesh the faces that are in view of the point

            @param pos Position vector
            @param max_fov Field of view to use to compute the faces in view
            @param edge_length Target edge length for the new faces
            @returns bool Success if true

            @author Shankar Kulumani
            @version 28 October 2018
        */
        bool remesh_faces_in_view(
                const Eigen::Ref<const Eigen::Vector3d>& pos,
                const double& max_fov=0.52,
                const double& edge_length=0.01);

        /** @fn Eigen::RowVector3d get_vertex(const Index& index) const
                
            Get the desired vertex given an index

            @param Index Index of vertex to get
            @returns vertex Eigen array of the specific vertex

            @author Shankar Kulumani
            @version 28 October 2018
        */
        template<typename Index>
        Eigen::RowVector3d get_vertex(const Index& index) const;
        
        /** @fn bool set_vertex(const Vertex_index& vd, 
         *              const Eigen::Ref<const Eigen::Vector3d>& vec)
                
            Update the position of a single vertex and recompute the affected
            mesh properties

            @param vd Vertex index to modify
            @param vec Eigen vector of the new vertex position
            @returns bool True if good

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool set_vertex(const Vertex_index& vd,
                const Eigen::Ref<const Eigen::Vector3d>& vec);
        
        /** @fn Eigen::RowVector3i get_face_vertices(const Index& index) const;
                
            Get the vertex indecies of all vertices in the face

            @param index A face index (basically an integer)
            @returns row a row vector of indices

            @author Shankar Kulumani
            @version 6 June 2018
        */
        template<typename Index>
        Eigen::RowVector3i get_face_vertices(const Index& index) const;

        // Getters for the property maps
        template<typename Index>
        Eigen::Vector3d get_face_normal(const Index& fd) const;

        template<typename Index>
        Eigen::Vector3d get_face_center(const Index& fd) const;

        Eigen::Matrix<double, Eigen::Dynamic, 3> get_face_center(const std::vector<Face_index>& face_vec) const;

        Eigen::Matrix<double, Eigen::Dynamic, 3> get_all_face_center( void ) const; 
        
        /** @fn Eigen::VectorXd get_all_face_area( void ) const
                
            Get the area of each face of the mesh

            @param void none
            @returns Vector of all face areas

            @author Shankar Kulumani
            @version 12 June 2018
        */
        Eigen::VectorXd get_all_face_area( void ) const;

        template<typename Index>
        Eigen::Matrix3d get_face_dyad(const Index& fd) const;
        
        template<typename Index>
        Eigen::Vector3d get_halfedge_normal(const Index& hd) const;
        
        template<typename Index>
        Eigen::Matrix3d get_edge_dyad(const Index& ed) const;
        
        template<typename Index>
        double get_edge_factor(const Index& ed) const;

        template<typename Index>
        double get_face_factor(const Index& ed) const;

        double get_sum_face_factor( void ) const;
        
        Mesh surface_mesh; /**< Surface_mesh object */
    private:
        
        /** @fn build_surface_mesh(const Eigen::Ref<const Eigen::MatrixXd>& V, const Eigen::Ref<const Eigen::MatrixXi>& F)
                
            Build the mesh and associated properties using the vertices and faces.

            @param V Input vertices nx3
            @param F Input faces mx3
            @returns none

            @author Shankar Kulumani
            @version 28 October 2018
        */
        void build_surface_mesh(
                const Eigen::Ref<const Eigen::MatrixXd>& V,
                const Eigen::Ref<const Eigen::MatrixXi>& F);
        
        /** @fn bool build_face_properties( void )
                
            Compute the face properties for all faces
            This basically calls compute_face_properties for each face

            @param void 
            @returns bool true if sucess

            @author Shankar Kulumani
            @version Date
        */
        bool build_face_properties( void );
        /** @fn bool update_face_properties(const std::vector<Face_index>& face_vec)
                
            Update the face properties, normal, center, and dyad, of faces listed
            in the vector

            @param face_vec std::vector holding the faces to modify
            @returns bool true if successful

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool update_face_properties(const std::vector<Face_index>& face_vec) ;
        
        /** @fn bool compute_face_properties(const Face_index& fd)
                
            Actually compute the face properties for the given face

            @param fd Face_index of face to modify
            @returns bool true if successful

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool compute_face_properties(const Face_index& fd);
        
        /** @fn bool build_halfedge_properties( void )
                
            Compute the halfedge normals for all halfedges. This basically
            calls compute_halfedge_properties on every halfedge

            @param void 
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 June 22018
        */
        bool build_halfedge_properties( void );
        
        /** @fn bool update_halfedge_properties(
         *      const std::vector<Halfedge_index>& halfedge_vec)
                
            Update the properties for each  halfedge in the vector by calling
            compute_halfedge_properties on each

            @param halfedge_vec std::vector of halfedges to update
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 Juen 2018
        */
        bool update_halfedge_properties(const std::vector<Halfedge_index>& halfedge_vec);

        /** @fn bool compute_halfedge_properties(const Halfedge_index& hd)
                
            Compute the halfedge normal 

            @param hd Halfedge_index
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool compute_halfedge_properties(const Halfedge_index& hd);
        
        /** @fn bool build_edge_properties( void )
                
            Compute edge properties for all edges. This just calls
            compute_edge_properties for each edge

            @param void 
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 June 2018
        */         
        bool build_edge_properties( void ); 

        /** @fn bool update_edge_properties(const std::vector<Edge_index>& edge_vec)
                
            Update the edge properties for each edge in the vector

            @param edge_vec vector of edge indices
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool update_edge_properties( const std::vector<Edge_index>& edge_vec);

        /** @fn bool compute_edge_properties(const Edge_index& ed)
                
            Cmopute the edge dyad for the given edge index

            @param ed Edge index
            @returns bool true if success

            @author Shankar Kulumani
            @version 9 June 2018
        */
        bool compute_edge_properties(const Edge_index& ed);

        /** @fn std::vector<Face_index> get_faces_with_vertex(const Vertex_index& vd)
                
            Get all face indices which contain this vertex vd

            @param vd Vertex index of the vertex in question
            @returns fd_vec Face index std::vector of all faces

            @author Shankar Kulumani
            @version 9 June 2018
        */
        std::vector<Face_index> get_faces_with_vertex(const Vertex_index& vd) const;
        
        /** @fn std::vector<Halfedge_index> get_halfedges_with_vertex(
         *                      const Vertex_index& vd) const
                
            Get all the halfedges with this vertex as a source or target

            @param vd Vertex index
            @returns halfedge_vec std:;vector of all halfedges associated with 
            this vertex

            @author Shankar Kulumani
            @version 9 June 2018
        */
        std::vector<Halfedge_index> get_halfedges_with_vertex(const Vertex_index& vd) const;
        
        /** @fn std::vector<Edge_index> get_edges_with_vertex(const Vertex_index& vd) const
                
            Get all edges that are associated with this vertex

            @param vd Vertex index
            @returns edge_vec vector of edge indices

            @author Shankar Kulumani
            @version 9 June 2018
        */
        std::vector<Edge_index> get_edges_with_vertex(const Vertex_index& vd) const;


};


#endif
