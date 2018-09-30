#include "mesh.hpp"
#include "polyhedron.hpp"

#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>

#include <tuple>
#include <assert.h>
#include <cmath>

// TODO COMPLETELY REMOVE POLYHEDRON
// Member methods
MeshData::MeshData(const Eigen::Ref<const Eigen::MatrixXd> &V,
        const Eigen::Ref<const Eigen::MatrixXi> &F) {
    this->build_surface_mesh(V, F);
}

void MeshData::build_surface_mesh(const Eigen::Ref<const Eigen::MatrixXd> & V,
        const Eigen::Ref<const Eigen::MatrixXi>& F) {
    // create some vertices
    std::vector<Vertex_index> vert_indices;
    // build the mesh
    // build vector of vertex descriptors

    for (int ii = 0; ii < V.rows(); ++ii) {
        Point p = Kernel::Point_3(V(ii, 0), V(ii, 1), V(ii, 2));
        Vertex_index v = this->surface_mesh.add_vertex(p);
        assert(surface_mesh.is_valid(v));
        vert_indices.push_back(v);
    }

    for (int ii = 0; ii < F.rows(); ++ii) {
        Point p1, p2, p3;
        p1 = Kernel::Point_3(V(F(ii, 0), 0), V(F(ii, 0), 1), V(F(ii, 0), 2));
        p2 = Kernel::Point_3(V(F(ii, 1), 0), V(F(ii, 1), 1), V(F(ii, 2), 2));
        p3 = Kernel::Point_3(V(F(ii, 2), 0), V(F(ii, 2), 1), V(F(ii, 2), 2));
        
        Vertex_index v1, v2, v3;

        v1 = vert_indices[F(ii, 0)];
        v2 = vert_indices[F(ii, 1)];
        v3 = vert_indices[F(ii, 2)];

        Face_index f = this->surface_mesh.add_face(v1, v2, v3);
        assert(surface_mesh.is_valid(f));
    }

    /* // can also loop over edges */
    /* for(Edge_index ed: surface_mesh.edges()) { */
    /*     // returns a vertex_index for the start/finish of the edge */
    /*     vds = source(ed, surface_mesh); */
    /*     vde = end(ed, surface_mesh); */
    /* } */

    assert(surface_mesh.is_valid());
    /* std::vector<std::string> props = surface_mesh.properties<Face_index>(); */
    
    /* BOOST_FOREACH(std::string p, props){ */
    /*     std::cout << p << std::endl; */
    /* } */

    /* std::cout << surface_mesh.properties()[0] << std::endl; */
    // face_normal, edge normal, halfedge normal, face dyad, edge dyad
    bool face_built, halfedge_built, edge_built;
    #pragma omp parallel if (false)
    {
        #pragma omp single
        {
            #pragma omp task depend(out:face_built)
            {face_built = build_face_properties();}

            #pragma omp task depend(out:halfedge_built)
            {halfedge_built = build_halfedge_properties();}

            #pragma omp task depend(in:face_built) depend(in:halfedge_built)
            {edge_built = build_edge_properties();}
        }
    }
}

bool MeshData::build_face_properties( void ) {
    // loop over all faces need to dereference the iterators but not the index
    for (Face_index fd: surface_mesh.faces() ){
        bool face_updated = compute_face_properties(fd);
        assert(face_updated);
    }
    return true;
}

bool MeshData::build_halfedge_properties( void ) {
    for (Halfedge_index hd: surface_mesh.halfedges()) {
        compute_halfedge_properties(hd);
    }
    return true;
}

bool MeshData::build_edge_properties( void ){
    for (Edge_index ed: surface_mesh.edges()) {
        compute_edge_properties(ed);
    }
    return true;
}

bool MeshData::build_edge_factor( const Eigen::Ref<const Eigen::Vector3d>& pos ) {
    // test if property map exists if not then create
    Mesh::Property_map<Edge_index, double> edge_factor;
    bool found, created;
    std::tie(edge_factor, found) = surface_mesh.add_property_map<
        Edge_index, double>("e:edge_factor");

    // loop over edges
    for (Edge_index ed : surface_mesh.edges()) {

        // get vertex of each edge endpoitn
        Vertex_index v1, v2;
        v1 = surface_mesh.vertex(ed, 0);
        v2 = surface_mesh.vertex(ed, 1);

        Eigen::Vector3d vec1, vec2;
        vec1 = get_vertex(v1);
        vec2 = get_vertex(v2);
        // subtract from state adn find norm
        double r1, r2;
        r1 = (vec1 - pos).norm();
        r2 = (vec2 - pos).norm();
        // find length of edge
        double e = (vec1 - vec2).norm();
        // take natural log
        edge_factor[ed] = std::log((r1 + r2 + e)/(r1 + r2 - e));
    }
    return true;
}

bool MeshData::build_face_factor(const Eigen::Ref<const Eigen::Vector3d>& pos) {
    Mesh::Property_map<Face_index, double> face_factor;
    bool created, found;
    std::tie(face_factor, found) = surface_mesh.property_map<
        Face_index, double>("f:face_factor");
    if (!found) {
        // face factor w property map
        std::tie(face_factor, created) = surface_mesh.add_property_map<Face_index, double>(
                "f:face_factor", 0);
        assert(created);
    }
    
    for (Face_index fd: surface_mesh.faces()) {
        Halfedge_index h1, h2, h3;
        h1 = surface_mesh.halfedge(fd);
        h2 = surface_mesh.next(h1);
        h3 = surface_mesh.next(h2);
        assert(surface_mesh.next(h3) == h1);

        Vertex_index v1, v2, v3;
        v1 = surface_mesh.source(h1);
        v2 = surface_mesh.source(h2);
        v3 = surface_mesh.source(h3);
        assert(surface_mesh.target(h1) == v2);
        assert(surface_mesh.target(h2) == v3);
        assert(surface_mesh.target(h3) == v1);

        // now extract the point into Eigen arrays
        Eigen::Vector3d vec1, vec2, vec3, r1, r2, r3;
        vec1 = get_vertex(v1);
        vec2 = get_vertex(v2);
        vec3 = get_vertex(v3);
        r1 = vec1 - pos;
        r2 = vec2 - pos;
        r3 = vec3 - pos;

        double num, den;
        num = r1.dot(r2.cross(r3));
        den = r1.norm() * r2.norm() * r3.norm() 
            + r1.norm() * r2.dot(r3) 
            + r2.norm() * r3.dot(r1)
            + r3.norm() * r1.dot(r2);

        face_factor[fd] = 2.0 * atan2(num, den);
    }

    return true;
}

void MeshData::update_mesh(const Eigen::Ref<const Eigen::MatrixXd> &V, 
        const Eigen::Ref<const Eigen::MatrixXi> &F) {
    // update the polyhedron and surface mesh
    // clear the mesh
    this->surface_mesh.clear();

    this->build_surface_mesh(V, F);
}

// MeshData Getters
Eigen::Matrix<double, Eigen::Dynamic, 3> MeshData::get_verts( void ) const {
    // extract out vertices from surface_mesh
    std::size_t num_v = surface_mesh.number_of_vertices();
    std::size_t row_index = 0;

    Eigen::Matrix<double, Eigen::Dynamic, 3> vertices(num_v, 3);

    for (Vertex_index vd: surface_mesh.vertices() ) {
        Eigen::Matrix<double, 1, 3> row; 
        row << surface_mesh.point(vd).x(), surface_mesh.point(vd).y(), surface_mesh.point(vd).z();
        
        vertices.row(row_index) = row;
        ++row_index;
    }

    return vertices;
}

Eigen::Matrix<int, Eigen::Dynamic, 3> MeshData::get_faces( void ) const {
    std::size_t num_f = surface_mesh.number_of_faces();
    std::size_t row_index = 0;

    Eigen::Matrix<int, Eigen::Dynamic, 3> faces(num_f, 3);

    for ( Face_index fd: surface_mesh.faces() ) {
        std::size_t col_index = 0;
        for (Vertex_index vd: vertices_around_face(surface_mesh.halfedge(fd), surface_mesh)) {
            faces(row_index, col_index) = (int)vd;
            ++col_index;
        }
        ++row_index;
    }

    return faces;
}

std::size_t MeshData::number_of_vertices( void ) const {
    return surface_mesh.number_of_vertices();
}

std::size_t MeshData::number_of_edges( void ) const {
    return surface_mesh.number_of_edges();
}

std::size_t MeshData::number_of_faces( void ) const {
    return surface_mesh.number_of_faces();
}

std::size_t MeshData::number_of_halfedges( void ) const {
    return surface_mesh.number_of_halfedges();
}

// Range iterators
Mesh::Vertex_range MeshData::vertices( void ) const {
    return surface_mesh.vertices();
}

Mesh::Face_range MeshData::faces( void ) const {
    return surface_mesh.faces();
}

Mesh::Edge_range MeshData::edges( void ) const {
    return surface_mesh.edges();
}

Mesh::Halfedge_range MeshData::halfedges( void ) const {
    return surface_mesh.halfedges();
}

template<typename Index>
Eigen::Vector3d MeshData::get_face_normal(const Index& fd_in) const {
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_unit_normal;
    bool found;
    std::tie(face_unit_normal, found) = surface_mesh.property_map<
        Face_index,  Eigen::Vector3d>("f:face_unit_normal");
    assert(found);
    Face_index fd(fd_in);
    return face_unit_normal[fd];
}

template<typename Index>
Eigen::Vector3d MeshData::get_face_center(const Index& fd_in) const {
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_center;
    bool found;
    std::tie(face_center, found) = surface_mesh.property_map<
        Face_index, Eigen::Vector3d>("f:face_center");
    assert(found);
    Face_index fd(fd_in);
    return face_center[fd];
}

Eigen::Matrix<double, Eigen::Dynamic, 3> MeshData::get_face_center(
        const std::vector<Face_index>& face_vec) const {
    
    Eigen::Matrix<double, Eigen::Dynamic, 3> face_centers(face_vec.size(), 3);
    
    std::size_t index = 0;
    for (Face_index fd: face_vec) {
        face_centers.row(index) = get_face_center(fd).transpose();
        ++index;
    }
    return face_centers;
}

Eigen::Matrix<double, Eigen::Dynamic, 3> MeshData::get_all_face_center( void ) const {
    Eigen::Matrix<double, Eigen::Dynamic, 3> face_centers(number_of_faces(), 3);

    std::size_t index = 0;
    for (Face_index fd: surface_mesh.faces()) {
        face_centers.row(index) = get_face_center(fd).transpose();
        ++index;
    }
    return face_centers;
}

Eigen::VectorXd MeshData::get_all_face_area( void ) const {
    Eigen::VectorXd face_area(number_of_faces(), 1);
    std::size_t index = 0;
    for (Face_index fd: surface_mesh.faces()) {
        face_area(index) = CGAL::Polygon_mesh_processing::face_area(fd, 
                surface_mesh);
        ++index;
    }
    return face_area;
}

template<typename Index>
Eigen::Matrix3d MeshData::get_face_dyad(const Index& fd_in) const {
    Mesh::Property_map<Face_index, Eigen::Matrix3d> face_dyad;
    bool found;
    std::tie(face_dyad, found) = surface_mesh.property_map<
        Face_index, Eigen::Matrix3d>("f:face_dyad");
    assert(found);
    Face_index fd(fd_in);
    return face_dyad[fd];
}

template<typename Index>
Eigen::Vector3d MeshData::get_halfedge_normal(const Index& hd_in) const {
    Mesh::Property_map<Halfedge_index, Eigen::Vector3d> halfedge_unit_normal;
    bool found;
    std::tie(halfedge_unit_normal, found) = surface_mesh.property_map<
        Halfedge_index, Eigen::Vector3d>("h:halfedge_unit_normal");
    Halfedge_index hd(hd_in);
    return halfedge_unit_normal[hd];
}

template<typename Index>
Eigen::Matrix3d MeshData::get_edge_dyad(const Index& ed_in) const {
    Mesh::Property_map<Edge_index, Eigen::Matrix3d> edge_dyad;
    bool found;
    std::tie(edge_dyad, found) = surface_mesh.property_map<
        Edge_index, Eigen::Matrix3d>("e:edge_dyad");
    Edge_index ed(ed_in);
    return edge_dyad[ed];
}

template<typename Index>
double MeshData::get_edge_factor(const Index& ed_in) const {
    Mesh::Property_map<Edge_index, double> edge_factor;
    bool found;
    std::tie(edge_factor, found) = surface_mesh.property_map<
        Edge_index, double>("e:edge_factor");
    Edge_index ed(ed_in);
    return edge_factor[ed];
}

template<typename Index>
double MeshData::get_face_factor(const Index& fd_in) const {
    Mesh::Property_map<Face_index, double> face_factor;
    bool found;
    std::tie(face_factor, found) = surface_mesh.property_map<
        Face_index, double>("f:face_factor");
    Face_index fd(fd_in);
    assert(found);
    return face_factor[fd];
}

double MeshData::get_sum_face_factor( void ) const {
    double w_face_sum = 0;
    for (Face_index fd : surface_mesh.faces()) {
        w_face_sum += get_face_factor(fd);    
    }
    return w_face_sum;
}

template<typename Index>
Eigen::RowVector3d MeshData::get_vertex(const Index& index) const {
    // form Vertex_index
    Vertex_index vd(index);
    
    // create an array and return
    Eigen::RowVector3d vertex;
    vertex(0) = surface_mesh.point(vd).x();
    vertex(1) = surface_mesh.point(vd).y();
    vertex(2) = surface_mesh.point(vd).z();

    return vertex;
}

bool MeshData::set_vertex(const Vertex_index& vd, 
        const Eigen::Ref<const Eigen::Vector3d>& vec) {
    
    Point p = Kernel::Point_3(vec(0), vec(1), vec(2));
    surface_mesh.point(vd) = p;

    // update the mesh properties associated with this vertex index
    std::vector<Face_index> face_vec = get_faces_with_vertex(vd);
    std::vector<Halfedge_index> halfedge_vec = get_halfedges_with_vertex(vd);
    std::vector<Edge_index> edge_vec = get_edges_with_vertex(vd);

    update_face_properties(face_vec);
    update_halfedge_properties(halfedge_vec);
    update_edge_properties(edge_vec);
}

bool MeshData::refine_faces(const std::vector<Face_index>& face_vec,
        std::vector<Face_index>& new_faces,
        std::vector<Vertex_index>& new_vertices,
        const int& density) {
    
    CGAL::Polygon_mesh_processing::refine(
            surface_mesh,
            face_vec,
            std::back_inserter(new_faces),
            std::back_inserter(new_vertices),
            CGAL::Polygon_mesh_processing::parameters::density_control_factor(density));

    surface_mesh.collect_garbage();
    build_face_properties();
    build_halfedge_properties();
    build_edge_properties();
    surface_mesh.collect_garbage();
    return true;
}

bool MeshData::remesh_faces(const std::vector<Face_index>& face_vec,
        const double& target_edge_length,
        const int& number_of_iterations) {

    CGAL::Polygon_mesh_processing::isotropic_remeshing(
            face_vec,
            target_edge_length,
            surface_mesh,
            CGAL::Polygon_mesh_processing::parameters::number_of_iterations(number_of_iterations));

    // now need to update all the face, edge, halfedge properties
    surface_mesh.collect_garbage();
    build_face_properties();
    build_halfedge_properties();
    build_edge_properties();
    surface_mesh.collect_garbage();
    return true;
}

std::vector<Face_index> MeshData::faces_in_fov(
        const Eigen::Ref<const Eigen::Vector3d>& pos,
        const double& max_fov) {
    std::vector<Face_index> faces_in_view;
    Eigen::Vector3d pos_uvec = pos.normalized();
    // loop over all the faces 
    for (Face_index fd : surface_mesh.faces()){
        // get a vertex in this face
        Halfedge_index hd = surface_mesh.halfedge(fd);
        Vertex_index vd = surface_mesh.source(hd);
        double angle = std::acos(pos_uvec.dot(get_vertex(vd).normalized()));
        if (angle < max_fov) {
            faces_in_view.push_back(fd);
        }
    }
    return faces_in_view;
}

std::vector<Vertex_index> MeshData::vertices_in_fov(
        const Eigen::Ref<const Eigen::Vector3d>& pos,
        const double& max_fov) {
    std::vector<Vertex_index> vertices_in_view;
    Eigen::Vector3d pos_uvec = pos.normalized();
    for (Vertex_index vd : surface_mesh.vertices()) {
        double angle = std::acos(pos_uvec.dot(get_vertex(vd).normalized()));
        if (angle < max_fov) {
            vertices_in_view.push_back(vd);
        }
    }
    return vertices_in_view;
}

Eigen::Matrix<double, Eigen::Dynamic, 3> MeshData::refine_faces_in_view(
        const Eigen::Ref<const Eigen::Vector3d>& pos,
        const double& max_fov) {
    // find current faces in fov
    std::vector<Face_index> faces_to_refine = faces_in_fov(pos, max_fov);
    // remesh those faces
    std::vector<Face_index> new_faces;
    std::vector<Vertex_index> new_vertices;
    refine_faces(faces_to_refine, new_faces, new_vertices, 4.0);
    /* remesh_faces(faces_to_refine, 0.01, 3); */
    // get the face centers
    Eigen::Matrix<double, Eigen::Dynamic, 3> new_face_centers;
    new_face_centers = get_face_center(new_faces);
    // return 
    return new_face_centers;
}

bool MeshData::remesh_faces_in_view(
        const Eigen::Ref<const Eigen::Vector3d>& pos,
        const double& max_fov, const double& edge_length) {
    std::vector<Face_index> faces_to_remesh = faces_in_fov(pos, max_fov);
    return remesh_faces(faces_to_remesh, edge_length, 3);
}

template<typename Index>
Eigen::RowVector3i MeshData::get_face_vertices(const Index& index) const {
    Face_index fd(index);

    Eigen::RowVector3i face_vertices;
    std::size_t col_index = 0;
    for (Vertex_index vd: vertices_around_face(surface_mesh.halfedge(fd), surface_mesh)) {
        face_vertices(col_index) = (int)vd;
        ++col_index;
    }
    return face_vertices;
}

std::vector<Face_index> MeshData::get_faces_with_vertex(const Vertex_index& vd) const {
    // loop around the vertex and get all the faces
    std::vector<Face_index> face_vec;

    for(Face_index fd : faces_around_target(surface_mesh.halfedge(vd), surface_mesh)) {
        face_vec.push_back(fd); 
    }
    return face_vec;
}

std::vector<Halfedge_index> MeshData::get_halfedges_with_vertex(
        const Vertex_index& vd) const {
    std::vector<Halfedge_index> halfedge_vec;

    for(Halfedge_index hd : halfedges_around_target(vd, surface_mesh)) {
        halfedge_vec.push_back(hd);
        halfedge_vec.push_back(surface_mesh.opposite(hd));
    }

    return halfedge_vec;
}

std::vector<Edge_index> MeshData::get_edges_with_vertex(const Vertex_index& vd) const {
    std::vector<Edge_index> edge_vec;

    for(Halfedge_index hd : halfedges_around_target(vd, surface_mesh)) {
        edge_vec.push_back(surface_mesh.edge(hd));
    }
    return edge_vec;
}

bool MeshData::update_face_properties(const std::vector<Face_index>& face_vec) {
    // loop over the face vec and update the properties of each
    for (Face_index fd : face_vec) {
        compute_face_properties(fd);
    }
    return true;
}

bool MeshData::compute_face_properties(const Face_index& fd) {
    // Build property maps for the surface mesh
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_unit_normal;
    bool created;
    std::tie( face_unit_normal, created ) 
        = surface_mesh.add_property_map<Face_index, Eigen::Vector3d>(
                "f:face_unit_normal", (Eigen::Vector3d() << 0, 0, 0).finished());
    /* assert(created); */
    
    // Center face property map
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_center;
    /* bool created; */
    std::tie(face_center, created) 
        = surface_mesh.add_property_map<Face_index, Eigen::Vector3d>(
                "f:face_center", (Eigen::Vector3d() << 0 ,0 ,0).finished());
    /* assert(created); */
    
    // Face dyad
    Mesh::Property_map<Face_index, Eigen::Matrix3d> face_dyad;
    std::tie(face_dyad, created)
        = surface_mesh.add_property_map<Face_index, Eigen::Matrix3d>(
                "f:face_dyad", Eigen::Matrix3d::Zero());
    /* assert(created); */ // true if new, false if exists but returned
    // need to consecutive vertices in face to get teh normal
    Halfedge_index h1, h2, h3;
    h1 = surface_mesh.halfedge(fd);
    h2 = surface_mesh.next(h1);
    h3 = surface_mesh.next(h2);
    assert(surface_mesh.next(h3) == h1);

    Vertex_index v1, v2, v3;
    v1 = surface_mesh.source(h1);
    v2 = surface_mesh.source(h2);
    v3 = surface_mesh.source(h3);
    assert(surface_mesh.target(h1) == v2);
    assert(surface_mesh.target(h2) == v3);
    assert(surface_mesh.target(h3) == v1);

    // now extract the point into Eigen arrays
    Eigen::Vector3d vec1, vec2, vec3;
    vec1 = get_vertex(v1);
    vec2 = get_vertex(v2);
    vec3 = get_vertex(v3);

    Eigen::Vector3d edge1, edge2;
    edge1 = vec2 - vec1;
    edge2 = vec3 - vec1;

    face_unit_normal[fd] = edge1.cross(edge2).normalized();
    face_center[fd] = 1.0 / 3.0 * (vec1 + vec2 + vec3);
    face_dyad[fd] = face_unit_normal[fd] * face_unit_normal[fd].transpose();

    assert(face_dyad[fd].isApprox(face_dyad[fd].transpose(), 1e-3));
    return true;
}

bool MeshData::update_halfedge_properties(const std::vector<Halfedge_index>& halfedge_vec) {
    for (Halfedge_index hd : halfedge_vec) {
        compute_halfedge_properties(hd); 
    }
    return true;
}

bool MeshData::compute_halfedge_properties(const Halfedge_index& hd) {
    Mesh::Property_map<Halfedge_index, Eigen::Vector3d> halfedge_unit_normal;
    bool created;
    std::tie(halfedge_unit_normal, created) 
        = surface_mesh.add_property_map<Halfedge_index, Eigen::Vector3d>(
                "h:halfedge_unit_normal", (Eigen::Vector3d() << 0, 0, 0).finished());
    
    /* assert(created); */
    
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_unit_normal;
    bool found;
    std::tie(face_unit_normal, found) 
        = surface_mesh.property_map<Face_index, Eigen::Vector3d>(
                "f:face_unit_normal");
    assert(found);

    Vertex_index vs, ve;
    vs = surface_mesh.source(hd);
    ve = surface_mesh.target(hd);

    Eigen::Vector3d vec_start, vec_end, vec_edge, face_normal;
    vec_start = get_vertex(vs);
    vec_end = get_vertex(ve);
    vec_edge = vec_end - vec_start;

    Face_index fd;
    fd = surface_mesh.face(hd);
    face_normal = face_unit_normal[fd];

    halfedge_unit_normal[hd] = vec_edge.cross(face_normal).normalized();
    
    return true;
}

bool MeshData::update_edge_properties(const std::vector<Edge_index>& edge_vec) {
    for(Edge_index ed : edge_vec) {
        compute_edge_properties(ed);
    }
    return true;
}

bool MeshData::compute_edge_properties(const Edge_index& ed) {
    // edge dyad
    Mesh::Property_map<Edge_index, Eigen::Matrix3d> edge_dyad;
    bool created;
    std::tie(edge_dyad, created)
        = surface_mesh.add_property_map<Edge_index, Eigen::Matrix3d>(
                "e:edge_dyad", Eigen::Matrix3d::Zero());
    /* assert(created); */
    
    // normal face property map
    Mesh::Property_map<Face_index, Eigen::Vector3d> face_unit_normal;
    bool found;
    std::tie(face_unit_normal, found) 
        = surface_mesh.property_map<Face_index, Eigen::Vector3d>(
                "f:face_unit_normal");
    /* assert(found); */

    // halfedge normal property map
    Mesh::Property_map<Halfedge_index, Eigen::Vector3d> halfedge_unit_normal;
    std::tie(halfedge_unit_normal, found)
        = surface_mesh.property_map<Halfedge_index, Eigen::Vector3d>(
                "h:halfedge_unit_normal");
    assert(found);

    Halfedge_index h1, h2;
    h1 = surface_mesh.halfedge(ed, 0);
    h2 = surface_mesh.halfedge(ed, 1);

    Face_index f1, f2;
    f1 = surface_mesh.face(h1);
    f2 = surface_mesh.face(h2);

    edge_dyad[ed] = face_unit_normal[f1] * halfedge_unit_normal[h1].transpose() 
        + face_unit_normal[f2] * halfedge_unit_normal[h2].transpose();
    // doesn't work for matrices close to zero
    /* assert((edge_dyad[ed] - edge_dyad[ed].transpose()).isApprox(Eigen::Matrix3d::Zero(), 1e-3)); */
    return true;
}

// Template Specialization
template Eigen::RowVector3d MeshData::get_vertex<std::size_t>(const std::size_t&) const;
template Eigen::RowVector3d MeshData::get_vertex<int>(const int&) const;
template Eigen::RowVector3d MeshData::get_vertex<Vertex_index>(const Vertex_index&) const;

template Eigen::RowVector3i MeshData::get_face_vertices<std::size_t>(const std::size_t&) const;
template Eigen::RowVector3i MeshData::get_face_vertices<int>(const int&) const;
template Eigen::RowVector3i MeshData::get_face_vertices<Face_index>(const Face_index&) const;

template Eigen::Vector3d MeshData::get_face_normal<Face_index>(const Face_index&) const;
template Eigen::Vector3d MeshData::get_face_center<Face_index>(const Face_index&) const;
template Eigen::Matrix3d MeshData::get_face_dyad<Face_index>(const Face_index&) const;

template Eigen::Vector3d MeshData::get_halfedge_normal<Halfedge_index>(const Halfedge_index&) const;

template Eigen::Matrix3d MeshData::get_edge_dyad<Edge_index>(const Edge_index&) const;

template double MeshData::get_edge_factor<Edge_index>(const Edge_index&) const;
template double MeshData::get_face_factor<Face_index>(const Face_index&) const;
