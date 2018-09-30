#ifndef CGAL_TYPES_H
#define CGAL_TYPES_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

// AABB tree for distance stuff
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// dD Spatial searching
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

typedef CGAL::Simple_cartesian<double>     Kernel;

// Basic primitive objects in CGAL
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Ray_3 Ray;

// polyhedron typedef
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>         Polyhedron;
typedef Polyhedron::Facet_iterator          Facet_iterator;
typedef Polyhedron::Vertex_iterator         Vertex_iterator;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

// Surface mesh typedef
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Face_index Face_index;
typedef Mesh::Edge_index Edge_index;
typedef Mesh::Halfedge_index Halfedge_index;

// AABB typedefs
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
typedef boost::optional< AABB_Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;
typedef boost::optional< AABB_Tree::Intersection_and_primitive_id<Plane>::Type > Plane_intersection;
typedef boost::optional< AABB_Tree::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;
typedef AABB_Tree::Primitive_id Primitive_id;

// dD spatial searching
typedef boost::graph_traits<Mesh>::vertex_descriptor Mesh_Point;
typedef boost::graph_traits<Mesh>::vertices_size_type size_type;

typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type Vertex_point_pmap;
typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
typedef CGAL::Search_traits_adapter<Mesh_Point,Vertex_point_pmap,Traits_base> KD_Traits;
typedef CGAL::Orthogonal_k_neighbor_search<KD_Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         KD_Tree;
typedef KD_Tree::Splitter                                                  Splitter;
typedef K_neighbor_search::Distance                                     Distance;

#endif
