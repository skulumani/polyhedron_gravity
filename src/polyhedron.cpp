#include "polyhedron.hpp"
#include "cgal_types.hpp"


void build_polyhedron_index(Polyhedron &P) {
    std::size_t ii = 0;
    for (Vertex_iterator vert = P.vertices_begin(); vert != P.vertices_end(); ++vert) {
        vert->id() = ii++; 
    }
    ii = 0; // reset the counter
    for (Facet_iterator facet = P.facets_begin(); facet != P.facets_end(); ++facet) {
        facet->id() = ii++;
    }

}

// TODO Add documentation - overload the () operator
template<typename HDS>
void Polyhedron_builder<HDS>::operator() (HDS &hds) {

    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    // create the cgal incremental builder
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    // initialize with #v, #f, #half-edges (optional)
    B.begin_surface(V.rows(), F.rows());

    // add all of the vertices
    for (int ii = 0; ii < V.rows(); ++ii) {
        B.add_vertex(Point(V(ii, 0), V(ii, 1), V(ii, 2)));
    }
    // add all of the faces
    for (int ii = 0; ii < F.rows(); ++ii) {
        B.begin_facet();
        B.add_vertex_to_facet(F(ii, 0));
        B.add_vertex_to_facet(F(ii, 1));
        B.add_vertex_to_facet(F(ii, 2));
        B.end_facet();
    }
    B.end_surface();
}

/* //TODO Add documentation */
/* void print_vertices(Polyhedron& P) { */
/*     CGAL::set_ascii_mode(std::cout); */
/*     std::cout << "#V: " << P.size_of_vertices() << " #F: " */
/*         << P.size_of_facets() << std::endl; */
/*     std::cout << "Printing all vertices" << std::endl; */
/*     std::copy (P.points_begin(), P.points_end(), std::ostream_iterator<Kernel::Point_3>( std::cout, "\n")); */
/* } */

/* //TODO Add documentation */
/* void loop_over_facets(Polyhedron &P) { */
/*     for ( Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f) { */
/*         Halfedge_facet_circulator v = f->facet_begin(); */
/*         std::cout << "Number of vertices around facet: " << CGAL::circulator_size(v) << std::endl; */
/*         do { */
/*             std::cout << "ID: " << v->vertex()->id() << " " << "Vertex: " << v->vertex()->point() << " "; */

/*         } while( ++v != f->facet_begin()); */

/*         std::cout << std::endl; */
/*     } */
    
/* } */

// TODO Add documentation
template<typename VectorType, typename IndexType>
void polyhedron_to_eigen(Polyhedron &P, Eigen::PlainObjectBase<VectorType> &V, Eigen::PlainObjectBase<IndexType> &F) {
    // loop over all the faces first
    
    // create some eigen arrays to store all the vertices
    const unsigned int num_v = P.size_of_vertices();
    const unsigned int num_f = P.size_of_facets();
    
    V.resize(num_v, 3);
    F.resize(num_f, 3);

    std::size_t row, col;
    row = 0;
    col = 0;
    /* loop_over_facets(P); */

    // Build V
    for (Vertex_iterator vert = P.vertices_begin(); vert != P.vertices_end(); ++vert) {
        V(row, 0)  = vert->point().x();
        V(row, 1)  = vert->point().y();
        V(row, 2)  = vert->point().z();
        row += 1;
    }

    // Build F
    row = 0;
    for ( Facet_iterator face = P.facets_begin(); face != P.facets_end(); ++face) {
        Halfedge_facet_circulator vert = face->facet_begin();
        col = 0;
        do {
            F(row, col) = vert->vertex()->id();
            col += 1;
        } while( ++vert != face->facet_begin());
        row += 1;
    }
}

void eigen_to_polyhedron(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Polyhedron &P) {
    Polyhedron_builder<HalfedgeDS> builder(V, F);
    P.delegate(builder);
    // loop and fill the eigen array
    build_polyhedron_index(P);
    /* CGAL_assertion(P.is_valid()); */
}


/* // TODO Add documentation - overload the () operator */
/* template<typename HDS> */
/* void Polyhedron_builder<HDS>::operator() (HDS &hds) { */

/*     typedef typename HDS::Vertex Vertex; */
/*     typedef typename Vertex::Point Point; */

/*     // create the cgal incremental builder */
/*     CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true); */
/*     // initialize with #v, #f, #half-edges (optional) */
/*     B.begin_surface(V.rows(), F.rows()); */

/*     // add all of the vertices */
/*     for (int ii = 0; ii < V.rows(); ++ii) { */
/*         B.add_vertex(Point(V(ii, 0), V(ii, 1), V(ii, 2))); */
/*     } */
/*     // add all of the faces */
/*     for (int ii = 0; ii < F.rows(); ++ii) { */
/*         B.begin_facet(); */
/*         B.add_vertex_to_facet(F(ii, 0)); */
/*         B.add_vertex_to_facet(F(ii, 1)); */
/*         B.add_vertex_to_facet(F(ii, 2)); */
/*         B.end_facet(); */
/*     } */
/*     B.end_surface(); */
/* } */

/* void print_polyhedron_stats(Polyhedron &P) { */
/*     std::cout << "Polyhedron is built in CGAL" << std::endl; */
/*     std::cout << "Valid : " << P.is_valid() << std::endl; */
/*     std::cout << "Vertices : " << P.size_of_vertices() << std::endl; */
/*     std::cout << "Faces : " << P.size_of_facets() << std::endl; */
/*     std::cout << "HalfEdges : " << P.size_of_halfedges() << std::endl; */
/* } */

/* // Explicit initialization of the template */
template void Polyhedron_builder<CGAL::HalfedgeDS_default<CGAL::Simple_cartesian<double>, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int> > >::operator()(CGAL::HalfedgeDS_default<CGAL::Simple_cartesian<double>, CGAL::I_Polyhedron_derived_items_3<CGAL::Polyhedron_items_with_id_3>, std::allocator<int> >&);

template void polyhedron_to_eigen<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >
(CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>, CGAL::Polyhedron_items_with_id_3, CGAL::HalfedgeDS_default, std::allocator<int> >&, 
 Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, 
 Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);

