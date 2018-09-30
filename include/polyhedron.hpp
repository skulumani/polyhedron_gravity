#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include "cgal_types.hpp"

#include <Eigen/Dense>

template<typename VectorType, typename IndexType>
void polyhedron_to_eigen(Polyhedron &P,
        Eigen::PlainObjectBase<VectorType> &V, Eigen::PlainObjectBase<IndexType> &F);


void eigen_to_polyhedron(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Polyhedron &P);

void build_polyhedron_index(Polyhedron &poly);

/* void print_polyhedron_stats(CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>, CGAL::Polyhedron_items_with_id_3> &P); */

// declaration for the polyhedron builder
template<typename HDS> 
class Polyhedron_builder : public CGAL::Modifier_base<HDS> {
    public:
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;

        Polyhedron_builder(const Eigen::MatrixXd &V_input, const Eigen::MatrixXi &F_input) : V(V_input), F(F_input) {}
    
        void operator() (HDS &hds);		
};

#endif
