/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DirichletBC2FEM.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "DirichletBC2FEM.h"

#include "GeoLib/Point.h"
#include "MeshLib/Tools/Tools.h"


namespace FemLib
{

DirichletBC2FEM::DirichletBC2FEM(const MeshLib::IMesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values)
{
    // pickup nodes on geometry
    MeshLib::findNodesOnGeometry(&msh, &geo, &vec_nodes);
    const size_t n_bc_nodes = vec_nodes.size();

    if (n_bc_nodes>0) {
        // set values
        vec_values.resize(n_bc_nodes);
        for (size_t i=0; i<n_bc_nodes; i++) {
            const GeoLib::Point* x = msh.getNodeCoordinatesRef(vec_nodes[i]);
            bc_func.eval(x->getData(), vec_values[i]);
        }
    } else {
        std::cout << "Dirichlet B.C. was not found." << std::endl;
    }
}

} //end
