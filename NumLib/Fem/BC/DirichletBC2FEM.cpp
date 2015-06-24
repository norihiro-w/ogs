/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "DirichletBC2FEM.h"

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "NumLib/Function/ITXFunction.h"

namespace NumLib
{

DirichletBC2FEM::DirichletBC2FEM(const MeshLib::Mesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values)
{
    // pickup nodes on geometry
	MeshGeoToolsLib::MeshNodeSearcher nodeSearcher(msh);
	vec_nodes = nodeSearcher.getMeshNodeIDs(geo);
    const size_t n_bc_nodes = vec_nodes.size();

    if (n_bc_nodes>0) {
        // set values
        vec_values.resize(n_bc_nodes);
        for (size_t i=0; i<n_bc_nodes; i++) {
            bc_func.eval(msh.getNode(vec_nodes[i])->getCoords(), vec_values[i]);
        }
    } else {
        std::cout << "Dirichlet B.C. was not found." << std::endl;
    }
}

} //end
