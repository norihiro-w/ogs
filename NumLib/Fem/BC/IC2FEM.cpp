/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */


#include "IC2FEM.h"

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "NumLib/Function/ITXFunction.h"

namespace NumLib
{

///
IC2FEM::IC2FEM(MeshGeoToolsLib::MeshNodeSearcher& nodeSearcher, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &ic_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values)
{
    auto* msh = nodeSearcher.getMesh();
    // pickup nodes on geometry
    vec_nodes = nodeSearcher.getMeshNodeIDs(geo);
    const size_t n_bc_nodes = vec_nodes.size();

    if (n_bc_nodes>0) {
        // set values
        vec_values.resize(n_bc_nodes);
        for (size_t i=0; i<n_bc_nodes; i++) {
            ic_func.eval(NumLib::TXPosition(NumLib::TXPosition::Node, vec_nodes[i], msh->getNode(vec_nodes[i])->getCoords()) , vec_values[i]);
        }
    } else {
        std::cout << "I.C. was not found." << std::endl;
    }
}


}
