/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FemIC.h"

#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/IVector.h"
#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "NumLib/Fem/BC/IC2FEM.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

namespace SolutionLib
{

///
void FemIC::addDistribution(const GeoLib::GeoObject* geo, const NumLib::ITXFunction* ic_func)
{
    assert(geo!=nullptr);
    assert(ic_func!=nullptr);
    _vec_geo.push_back(geo);
    _vec_func.push_back(ic_func);
}

/// setup 
void FemIC::setup(MathLib::IVector &u0) const
{
    if (_vec_geo.size()==0)
        WARN("***WARN: IC not found.");

    for (size_t i=0; i<_vec_geo.size(); i++) {
        std::vector<size_t> vec_node_id;
        std::vector<double> vec_node_value;
        NumLib::IC2FEM ic2fem(*_mshNodeSearcher, *_vec_geo[i], *_vec_func[i], vec_node_id, vec_node_value);
        for (size_t j=0; j<vec_node_id.size(); j++) {
            u0.set(vec_node_id[j], vec_node_value[j]);
        }
    }
    
}

}
