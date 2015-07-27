/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "NeumannBC2FEM.h"

#include <map>

#include "MathLib/DataType.h"

#include "GeoLib/GeoObject.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

#include "NumLib/Fem/Tools/IFeObjectContainer.h"
#include "NumLib/Fem/Integration/IIntegration.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Function/TXPosition.h"

namespace NumLib
{

NeumannBC2FEM::NeumannBC2FEM(
        MeshGeoToolsLib::MeshNodeSearcher& nodeSearcher, MeshGeoToolsLib::BoundaryElementsSearcher& beSearcher,
        const double &current_time,
        IFeObjectContainer &feObjects, const GeoLib::GeoObject &_geo,
        const NumLib::ITXFunction &_bc_func,
        std::vector<size_t> &_vec_nodes, std::vector<double> &_vec_values)
{
    auto const* msh = nodeSearcher.getMesh();
    // pickup nodes on geo
    _vec_nodes = nodeSearcher.getMeshNodeIDs(_geo);

    // distribute to RHS
    _vec_values.resize(_vec_nodes.size());
    switch (_geo.getGeoType())
    {
        case  GeoLib::GEOTYPE::POINT:
        {
            // no need to integrate
            // get discrete values at nodes
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                auto* x = msh->getNode(_vec_nodes[i])->getCoords();
                NumLib::TXPosition pos(current_time, x);
                _bc_func.eval(pos, _vec_values[i]);
            }
            break;
        }
        case GeoLib::GEOTYPE::POLYLINE:
        {
            // find edge elements on the geo
            std::vector<MeshLib::Element*> vec_edge_eles;
            vec_edge_eles = beSearcher.getBoundaryElements(_geo);
            // for each edge elements found
            std::map<size_t, double> map_nodeId2val;
            for (auto edge : vec_edge_eles) {
                for (std::size_t i=0; i<edge->getNNodes(); i++)
                    map_nodeId2val[edge->getNodeIndex(i)] = 0.0;
            }
            for (auto e : vec_edge_eles)
            {
                //e->setCurrentOrder(msh.getCurrentOrder());
                const size_t edge_nnodes = e->getNNodes();
                // set values at nodes
                MathLib::LocalVector  nodal_val(edge_nnodes);
                for (size_t i_nod=0; i_nod<edge_nnodes; i_nod++) {
                    NumLib::TXPosition pos(current_time, e->getNode(i_nod)->getCoords());
                    double v = .0;
                    _bc_func.eval(pos, v);
                    nodal_val[i_nod] = v;
                }
                // compute integrals
                IFiniteElement *fe_edge = feObjects.getFeObject(*e);
                fe_edge->setMeshElement(*e);
                auto q = NumLib::getIntegrationMethod(e->getCellType());
                MathLib::LocalMatrix M = MathLib::LocalMatrix::Zero(edge_nnodes, edge_nnodes);
                //double x_ref[3];
                NumLib::DynamicShapeMatrices shapeMat(e->getDimension(), msh->getDimension(), e->getNNodes());
                for (size_t j=0; j<q->getNPoints(); j++) {
                    shapeMat.setZero();
                    auto wp = q->getWeightedPoint(j);
                    fe_edge->computeShapeFunctionsd(wp.getCoords(), shapeMat);
                    M.noalias() += (shapeMat.N * shapeMat.N.transpose()) * shapeMat.detJ * wp.getWeight();
                }
                delete q;
                MathLib::LocalVector result = M * nodal_val;
                // add into RHS values
                for (size_t k=0; k<edge_nnodes; k++) {
                    map_nodeId2val[e->getNodeIndex(k)] += result[k];
                }
            }
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                _vec_values[i] = map_nodeId2val[_vec_nodes[i]];
            }
            break;
        }
        default:
            throw "Given GeoType is not supported yet in NeumannBC2FEM";
            break;
    }
}

}
