/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDirichletBC.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemDirichletBC.h"

#include <logog/include/logog.hpp>

#include "NumLib/Fem/BC/DirichletBC2FEM.h"

namespace SolutionLib
{

FemDirichletBC::FemDirichletBC(MeshGeoToolsLib::MeshNodeSearcher* mshNodeSearcher, const GeoLib::GeoObject* geo, NumLib::ITXFunction* bc_func)
    : _mshNodeSearcher(mshNodeSearcher), _geo(geo), _bc_func(bc_func)
{
    _is_transient = !bc_func->isTemporallyConst();
    _do_setup = true;
}

FemDirichletBC::FemDirichletBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
    : _mshNodeSearcher(nullptr), _geo(nullptr), _bc_func(nullptr), _vec_nodes(vec_node_id), _vec_values(vec_node_values)
{
    _is_transient = false;
    _do_setup = false;
}

/// setup B.C.
void FemDirichletBC::setup(NumLib::PolynomialOrder order)
{
    if (!_do_setup) return;
    if (!_is_transient) _do_setup = false;

    //_msh->setCurrentOrder(order);
    NumLib::DirichletBC2FEM convert(*_mshNodeSearcher, *_geo, *_bc_func, _vec_nodes, _vec_values);

    if (_vec_nodes.size()==0)
        INFOa("***INFO: No Dirichlet BC found in FemDirichletBC::setup()");

    if (!_is_transient)
        _do_setup = false;
}

void FemDirichletBC::write(std::ostream &os) const
{
    os << "# Dirichlet BC\n";
    for (std::size_t i=0; i<_vec_nodes.size(); i++)
        os << _vec_nodes[i] << " " << _vec_values[i] << "\n";
}

}
