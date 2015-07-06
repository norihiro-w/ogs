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

FemDirichletBC::FemDirichletBC(const MeshLib::Mesh* msh, const GeoLib::GeoObject* geo, NumLib::ITXFunction* bc_func)
    : _msh(msh), _geo(geo), _bc_func(bc_func)
{
    _is_transient = !bc_func->isTemporallyConst();
    _do_setup = true;
}

FemDirichletBC::FemDirichletBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
    : _msh(NULL), _geo(NULL), _bc_func(NULL), _vec_nodes(vec_node_id), _vec_values(vec_node_values)
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
    NumLib::DirichletBC2FEM convert(*_msh, *_geo, *_bc_func, _vec_nodes, _vec_values);

    if (_vec_nodes.size()==0)
        INFO("***INFO: No Dirichlet BC found in FemDirichletBC::setup()");

    if (!_is_transient)
        _do_setup = false;
}

}
