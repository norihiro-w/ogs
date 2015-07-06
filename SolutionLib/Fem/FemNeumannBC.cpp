/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNeumannBC.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemNeumannBC.h"

#include "logog/include/logog.hpp"

#include "NumLib/Fem/BC/NeumannBC2FEM.h"

namespace SolutionLib
{

FemNeumannBC::FemNeumannBC(const MeshLib::Mesh *msh, NumLib::IFeObjectContainer* feObjects, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
: _msh(msh), _feObjects(feObjects), _geo(geo), _bc_func(func->clone()), _t(.0)
{
    _is_transient = !_bc_func->isTemporallyConst();
    _do_setup = true;
}

FemNeumannBC::FemNeumannBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
: _msh(NULL), _feObjects(NULL), _geo(NULL), _bc_func(NULL), _vec_nodes(vec_node_id), _vec_values(vec_node_values), _t(.0)
{
    _is_transient = false;
    _do_setup = false;
}

FemNeumannBC::FemNeumannBC(const FemNeumannBC &src)
: _msh(src._msh), _feObjects(src._feObjects), _geo(src._geo), _bc_func(src._bc_func->clone()),
  _vec_nodes(src._vec_nodes), _vec_values(src._vec_values), _t(src._t),
  _is_transient(src._is_transient), _do_setup(src._do_setup)
{
}

FemNeumannBC::~FemNeumannBC()
{
    delete _bc_func;
}

FemNeumannBC* FemNeumannBC::clone() const
{
    return new FemNeumannBC(*this);
}

void FemNeumannBC::initCurrentTime(double t)
{
    if (_is_transient && _t != t) {
        _t = t;
        _do_setup = true;
    }
}

/// setup BC.
void FemNeumannBC::setup(NumLib::PolynomialOrder order)
{
    if (!_do_setup)
        return;
    _do_setup = false;

    //_msh->setCurrentOrder(order);
    _vec_nodes.clear();
    _vec_values.clear();
    NumLib::NeumannBC2FEM convert(*_msh, _t, *_feObjects, *_geo, *_bc_func, _vec_nodes, _vec_values);
    if (_vec_nodes.size()==0)
        INFO("***INFO: No Neumann BC found in FemDirichletBC::setup()");
}

}
