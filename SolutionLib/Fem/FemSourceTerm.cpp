/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FemSourceTerm.h"

#include "MeshLib/Mesh.h"

namespace SolutionLib
{

FemSourceTerm::FemSourceTerm(MeshGeoToolsLib::MeshNodeSearcher* mshNodeSearcher, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
: _mshNodeSearcher(mshNodeSearcher), _geo(geo), _bc_func(func->clone()), _t(.0)
{
    _is_transient = !_bc_func->isTemporallyConst();
    _do_setup = true;
}

FemSourceTerm::FemSourceTerm(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
: _mshNodeSearcher(nullptr), _geo(nullptr), _bc_func(nullptr), _vec_nodes(vec_node_id), _vec_values(vec_node_values), _t(.0)
{
    _is_transient = false;
    _do_setup = false;
}

FemSourceTerm::FemSourceTerm(const FemSourceTerm &src)
: _mshNodeSearcher(src._mshNodeSearcher), _geo(src._geo), _bc_func(src._bc_func->clone()), _vec_nodes(src._vec_nodes),
  _vec_values(src._vec_values), _t(src._t), _is_transient(src._is_transient),
  _do_setup(src._do_setup)
{
}

FemSourceTerm::~FemSourceTerm()
{
    delete _bc_func;
}

FemSourceTerm* FemSourceTerm::clone() const
{
    return new FemSourceTerm(*this);
}

void FemSourceTerm::initCurrentTime(double t)
{
    if (_is_transient && _t != t) {
        _t = t;
        _do_setup = true;
    }
}

/// setup
void FemSourceTerm::setup(NumLib::PolynomialOrder order)
{
    if (!_do_setup) return;
    _do_setup = false;

    _vec_nodes.clear();
    auto* msh = _mshNodeSearcher->getMesh();
    //_msh->setCurrentOrder(order);
    _vec_nodes = _mshNodeSearcher->getMeshNodeIDs(*_geo);
    _vec_values.resize(_vec_nodes.size());
    for (size_t i=0; i<_vec_nodes.size(); i++) {
        NumLib::TXPosition pos(_t, msh->getNode(_vec_nodes[i])->getCoords());
        _bc_func->eval(pos, _vec_values[i]);
    }

}

}
