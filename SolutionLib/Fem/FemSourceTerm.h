/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemSourceTerm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
//#include "MeshLib/Tools/Tools.h"

#include "NumLib/Function/ITXFunction.h"

#include "IFemNeumannBC.h"

namespace SolutionLib
{

/**
 * 
 */
class FemSourceTerm : public IFemNeumannBC 
{
public:
    /// 
    FemSourceTerm(const MeshLib::Mesh *msh, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func);

    ///
    FemSourceTerm(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values);

    ///
    FemSourceTerm(const FemSourceTerm &src);

    ///
    virtual ~FemSourceTerm();

    /// clone this object
    FemSourceTerm* clone() const;

    virtual void initCurrentTime(double t);

    /// setup  
    virtual void setup(NumLib::PolynomialOrder order);

    /// get a list of boundary condition nodes
    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};

    /// get a list of boundary condition values
    std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    // node id, var id, value
    const MeshLib::Mesh* _msh;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    double _t;
    bool _is_transient;
    bool _do_setup;
};

}
