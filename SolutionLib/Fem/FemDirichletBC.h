/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDirichletBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <iostream>
#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Fem/PolynomialOrder.h"

namespace SolutionLib
{

/**
 * \brief DirichletBC data for FEM
 *
 * - mesh
 * - geometry
 * - BC function
 */
class FemDirichletBC
{
public:
    /**
     *
     * @param msh       Pointer to a mesh object
     * @param geo       Pointer to a geometric object representing locations of a boundary condition
     * @param bc_func   Pointer to a temporal-spatial function giving boundary values
     */
    FemDirichletBC(MeshGeoToolsLib::MeshNodeSearcher* mshNodeSearcher, const GeoLib::GeoObject* geo, NumLib::ITXFunction* bc_func);

    /**
     * 
     * @param vec_node_id
     * @param vec_node_values
     */
    FemDirichletBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values);

    ///
    virtual ~FemDirichletBC() {}

    /// setup B.C.
    void setup(NumLib::PolynomialOrder order);

    ///
    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;}

    ///
    std::vector<double>& getListOfBCValues() {return _vec_values;}

    ///
    bool isTransient() const {return _is_transient;}

    void write(std::ostream &os = std::cout) const;

private:
    MeshGeoToolsLib::MeshNodeSearcher* _mshNodeSearcher;
    const GeoLib::GeoObject* _geo;
    NumLib::ITXFunction* _bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;
};

inline std::ostream& operator<< (std::ostream &os, const FemDirichletBC &p)
{
    p.write (os);
    return os;
}

}