/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Fem/Tools/IFeObjectContainer.h"

#include "IFemNeumannBC.h"


namespace SolutionLib
{

/**
 * \brief Neumann BC class for a variable
 *
 */
class FemNeumannBC : public IFemNeumannBC
{
public:
    /// 
    FemNeumannBC(const MeshLib::Mesh *msh, NumLib::IFeObjectContainer* feObjects, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func);

    ///
    FemNeumannBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values);

    /**
     * Copy constructor
     * @param src
     */
    FemNeumannBC(const FemNeumannBC &src);

    ///
    virtual ~FemNeumannBC();

    /// clone this object
    virtual FemNeumannBC* clone() const;

    /// setup B.C.
    /// \param order Polynomial order
    virtual void setup(NumLib::PolynomialOrder order);

    /**
     * set current time
     * @param t
     */
    virtual void initCurrentTime(double t);

    /// get a list of boundary condition nodes
    virtual std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};

    /// get a list of boundary condition values
    virtual std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    const MeshLib::Mesh* _msh;
    NumLib::IFeObjectContainer* _feObjects;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    double _t;
    bool _is_transient;
    bool _do_setup;

};

}
