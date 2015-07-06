/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFemNeumannBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/Fem/PolynomialOrder.h"

namespace SolutionLib
{

/**
 * \brief Common interface of Neumann BC and ST
 */
class IFemNeumannBC
{
public:

    ///
    virtual ~IFemNeumannBC() {};

    /// setup B.C.
    /// \param order Polynomial order
    virtual void setup(NumLib::PolynomialOrder order) = 0;

    /**
     * set current time
     * @param t
     */
    virtual void initCurrentTime(double t) = 0;

    /// get a list of boundary condition nodes
    virtual std::vector<size_t>& getListOfBCNodes() = 0;

    /// get a list of boundary condition values
    virtual std::vector<double>& getListOfBCValues() = 0;

protected:
    IFemNeumannBC() {};
};

}
