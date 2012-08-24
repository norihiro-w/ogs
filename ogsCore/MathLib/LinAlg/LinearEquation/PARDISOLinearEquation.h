/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PARDISOLinearEquation.h
 *
 * Created on 2012-08-23 by Norihiro Watanabe
 */

#ifdef USE_MKL_PARDISO

#pragma once

#include "AbstractCRSLinearEquation.h"

namespace MathLib
{

class PARDISOLinearEquation : public AbstractCRSLinearEquation<signed>
{
public:
    virtual ~PARDISOLinearEquation();

    void setOption(const BaseLib::Options &/*option*/) {};

protected:
    void solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x);

};

}

#endif

