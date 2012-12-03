/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NRErrorAbsResMNorm.h
 *
 * Created on 2012-11-30 by Haibing Shao
 */

#pragma once

#include <cmath>
#include <algorithm>

#include "logog.hpp"

#include "MathLib/LinAlg/VectorNorms.h"


namespace MathLib
{

/**
 * \brief Check norm-max of residual. 
 */
class NRErrorAbsResMNorm
{
public:
    NRErrorAbsResMNorm() : _itr_count(0) {};

    template<class T_D0>
    inline double error(const T_D0* r, const T_D0* /* dx */, const T_D0* /* x */)
    {
        double abs_mnorm_r = norm_max(*r, r->size());

        return calc(abs_mnorm_r);
    }

private:
    inline double calc(double abs_mnorm_r)
    {
        INFO("-> %d: ||r||_inf=%1.3e ", _itr_count, abs_mnorm_r);

        _itr_count++;

        return abs_mnorm_r;
    }

private:
    size_t _itr_count;
};

}
