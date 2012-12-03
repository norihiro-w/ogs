/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NRPrePostDummy.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <cmath>
#include <algorithm>

namespace MathLib
{

/**
 * \brief Dummy class for pre-post processign in Newton iterations
 *
 */
class NRIterationStepInitializerDummy
{
public:
    /// Default setting
    NRIterationStepInitializerDummy() {};

    /**
     * do pre-post processing
     *
     * \param dx        solution increment
     * \param x_new     new solution
     */
	template<class T_X, class F_RESIDUALS, class F_DX>
    void pre_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
    }

	template<class T_X, class F_RESIDUALS, class F_DX>
    void post_process(const T_X &/*dx*/, const T_X &/*x_new*/, F_RESIDUALS &/*f_residuals*/, F_DX &/*f_dx*/)
    {
    }
};


}
