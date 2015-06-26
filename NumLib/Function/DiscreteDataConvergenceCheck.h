/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <string>

#include "NumLib/Coupling/Algorithm/IConvergenceCheck.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

namespace NumLib
{

/**
 * \brief Convergence check for iterative calculation of discrete data
 *
 */
class DiscreteDataConvergenceCheck : public NumLib::IConvergenceCheck
{
public:
    ///
    DiscreteDataConvergenceCheck() { };

    ///
    virtual ~DiscreteDataConvergenceCheck() {};

    ///
    virtual bool isConverged(
            NumLib::UnnamedParameterSet& vars_prev,
            NumLib::UnnamedParameterSet& vars_current,
            double eps, double &v_diff)
    {
        for (size_t i=0; i<vars_prev.size(); i++) {
            const NumLib::ITXDiscreteFunction* f_fem_prev = vars_prev.get<NumLib::ITXDiscreteFunction>(i);
            const NumLib::ITXDiscreteFunction* f_fem_cur = vars_current.get<NumLib::ITXDiscreteFunction>(i);
            v_diff = f_fem_prev->diff_norm(*f_fem_cur, MathLib::VecNormType::INFINITY_N);

            if (v_diff>eps) {
                return false;
            }
        }

        return true;
    }
};

}
