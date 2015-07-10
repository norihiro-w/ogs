/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NewtonRaphson.h"
#include <limits>

namespace MathLib
{

namespace Nonlinear
{

NewtonRaphson::NewtonRaphson()
: _normType(VecNormType::INFINITY_N),
  _r_abs_tol(std::numeric_limits<double>::max()), _r_rel_tol(1e-6), _dx_rel_tol(.0),
  _max_itr(25), _printErrors(true), _n_iterations(0), _r_abs_error(.0), _r_rel_error(.0), _dx_rel_error(.0),
  _J(nullptr), _r(nullptr), _b(nullptr), _option(nullptr)
{
}

NewtonRaphson::NewtonRaphson(IMatrix* J, IVector* r, boost::property_tree::ptree const*const option)
: _normType(VecNormType::INFINITY_N),
  _r_abs_tol(std::numeric_limits<double>::max()), _r_rel_tol(1e-6), _dx_rel_tol(.0),
  _max_itr(25), _printErrors(true), _n_iterations(0), _r_abs_error(.0), _r_rel_error(.0), _dx_rel_error(.0),
  _J(J), _r(r), _b(r->duplicate()), _option(option)
{
}

NewtonRaphson::~NewtonRaphson()
{
    delete _b;
}

} // Nonlinear
} //end
