/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <string>

#include <boost/property_tree/ptree.hpp>

#include "NonlinearSolverOption.h"
#include "DiscreteNonlinearSolverFactory.h"

namespace NumLib
{

template <
    class F_LINEAR, 
    class F_R, 
    class F_J,
    class T_NL_FACTORY=DiscreteNonlinearSolverFactory
>
class TemplateDiscreteNonlinearSolver
{
public:
    TemplateDiscreteNonlinearSolver(F_LINEAR* f_l, F_R* f_r, F_J* f_J, MathLib::IMatrix* J, MathLib::IVector* r)
    : _f_l(f_l), _f_r(f_r), _f_J(f_J), _solver(nullptr), _nl_factory(new T_NL_FACTORY), _J(J), _r(r)
    {
    }

    virtual ~TemplateDiscreteNonlinearSolver()
    {
        delete _solver;
        delete _nl_factory;
    }

    void setOption(const boost::property_tree::ptree &option)
    {
        auto op = option.get_child_optional("NonlinearSolver");
        if (!op) return;

        if (auto v = op->get_optional<std::string>("solver_type"))
            _option.solver_type = _option.getSolverType(*v);
        if (auto v = op->get_optional<double>("error_tolerance"))
            _option.error_tolerance = *v;
        if (auto v = op->get_optional<long>("max_iteration_step"))
            _option.max_iteration = *v;
    }

    void setOption(const NonlinerSolverOption &option) { _option = option; }

    NonlinerSolverOption &getOption() const { return _option; }

    void solve(MathLib::IVector &x)
    {
        if (_solver==nullptr)
            _solver = _nl_factory->create(_option, _f_l, _f_r, _f_J, _J, _r);
        _solver->solve(x);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(TemplateDiscreteNonlinearSolver);

private:
    NonlinerSolverOption _option;
    F_LINEAR* _f_l;
    F_R* _f_r;
    F_J* _f_J;
    INonlinearSolver* _solver;
    T_NL_FACTORY* _nl_factory;
    MathLib::IMatrix* _J;
    MathLib::IVector* _r;
};

} //end

