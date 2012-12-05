/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemResidualEQS.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TransientAssembler/TransientElementWiseVectorUpdater.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"

#include "FemVariable.h"
#include "SolutionLib/DataType.h"

namespace SolutionLib
{

/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <
    class T_DIS_SYS,
    class T_LOCAL_RESIDUAL_ASSEMBLER
    >
class TemplateTransientResidualFEMFunction
    : public NumLib::TemplateFunction<SolutionVector, SolutionVector>
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef FemVariable MyVariable;
    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;
    typedef typename NumLib::TransientElementWiseVectorUpdater<UserLocalResidualAssembler> MyUpdater;
    typedef typename T_DIS_SYS::template MyVectorAssembler<double, MyUpdater>::type MyGlobalAssembler;

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TemplateTransientResidualFEMFunction(MyDiscreteSystem* sys, const std::vector<MyVariable*> &list_var, DiscreteLib::DofEquationIdTable* dofManager, UserLocalResidualAssembler* asssembler)
        : _dis_sys(sys), _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0), _list_var(list_var)
    {
    };

    ///
    virtual ~TemplateTransientResidualFEMFunction() {};

    ///
    NumLib::TemplateFunction<SolutionVector,SolutionVector>* clone() const
    {
        return new TemplateTransientResidualFEMFunction
                    <
                        T_DIS_SYS,
                        UserLocalResidualAssembler
                    >(_dis_sys, _list_var, _dofManager, _local_assembler);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const SolutionVector* u_n0, SolutionVector* st)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
        this->_st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void eval(const SolutionVector &u_n1, SolutionVector &r);

	const NumLib::TimeStep* getTimeStepObj(void)
	{
		return _t_n1; 	
	}

private:
    MyDiscreteSystem *_dis_sys;
    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    const NumLib::TimeStep* _t_n1;
    const SolutionVector* _u_n0;
    SolutionVector* _st;
    std::vector<MyVariable*> _list_var;
};


template <class T1, class T2>
void TemplateTransientResidualFEMFunction<T1,T2>::eval(const SolutionVector &u_n1, SolutionVector &r)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const SolutionVector *u_n = this->_u_n0;
    size_t msh_id = _dis_sys->getMesh()->getID();

    // assembly
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    MyUpdater updater(&t_n1, msh, _dofManager, u_n, &u_n1, _local_assembler);
    MyGlobalAssembler assembler(&updater);
    r = .0;
    r.construct(*_dofManager, assembler);


    // add source/sink term (i.e. RHS in linear equation)
    if (_st!=0) {
        r += *_st;
    }

    // set residuals to zero for Dirichlet BC
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            DiscreteLib::convertToEqsValues(*_dofManager, i, msh_id, bc1->getListOfBCNodes(), bc1->getListOfBCValues(), list_bc1_eqs_id, list_bc1_val);
        }
    }
    r.setSubvector(list_bc1_eqs_id, .0);
}
} //end
