/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/Function/FemNodalFunction.h"

#include "IElementAssembler.h"
#include "FemVariable.h"

namespace SolutionLib
{

/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
class TransientFEMResidualFunction
    : public NumLib::TemplateFunction<MathLib::IVector, MathLib::IVector>
{
public:

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TransientFEMResidualFunction(
    		MeshLib::Mesh* msh, const std::vector<FemVariable*> &list_var,
			AssemblerLib::LocalToGlobalIndexMap* dofManager, IElementAssembler* asssembler)
        : _msh(msh), _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0), _list_var(list_var)
    {
    }

    ///
    virtual ~TransientFEMResidualFunction() {}

    ///
    NumLib::TemplateFunction<MathLib::IVector,MathLib::IVector>* clone() const
    {
        return new TransientFEMResidualFunction
                    (_msh, _list_var, _dofManager, _local_assembler);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const MathLib::IVector* u_n0, MathLib::IVector* st)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
        this->_st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void operator()(const MathLib::IVector &u_n1, MathLib::IVector &r);

	const NumLib::TimeStep* getTimeStepObj(void)
	{
		return _t_n1; 	
	}

private:
	MeshLib::Mesh* _msh;
	IElementAssembler *_local_assembler;
    AssemblerLib::LocalToGlobalIndexMap* _dofManager;
    const NumLib::TimeStep* _t_n1;
    const MathLib::IVector* _u_n0;
    MathLib::IVector* _st;
    std::vector<FemVariable*> _list_var;
};


} //end
