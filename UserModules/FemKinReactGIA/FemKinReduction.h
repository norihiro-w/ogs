/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemKinReduction.h
 *
 * Created on 2012-09-20 by Haibing Shao
 */

#ifndef FEM_KIN_REDUCTION_H
#define FEM_KIN_REDUCTION_H

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "SolutionLib/Core/TimeSteppingProblem.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "ChemLib/chemReductionKin.h"

namespace SolutionLib
{

/**
 * \brief IVBV problems for FEM
 *
 * This class contains
 * - Variables
 * - Equations
 * - Reference to discrete system
 *
 * \tparam T_FEM_EQUATION   FEM equation
 */
template < class T_DIS_SYS >
class FemKinReduction : public TimeSteppingProblem
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef FemVariable MyVariable;

    /**
      * constructor
      */
    explicit FemKinReduction(MyDiscreteSystem* dis, ogsChem::chemReductionKin* myReductionKin )
     : _discrete_system(dis), _ReductionKin(myReductionKin)
    {
    }

    /**
      * destructor
      */
    virtual ~FemKinReduction()
    {
        BaseLib::releaseObjectsInStdVector(_variables);
        _ReductionKin = NULL;
        _discrete_system = NULL;
    }

    /**
      * get this discrete system
      */
    MyDiscreteSystem* getDiscreteSystem() {return _discrete_system;};

    /**
      * create FE approximation field
      */
    MyVariable* addVariable(const std::string &name)
    {
        _variables.push_back(new MyVariable(_variables.size(), name));
        return _variables.back();
    }

    /**
      * get a variable
      */
    MyVariable* getVariable(size_t i) const { return _variables[i]; }

    /**
      * get the number of variables
      */
    size_t getNumberOfVariables() const { return _variables.size(); }

	/**
      * get the pointer to reduction scheme
      */
	ogsChem::chemReductionKin* getReductionScheme() const { return _ReductionKin; }

private:
    DISALLOW_COPY_AND_ASSIGN(FemKinReduction);

private:
    /**
      * discretization 
      */ 
    MyDiscreteSystem* _discrete_system;

    /**
      * vector of variables
      */ 
    std::vector<MyVariable*> _variables;

    /**
      * pointer to reduction scheme
      */ 
	ogsChem::chemReductionKin* _ReductionKin;

};

} //end

#endif  // end of ifndef
