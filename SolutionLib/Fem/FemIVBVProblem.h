/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/CodingTools.h"
#include "NumLib/Fem/Tools/LagrangeFeObjectContainer.h"
#include "SolutionLib/Core/TimeSteppingProblem.h"

#include "IElementAssembler.h"
#include "FemVariable.h"

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
 */
class FemIVBVProblem : public TimeSteppingProblem
{
public:
    typedef FemVariable MyVariable;

    ///
    FemIVBVProblem(const MeshLib::Mesh &msh, IElementAssembler* assembler)
     : _msh(msh), _assembler(assembler)
    {
    }

    ///
    virtual ~FemIVBVProblem()
    {
        for (auto p : _variables) delete p;
        delete _assembler;
    }

    /// create FE approximation field
    MyVariable* addVariable(const std::string &name, NumLib::PolynomialOrder initial_order = NumLib::PolynomialOrder::Linear)
    {
        MyVariable* var = new MyVariable(_variables.size(), name, initial_order);
        NumLib::LagrangeFeObjectContainer* feContainer = new NumLib::LagrangeFeObjectContainer(&_msh);
        var->setFeObjectContainer(feContainer);
        _variables.push_back(var);
        return var;
    }

    /// get a variable
    MyVariable* getVariable(size_t i) const { return _variables[i]; }

    /// get the number of variables
    size_t getNumberOfVariables() const { return _variables.size(); }

    /// get this equation
    IElementAssembler* getAssembler() const {return _assembler;}

    const MeshLib::Mesh& getMesh() const {return _msh;}

private:
    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);

private:
    const MeshLib::Mesh& _msh;
    std::vector<MyVariable*> _variables;
    IElementAssembler* _assembler;
};

} //end
