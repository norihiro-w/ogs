
#pragma once

#include "IProblem.h"

#include <vector>

#include "Base/CodingTools.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"


namespace SolutionLib
{

/**
 * \brief IVBV problems using FEM
 *
 *- Variables
 *- Equations (local assembly)
 *- IC
 *- BC
 */
template<class T_ASSEMBLY>
class FemIVBVProblem : public AbstractMeshBasedDiscreteIVBVProblem
{
public:
    FemIVBVProblem(DiscreteLib::DiscreteSystem &dis, MeshLib::IMesh &msh, T_ASSEMBLY &user_assembly)
        : AbstractMeshBasedDiscreteIVBVProblem(msh), _user_assembly(user_assembly), _discrete_system(&dis)
    {
        Base::zeroObject(_map_var, _map_ic);
    }

    virtual ~FemIVBVProblem()
    {
        Base::releaseObject(_map_var, _map_ic);
        Base::releaseObjectsInStdVector(_map_bc1);
        Base::releaseObjectsInStdVector(_map_bc2);
    }

    /// get the number of variables
    size_t getNumberOfVariables() const { return 1; }

    /// create FE approximation field
    size_t createField(FemLib::PolynomialOrder::type order)
    {
        if (_map_var==0) {
            _map_var = new FemLib::FemNodalFunctionScalar(*_discrete_system, *getMesh(), order);
        }
        return 0;
    }

    /// get the FE field
    FemLib::FemNodalFunctionScalar* getField(size_t) const
    {
        return _map_var;
    }

    void setIC(int, SpatialFunction& ic)
    {
        _map_ic = ic.clone();
    }

    SpatialFunction* getIC(int) const
    {
        return _map_ic;
    };

    void addDirichletBC(int, GeoLib::GeoObject &geo, bool is_transient, SpatialFunction& bc1)
    {
        addDirichletBC(*new FemLib::FemDirichletBC<double>(_map_var, &geo, is_transient, &bc1, new FemLib::DiagonalizeMethod()));
    }


    size_t getNumberOfDirichletBC(int n=0) const {return _map_bc1.size();};

    SpatialFunction* getDirichletBC(int, int bc_id) const
    {
        return _map_bc1[bc_id];
    };

    FemLib::FemDirichletBC<double>* getFemDirichletBC(int bc_id) const 
    {
        return _map_bc1[bc_id];
    };

    void addNeumannBC(int, GeoLib::GeoObject &geo, bool is_transient, SpatialFunction& bc2)
    {
        addNeumannBC(*new FemLib::FemNeumannBC<double, double>(_map_var, &geo, is_transient, &bc2));
    }

    size_t getNumberOfNeumannBC(int n=0) const {return _map_bc2.size();};

    SpatialFunction* getNeumannBC(int, int bc_id) const
    {
        return _map_bc2[bc_id];
    };

    FemLib::FemNeumannBC<double, double>* getFemNeumannBC(int bc_id) const 
    {
        return _map_bc2[bc_id];
    };

    T_ASSEMBLY& getElementAssemlby()
    {
        return _user_assembly;
    }

private:
    T_ASSEMBLY _user_assembly;
    DiscreteLib::DiscreteSystem* _discrete_system;
    FemLib::FemNodalFunctionScalar* _map_var;
    SpatialFunction* _map_ic;
    std::vector<FemLib::FemDirichletBC<double>*> _map_bc1;
    std::vector<FemLib::FemNeumannBC<double, double>*> _map_bc2;


    void addDirichletBC(FemLib::FemDirichletBC<double>& bc1)
    {
        _map_bc1.push_back(&bc1);
    }

    void addNeumannBC(FemLib::FemNeumannBC<double, double>& bc2)
    {
        _map_bc2.push_back(&bc2);
    }

    DISALLOW_COPY_AND_ASSIGN(FemIVBVProblem);
};



} //end
