/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FeHeatTransportAssembler.h"

#include <iomanip>
#include <boost/any.hpp>

#include "MeshLib/ElementCoordinatesMappingLocal.h"

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"
#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "NumLib/Fem/Tools/LagrangeFeObjectContainer.h"

#include "MaterialLib/FluidModel.h"
#include "MaterialLib/PorousMediumModel.h"
#include "THMCLib/Ogs6FemData.h"


FeHeatTransportAssembler::FeHeatTransportAssembler(const NumLib::IFeObjectContainer &feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
: _feObjects(feObjects), _global_coords(problem_coordinates),
  _e(nullptr), _dof(nullptr), fe(nullptr), /*fe_data(nullptr),*/ _pm_model(nullptr), _solid_model(nullptr), _fluid_model(nullptr)
{
}


void FeHeatTransportAssembler::reset(const MeshLib::Element &e)
{
    _e = &e;
    const size_t mat_id = e.getValue();
    auto* projectData = THMCLib::Ogs6FemData::getInstance();
    _pm_model = projectData->list_pm[mat_id];
    _solid_model = projectData->list_solid[mat_id];
    _fluid_model = projectData->list_fluid[0];
    _sh.resize(e.getDimension(), _global_coords.getDimension(), e.getNNodes());
    fe = _feObjects.getFeObject(*_e);
    auto &mshData = projectData->fe_mesh_data[0];
    THMCLib::updateFeData(mshData, *fe, fe_data);
}

void FeHeatTransportAssembler::reset(const MeshLib::Element &e, const AssemblerLib::LocalToGlobalIndexMap &dof)
{
    reset(e);
    _dof = &dof;
    const size_t n_dof = _dof->dofSize();
    M.resize(n_dof, n_dof);
    K.resize(n_dof, n_dof);
    F.resize(n_dof);
}

void FeHeatTransportAssembler::linear(const NumLib::TimeStep &time,
        const MathLib::LocalVector &T1, const MathLib::LocalVector &T0,
        MathLib::LocalMatrix &localA, MathLib::LocalVector &localRHS)
{
    fe_data.T1 = T1;
    fe_data.T0 = T0;
    //-----------------------------------------
    // Numerical integration
    //-----------------------------------------
    M.setZero();
    K.setZero();
    F.setZero();
    auto& q = fe->getIntegrationMethod();
    for (size_t j=0; j<q.getNPoints(); j++)
    {
        auto wp = q.getWeightedPoint(j);
        fe->computeShapeFunctionsd(wp.getCoords(), _sh);
        const MathLib::LocalMatrix &N_T = _sh.N;
        const MathLib::LocalMatrix &dN_T = _sh.dNdx;
        const MathLib::LocalMatrix &W = N_T;
        const MathLib::LocalMatrix &dW = dN_T;

        //-----------------------------------------
        // Gauss point values
        //-----------------------------------------
        MaterialLib::StateVariables var(THMCLib::getStateVariables(fe_data, _sh));
        auto f = (*_fluid_model)(var);
        auto s = (*_solid_model)(var, _global_coords.getDimension(), *_e);
        auto pm = (*_pm_model)(var, _global_coords.getDimension(), *_e, s, f);
        const MathLib::LocalVector &v = fe_data.v1[j];
        const double fac =  pm.geo_area * _sh.detJ * wp.getWeight();

        //-----------------------------------------
        // Evaluation of matrices
        //-----------------------------------------
        M.noalias() += W * pm.Cp * N_T.transpose() * fac;
        K.noalias() += dW.transpose() * pm.lamba * dN_T * fac;
        K.noalias() += N_T * f.rho * f.cp * v.transpose() * dW * fac;
    }

    const double theta = 1.0;
    localA.noalias() = 1/time.dt() * M + theta * K;
    localRHS.noalias() = (1/time.dt() * M - (1-theta) * K)*fe_data.T0 + F;
}

void FeHeatTransportAssembler::residual(const NumLib::TimeStep &t, const MathLib::LocalVector &p1,
        const MathLib::LocalVector &p0,  MathLib::LocalVector &residual)
{
    const size_t n_dof = _dof->dofSize();
    A.resize(n_dof, n_dof);
    RHS.resize(n_dof);
    linear(t, p1, p0, A, RHS);
    residual.noalias() = A*p1 - RHS;
}

void FeHeatTransportAssembler::jacobian(const NumLib::TimeStep &ts, const MathLib::LocalVector &T1,
        const MathLib::LocalVector &/*T0*/,  MathLib::LocalMatrix &localJ)
{
    fe_data.T1 = T1;
    //-----------------------------------------
    // Numerical integration
    //-----------------------------------------
    M.setZero();
    K.setZero();
    auto& q = fe->getIntegrationMethod();
    for (size_t j=0; j<q.getNPoints(); j++)
    {
        auto wp = q.getWeightedPoint(j);
        fe->computeShapeFunctionsd(wp.getCoords(), _sh);
        const MathLib::LocalMatrix &N_T = _sh.N;
        const MathLib::LocalMatrix &dN_T = _sh.dNdx;
        const MathLib::LocalMatrix &W = N_T;
        const MathLib::LocalMatrix &dW = dN_T;

        //-----------------------------------------
        // Gauss point values
        //-----------------------------------------
        MaterialLib::StateVariables var(THMCLib::getStateVariables(fe_data, _sh));
        auto f = (*_fluid_model)(var);
        auto s = (*_solid_model)(var, _global_coords.getDimension(), *_e);
        auto pm = (*_pm_model)(var, _global_coords.getDimension(), *_e, s, f);
        const MathLib::LocalVector &v = fe_data.v1[j];
        const double fac =  pm.geo_area * _sh.detJ * wp.getWeight();

        //-----------------------------------------
        // Evaluation of matrices
        //-----------------------------------------
        M.noalias() += W * pm.Cp * N_T.transpose() * fac;
        K.noalias() += dW.transpose() * pm.lamba * dN_T * fac;
        K.noalias() += N_T * f.rho * f.cp * v.transpose() * dW * fac;
    }

    const double theta = 1.0;
    localJ.noalias() = 1/ts.dt() * M + theta * K;
}


