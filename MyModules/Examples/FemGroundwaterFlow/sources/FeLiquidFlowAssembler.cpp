/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FeLiquidFlowAssembler.h"

#include <boost/any.hpp>

#include "MeshLib/ElementCoordinatesMappingLocal.h"

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"
#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "NumLib/Fem/Tools/LagrangeFeObjectContainer.h"

#include "MaterialLib/FluidModel.h"
#include "MaterialLib/PorousMediumModel.h"
#include "THMCLib/Ogs6FemData.h"


FeLiquidFlowAssembler::FeLiquidFlowAssembler(const NumLib::IFeObjectContainer &feObjects, const MeshLib::CoordinateSystem &problem_coordinates)
: _feObjects(feObjects), _global_coords(problem_coordinates), _vec_g(MathLib::LocalVector::Zero(_global_coords.getDimension())),
  _e(nullptr), /*_dof(nullptr),*/ fe(nullptr), /*fe_data(nullptr),*/ _pm_model(nullptr), _fluid_model(nullptr)
{
    if (_global_coords.hasZ())
        _vec_g[_global_coords.getIndexOfZ()] = -9.81;
}


void FeLiquidFlowAssembler::reset(const MeshLib::Element &e)
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
    const size_t n_dof = e.getNNodes();
    M.resize(n_dof, n_dof);
    K.resize(n_dof, n_dof);
    F.resize(n_dof);
}

void FeLiquidFlowAssembler::linear(const NumLib::TimeStep &time,
        const MathLib::LocalVector &p1, const MathLib::LocalVector &p0,
        MathLib::LocalMatrix &localA, MathLib::LocalVector &localRHS)
{
    fe_data.p1 = p1;
    fe_data.p0 = p0;
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
        const MathLib::LocalMatrix &N_p = _sh.N;
        const MathLib::LocalMatrix &dN_p = _sh.dNdx;
        const MathLib::LocalMatrix &W = N_p;
        const MathLib::LocalMatrix &dW = dN_p;

        //-----------------------------------------
        // Gauss point values
        //-----------------------------------------
        MaterialLib::StateVariables var(THMCLib::getStateVariables(fe_data, _sh));
        auto s = (*_solid_model)(var, _global_coords.getDimension(), *_e);
        auto f = (*_fluid_model)(var);
        auto pm = (*_pm_model)(var, _global_coords.getDimension(), *_e, s, f);
        const double fac =  pm.geo_area * _sh.detJ * wp.getWeight();

        //-----------------------------------------
        // Evaluation of matrices
        //-----------------------------------------
        M.noalias() += W * pm.Ss * N_p.transpose() * fac;
        K.noalias() += dW.transpose() * pm.k / f.mu * dN_p * fac;
        F.noalias() += dW.transpose() * pm.k / f.mu * f.rho * _vec_g * fac;
    }
//    std::cout << "# Element " << _e->getID() << std::endl;
//    std::cout << "M=\n" << M << std::endl;
//    std::cout << "K=\n" << K << std::endl;
//    std::cout << "F=\n" << F << std::endl;

    const double theta = 1.0;
    localA.noalias() = 1/time.dt() * M + theta * K;
    localRHS.noalias() = (1/time.dt() * M - (1-theta) * K)*fe_data.p0 + F;
}

void FeLiquidFlowAssembler::residual(const NumLib::TimeStep &t, const MathLib::LocalVector &p1,
        const MathLib::LocalVector &p0,  MathLib::LocalVector &residual)
{
    const size_t n_dof = _e->getNNodes();
    A.resize(n_dof, n_dof);
    RHS.resize(n_dof);
    linear(t, p1, p0, A, RHS);
    residual.noalias() = A*p1 - RHS;
    //std::stringstream ss;
    //ss << "# Element " << _e->getID() << std::endl;
    //ss << "p1 = " << p1.transpose() << std::endl;
    //ss << "p0 = " << p0.transpose() << std::endl;
    //ss << "r = " << residual.transpose() << std::endl;
    //INFO("%s", ss.str().data());
}

void FeLiquidFlowAssembler::jacobian(const NumLib::TimeStep &ts, const MathLib::LocalVector &p1,
        const MathLib::LocalVector &/*p0*/,  MathLib::LocalMatrix &localJ)
{
    fe_data.p1 = p1;
    //-----------------------------------------
    // Numerical integration
    //-----------------------------------------
    M.setZero();
    K.setZero();
    auto &q = fe->getIntegrationMethod();
    for (size_t j=0; j<q.getNPoints(); j++)
    {
        auto wp = q.getWeightedPoint(j);
        fe->computeShapeFunctionsd(wp.getCoords(), _sh);
        const MathLib::LocalMatrix &N_p = _sh.N;
        const MathLib::LocalMatrix &dN_p = _sh.dNdx;
        const MathLib::LocalMatrix &W = N_p;
        const MathLib::LocalMatrix &dW = dN_p;

        //-----------------------------------------
        // Gauss point values
        //-----------------------------------------
        MaterialLib::StateVariables var(THMCLib::getStateVariables(fe_data, _sh));
        auto s = (*_solid_model)(var, _global_coords.getDimension(), *_e);
        auto f = (*_fluid_model)(var);
        auto pm = (*_pm_model)(var, _global_coords.getDimension(), *_e, s, f);

        //-----------------------------------------
        // Evaluation of FEM equations
        //-----------------------------------------
        const double fac =  pm.geo_area * _sh.detJ * wp.getWeight();
        M.noalias() += W * pm.Ss * N_p.transpose() * fac;
        K.noalias() += dW.transpose() * pm.k / f.mu * dN_p * fac;
    }

    double theta = 1.0;
    localJ.noalias() = 1./ts.dt() * M + theta * K;
}

void FeLiquidFlowAssembler::velocity(const MathLib::LocalVector &p, std::vector<MathLib::LocalVector> &gp_v)
{
    auto &q = fe->getIntegrationMethod();
    gp_v.resize(q.getNPoints());
    for (size_t ip=0; ip<q.getNPoints(); ip++)
    {
        auto wp = q.getWeightedPoint(ip);
        fe->computeShapeFunctionsd(wp.getCoords(), _sh);
        const MathLib::LocalMatrix &dN_p = _sh.dNdx;

        MaterialLib::StateVariables var(THMCLib::getStateVariables(fe_data, _sh));
        auto s = (*_solid_model)(var, _global_coords.getDimension(), *_e);
        auto f = (*_fluid_model)(var);
        auto pm = (*_pm_model)(var, _global_coords.getDimension(), *_e, s, f);

        gp_v[ip].noalias() = - pm.k / f.mu * (dN_p * p - f.rho * _vec_g);
    }

}

