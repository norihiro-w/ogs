/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementVelocity.hpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

//#include "ElementVelocity.h"

#include "logog.hpp"
#include "MathLib/Vector.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "MaterialLib/PorousMedia.h"

#include "Ogs6FemData.h"

template <class T>
bool FunctionPressureToElementVelocity<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOption<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);

    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);
    _dis = dis;
    _vel = new MyIntegrationPointFunctionVector();
    _vel->initialize(dis);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Velocity), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Velocity, _vel);

    return true;
}

template <class T>
void FunctionPressureToElementVelocity<T>::accept(const NumLib::TimeStep &/*time*/)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Velocity), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);
};

template <class T>
int FunctionPressureToElementVelocity<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Calculating Darcy velocity within elements from fluid pressure...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    MyNodalFunctionScalar *f_p = (MyNodalFunctionScalar*)getInput(Pressure);
    MyIntegrationPointFunctionVector *vel = _vel;;

    MaterialLib::Fluid* fluidphase = Ogs6FemData::getInstance()->list_fluid[0];

    FemLib::LagrangianFeObjectContainer* feObjects = _feObjects;
    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElemenet(i_e);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());
        size_t mat_id = e->getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        // fluid
        double mu = .0;
        fluidphase->dynamic_viscosity->eval(e_pos, mu);
        double rho_f = .0;
        fluidphase->density->eval(e_pos, rho_f);
        // media
        double k;
        pm->permeability->eval(e_pos, k);
        double k_mu;
        k_mu = k / mu;

        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        NumLib::LocalVector local_p(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_p[j] = f_p->getValue(e->getNodeID(j));
        // for each integration points
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        double r[3] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        vel->setNumberOfIntegationPoints(i_e, n_gp);
        NumLib::LocalVector xi(e->getNumberOfNodes());
        NumLib::LocalVector yi(e->getNumberOfNodes());
        NumLib::LocalVector zi(e->getNumberOfNodes());
        for (size_t i=0; i<e->getNumberOfNodes(); i++) {
            const GeoLib::Point* pt = msh->getNodeCoordinatesRef(e->getNodeID(i));
            xi[i] = (*pt)[0];
            yi[i] = (*pt)[1];
            zi[i] = (*pt)[2];
        }
        NumLib::LocalVector q(3);
        for (size_t ip=0; ip<n_gp; ip++) {
            q *= .0;
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            const NumLib::LocalMatrix* dN = fe->getGradBasisFunction();
            NumLib::LocalMatrix* N = fe->getBasisFunction();
            std::vector<double> xx(3, .0);
            NumLib::LocalVector tmp_v;
            tmp_v = (*N) * xi;
            xx[0] = tmp_v[0];
            tmp_v = (*N) * yi;
            xx[1] = tmp_v[0];
            tmp_v = (*N) * zi;
            xx[2] = tmp_v[0];
            NumLib::TXPosition pos(&xx[0]);

            static_cast<NumLib::LocalVector>(q.head(msh->getDimension())) = (*dN) * local_p * (-1.0) * k_mu;
            //TODO grav

            vel->setIntegrationPointValue(i_e, ip, q);
        }
    }
    setOutput(Velocity, vel);
    return 0;
}

