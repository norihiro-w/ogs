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
bool FunctionHeadToElementVelocity<T>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    size_t tim_id = option.getOptionAsNum<size_t>("TimeGroupID");
    _tim = femData->list_tim[tim_id];
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);
    _dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _vel = new MyIntegrationPointFunctionVector();

    // initialize integration point velocity
    _vel->initialize(_dis);
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        const MeshLib::IElement* e = msh->getElement(i_e);
        FemLib::IFiniteElement *fe = _feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        _vel->setNumberOfIntegationPoints(i_e, n_gp);
        MathLib::LocalVector q = MathLib::LocalVector::Zero(3);
        for (size_t ip=0; ip<n_gp; ip++) {
            _vel->setIntegrationPointValue(i_e, ip, q);
        }
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Velocity), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    _vel_3d = new My3DIntegrationPointFunctionVector(_vel, msh->getGeometricProperty()->getCoordinateSystem());
    this->setOutput(Velocity, _vel_3d); //_vel_3d

    return true;
}

template <class T>
void FunctionHeadToElementVelocity<T>::accept(const NumLib::TimeStep &time)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Velocity), _dis->getMesh()->getID(), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);

    _tim->accept(time.getTime());
};

template <class T>
int FunctionHeadToElementVelocity<T>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Calculating Darcy velocity within elements from hydraulic head...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    const MeshLib::CoordinateSystem coord = msh->getGeometricProperty()->getCoordinateSystem();
    FemLib::LagrangeFeObjectContainer* feObjects = _feObjects;
    const MyNodalFunctionScalar* head = (MyNodalFunctionScalar*)getInput(Head);
    MyIntegrationPointFunctionVector *vel = _vel;;

    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        //____________________________________________
        // collect element information
        const MeshLib::IElement* e = msh->getElement(i_e);
        const MeshLib::ElementCoordinatesMappingLocal* ele_local_coord
            = (MeshLib::ElementCoordinatesMappingLocal*)e->getMappedCoordinates();
        const MathLib::LocalMatrix &matR = ele_local_coord->getRotationMatrixToOriginal();
        const size_t mat_id = e->getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        MathLib::LocalVector local_h(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_h[j] = head->getValue(e->getNodeID(j));

        //____________________________________________
        // define FE object
        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();

        //____________________________________________
        // for each integration points
        double r[3] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        vel->setNumberOfIntegationPoints(i_e, n_gp);
        MathLib::LocalVector q = MathLib::LocalVector::Zero(3);
        double xx[3] = {.0};
        for (size_t ip=0; ip<n_gp; ip++) {
            // calculate shape functions etc at this integration point
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            const MathLib::LocalMatrix* dN = fe->getGradBasisFunction();
            //MathLib::LocalMatrix* N = fe->getBasisFunction();
            fe->getRealCoordinates(xx);

            // calculate material parameter at this integration point
            NumLib::TXPosition pos(NumLib::TXPosition::Element, e->getID(), &xx[0]);
            double k;
            pm->hydraulic_conductivity->eval(pos, k);
            MathLib::LocalMatrix local_k = MathLib::LocalMatrix::Identity(e->getDimension(), e->getDimension());
            local_k *= k;
            MathLib::LocalMatrix global_k;
            if (e->getDimension() < coord.getDimension()) {
                MathLib::LocalMatrix local2 = MathLib::LocalMatrix::Zero(coord.getDimension(), coord.getDimension());
    //            local2.topLeftCorner(local_k_mu.rows(), local_k_mu.cols()) = local_k_mu;
                local2.block(0, 0, local_k.rows(), local_k.cols()) = local_k.block(0, 0, local_k.rows(), local_k.cols());
                global_k = matR * local2 * matR.transpose();
            } else {
                global_k = local_k;
            }

            // calculate Darcy velocity
            q = - global_k * (*dN) * local_h;

            // update array
            vel->setIntegrationPointValue(i_e, ip, q);
        }
    }

    _vel_3d->resetVectorFunction(vel);
    setOutput(Velocity, _vel_3d);

    return 0;
}

