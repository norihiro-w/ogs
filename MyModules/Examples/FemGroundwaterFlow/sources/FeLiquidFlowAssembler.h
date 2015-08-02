/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/DataType.h"

#include "MeshLib/CoordinateSystem.h"
#include "NumLib/Fem/Tools/IFeObjectContainer.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "SolutionLib/Fem/IElementAssembler.h"
#include "THMCLib/FeElementData.h"

namespace NumLib
{
class LagrangeFeObjectContainer;
class IFiniteElement;
}

namespace MaterialLib
{
struct PorousMediumModel;
struct FluidModel;
struct SolidModel;
}

/**
 * \brief Local assembly class for time-ODE GW equation
 */
class FeLiquidFlowAssembler: public SolutionLib::IElementAssembler
{
public:
    FeLiquidFlowAssembler(const NumLib::IFeObjectContainer &feObjects,
            const MeshLib::CoordinateSystem &problem_coordinates);

    virtual ~FeLiquidFlowAssembler() {}

    void reset(const MeshLib::Element &e) override;

    void linear(const NumLib::TimeStep &/*time*/, const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,
            MathLib::LocalMatrix &localA, MathLib::LocalVector &localRHS) override;

    void residual(const NumLib::TimeStep &/*time*/,
            const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,
            MathLib::LocalVector &localR) override;

    void jacobian(const NumLib::TimeStep &/*time*/,
            const MathLib::LocalVector &/*u1*/, const MathLib::LocalVector &/*u0*/,
            MathLib::LocalMatrix &localJ) override;

    void velocity(const MathLib::LocalVector &p, std::vector<MathLib::LocalVector> &gp_v);

private:
    const NumLib::IFeObjectContainer &_feObjects;
    const MeshLib::CoordinateSystem _global_coords;
    MathLib::LocalVector _vec_g;
    const MeshLib::Element* _e;
    NumLib::IFiniteElement* fe;
    THMCLib::FeElementData fe_data;
    MaterialLib::PorousMediumModel* _pm_model;
    MaterialLib::SolidModel* _solid_model;
    MaterialLib::FluidModel* _fluid_model;
    NumLib::DynamicShapeMatrices _sh;
    MathLib::LocalMatrix M;
    MathLib::LocalMatrix K;
    MathLib::LocalVector F;
    MathLib::LocalMatrix A;
    MathLib::LocalVector RHS;
};

