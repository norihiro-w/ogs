/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GWTools.h"

#include <logog/include/logog.hpp>

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"
#include "THMCLib/Ogs6FemData.h"
#include "THMCLib/MaterialLib/FluidModel.h"
#include "THMCLib/MaterialLib/PorousMediumModel.h"
#include "FeLiquidFlowAssembler.h"

void getDarcyVelocity(const NumLib::FemNodalFunctionScalar &f_p, NumLib::FEMIntegrationPointFunctionVector &vel)
{
    auto &msh = f_p.getMesh();
//    auto coord = msh.getCoordinateSystem();
    auto feObjects = f_p.getFeObjectContainer();

//    MaterialLib::FluidModel* fluid_model = THMCLib::Ogs6FemData::getInstance()->list_fluid[0];
//    MathLib::LocalVector vec_g = MathLib::LocalVector::Zero(coord.getDimension());
//    if (coord.hasZ())
//        vec_g[coord.getIndexOfZ()] = -9.81;
    FeLiquidFlowAssembler assembler(*feObjects, msh.getCoordinateSystem());

    for (size_t i_e=0; i_e<msh.getNElements(); i_e++)
    {
        auto e = msh.getElement(i_e);
        MathLib::LocalVector local_p(e->getNNodes());
        for (size_t j=0; j<e->getNNodes(); j++)
            local_p[j] = f_p.getValue(e->getNodeIndex(j));

        assembler.reset(*e);
        std::vector<MathLib::LocalVector> gp_vel;
        assembler.velocity(local_p, gp_vel);

        vel.setNumberOfIntegationPoints(i_e, gp_vel.size());
        for (size_t ip=0; ip<gp_vel.size(); ip++) {
            vel.setIntegrationPointValue(i_e, ip, gp_vel[ip]);
        }
    }
}

