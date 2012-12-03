/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemVariableBuilderXi.cpp
 *
 * Created on 2012-11-16 by Haibing Shao
 */

#include "FemVariableBuilderXi.h"

#include <vector>
#include "logog.hpp"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "SolutionLib/Fem/IFemNeumannBC.h"
#include "SolutionLib/Fem/FemSourceTerm.h"

void FemVariableBuilderXi::doit(const std::string &given_var_name, const BaseLib::Options &option, const MeshLib::IMesh* msh, const GeoLib::GEOObjects *geo, const std::string &geo_unique_name, FemLib::LagrangianFeObjectContainer* _feObjects, SolutionLib::FemVariable* var) const
{
    NumLib::TXFunctionBuilder f_builder;
    // IC
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(msh);
    const BaseLib::Options* opICList = option.getSubGroup("ICList");
    std::vector<const BaseLib::Options*> vec_opIC = opICList->getSubGroupList("IC");
    for (size_t i=0; i<vec_opIC.size(); i++)
    {
        const BaseLib::Options* opIC = vec_opIC[i];
        std::string var_name = opIC->getOption("Variable");
        if (var_name.compare(given_var_name)!=0) continue;
        std::string geo_type = opIC->getOption("GeometryType");
        std::string geo_name = opIC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
        std::string dis_name = opIC->getOption("DistributionType");
        double dis_v = opIC->getOptionAsNum<double>("DistributionValue");
        NumLib::ITXFunction* f_ic =  f_builder.create(*opIC);
        var_ic->addDistribution(geo_obj, f_ic);
    }
    var->setIC(var_ic);
    // BC
    const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    std::vector<const BaseLib::Options*> vec_opBC = opBCList->getSubGroupList("BC");
    for (size_t i=0; i<vec_opBC.size(); i++)
    {
        const BaseLib::Options* opBC = vec_opBC[i];
        std::string var_name = opBC->getOption("Variable");
        if (var_name.compare(given_var_name)!=0) continue;
        std::string geo_type = opBC->getOption("GeometryType");
        std::string geo_name = opBC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
        std::string dis_name = opBC->getOption("DistributionType");
        double dis_v = opBC->getOptionAsNum<double>("DistributionValue");
        NumLib::ITXFunction* f_bc =  f_builder.create(*opBC);
        var->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
    }
    // ST
    const BaseLib::Options* opSTList = option.getSubGroup("STList");
    std::vector<const BaseLib::Options*> vec_opST = opSTList->getSubGroupList("ST");
    for (size_t i=0; i<vec_opST.size(); i++)
    {
        const BaseLib::Options* opST = vec_opST[i];
        std::string var_name = opST->getOption("Variable");
        if (var_name.compare(given_var_name)!=0) continue;
        std::string geo_type = opST->getOption("GeometryType");
        std::string geo_name = opST->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = geo->searchGeoByName(geo_unique_name, geo_type, geo_name);
        std::string st_type = opST->getOption("STType");
        std::string dis_name = opST->getOption("DistributionType");
        double dis_v = opST->getOptionAsNum<double>("DistributionValue");
        if (st_type.compare("NEUMANN")==0) {
            dis_v *= -1; // user set inflow as positive sign but internally negative
        }
        NumLib::ITXFunction* f_st =  f_builder.create(*opST);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            var->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }
}
