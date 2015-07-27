/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FemVariableBuilder.h"

#include <vector>
#include <cassert>

#include <logog/include/logog.hpp>

#include "GeoLib/GEOObjects.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

#include "NumLib/Fem/Tools/IFeObjectContainer.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXCompositFunction.h"
#include "NumLib/Function/TXFunctionConstant.h"
#include "NumLib/Function/Multiplication.h"

#include "SolutionLib/Fem/IFemNeumannBC.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "SolutionLib/Fem/FemVariable.h"

namespace SolutionLib
{

void FemVariableBuilder::doit(const std::string &given_var_name,
        boost::property_tree::ptree const& option,
        MeshGeoToolsLib::MeshNodeSearcher* mshNodeSearch,
        MeshGeoToolsLib::BoundaryElementsSearcher* beSearch,
        const GeoLib::GEOObjects *geo,
        const std::string &geo_unique_name, NumLib::IFeObjectContainer* _feObjects,
        SolutionLib::FemVariable* var) const
{
    // IC
    SolutionLib::FemIC* var_ic = new SolutionLib::FemIC(mshNodeSearch);
    auto opICList = option.get_child_optional("ICList");
    if (opICList) {
        auto range = opICList->equal_range("IC");
        for (auto it=range.first; it!=range.second; ++it)
        {
            std::string var_name = it->second.get<std::string>("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = it->second.get<std::string>("GeometryType");
            std::string geo_name = it->second.get<std::string>("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->getGeoObject(geo_unique_name, GeoLib::convertGeoType(geo_type), geo_name);
            //assert(opIC->hasOption("DistributionType"));
            //auto &opDistribution = it->second.get_child("Distribution");
            NumLib::ITXFunction* f_ic = NumLib::TXFunctionBuilder::create(it->second);
            var_ic->addDistribution(geo_obj, f_ic);
        }
    }
    var->setIC(var_ic);

    // BC
    auto opBCList = option.get_child_optional("BCList");
    if (opBCList) {
        auto range = opBCList->equal_range("BC");
        for (auto it=range.first; it!=range.second; ++it)
        {
            std::string var_name = it->second.get<std::string>("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = it->second.get<std::string>("GeometryType");
            std::string geo_name = it->second.get<std::string>("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->getGeoObject(geo_unique_name, GeoLib::convertGeoType(geo_type), geo_name);
//            assert(opBC->hasOption("DistributionType"));
            //auto &opDistribution = it->second.get_child("Distribution");
            NumLib::ITXFunction* f_bc = NumLib::TXFunctionBuilder::create(it->second);
            var->addDirichletBC(new SolutionLib::FemDirichletBC(mshNodeSearch, geo_obj, f_bc));
        }
    }

    // ST
    auto opSTList = option.get_child_optional("STList");
    if (opSTList) {
        auto range = opSTList->equal_range("ST");
        for (auto it=range.first; it!=range.second; ++it)
        {
            std::string var_name = it->second.get<std::string>("Variable");
            if (var_name.compare(given_var_name)!=0) continue;
            std::string geo_type = it->second.get<std::string>("GeometryType");
            std::string geo_name = it->second.get<std::string>("GeometryName");
            const GeoLib::GeoObject* geo_obj = geo->getGeoObject(geo_unique_name, GeoLib::convertGeoType(geo_type), geo_name);
            std::string st_type = it->second.get<std::string>("STType");
            //assert(it->second.hasOption("DistributionType"));
            //auto &opDistribution = it->second.get_child("Distribution");
            NumLib::ITXFunction* f_st = NumLib::TXFunctionBuilder::create(it->second);
            if (st_type.compare("NEUMANN")==0) {
                // user set inflow as positive sign but internally negative
                f_st = new NumLib::TXCompositFunction
                <
                    NumLib::ITXFunction, NumLib::TXFunctionConstant,
                    NumLib::Multiplication
                >(f_st, new NumLib::TXFunctionConstant(-1.));
            }
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(mshNodeSearch, beSearch, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(mshNodeSearch, geo_obj, f_st);
            }
            var->addNeumannBC(femSt);
        }
    }
}

}// SolutionLib
