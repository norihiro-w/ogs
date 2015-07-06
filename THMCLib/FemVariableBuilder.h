/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemVariableBuilder.h
 *
 * Created on 2012-10-19 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <boost/property_tree/ptree.hpp>

namespace GeoLib
{
class GEOObjects;
}

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class IFeObjectContainer;
}

namespace SolutionLib
{
class FemVariable;

/**
 * \brief FemVariable builder based on BaseLib::Options
 */
class FemVariableBuilder
{
public:
    /**
     * build a variable, i.e. setting IC, BC and ST
     *
     * @param option            Options
     * @param msh               Pointer to a Mesh
     * @param geo               Pointer to geometric objects
     * @param geo_unique_name   Geometry name
     * @param _feObjects        Pointer to FE objects
     * @param var               Pointer to a variable to be configured
     */
    void doit(const std::string &given_var_name, boost::property_tree::ptree const& option,
              const MeshLib::Mesh* msh, const GeoLib::GEOObjects *geo, const std::string &geo_unique_name,
              NumLib::IFeObjectContainer* _feObjects, SolutionLib::FemVariable* var) const;
};

} // SolutionLib
