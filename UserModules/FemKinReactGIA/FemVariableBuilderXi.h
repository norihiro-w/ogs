/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemVariableBuilderXi.h
 *
 * Created on 2012-11-16 by Haibing Shao
 */

#ifndef FEM_VARIABLE_BUILDER_XI_H
#define FEM_VARIABLE_BUILDER_XI_H

#include <string>
#include "BaseLib/Options.h"
#include "GeoLib/GEOObjects.h"
#include "MeshLib/Core/IMesh.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "FemVariableBuilder.h"

/**
 * \brief FemVariable builder for xi based on BaseLib::Options
 * overwriting the original doit function
 */
class FemVariableBuilderXi : FemVariableBuilder
{
public:
    /**
     * build the variable eta, i.e. setting IC, BC and ST
     *
     * @param option            Options
     * @param msh               Pointer to a Mesh
     * @param geo               Pointer to geometric objects
     * @param geo_unique_name   Geometry name
     * @param _feObjects        Pointer to FE objects
     * @param var               Pointer to a variable to be configured
     */
    void doit(const std::string &given_var_name, const BaseLib::Options &option, const MeshLib::IMesh* msh, const GeoLib::GEOObjects *geo, const std::string &geo_unique_name, FemLib::LagrangianFeObjectContainer* _feObjects, SolutionLib::FemVariable* var) const;
};

#endif