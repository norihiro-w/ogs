/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FunctionParameter.h"

#include <unordered_map>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Function");

    std::vector<std::string> vec_expressions;

    //! \ogs_file_param{prj__parameters__parameter__Function__expression}
    for (auto p : config.getConfigSubtreeList("expression"))
    {
        std::string const expression_str = p.getValue<std::string>();
        vec_expressions.emplace_back(expression_str);
    }

    std::unordered_map<std::string,double> map_extra_variables;
    //! \ogs_file_param{prj__parameters__parameter__Function__extra_variable}
    for (auto configExtraVar : config.getConfigSubtreeList("extra_variable"))
    {
        //! \ogs_file_param{prj__parameters__parameter__Function__extra_variable_name}
        auto const var_name = configExtraVar.getConfigParameter<std::string>("name");
        //! \ogs_file_param{prj__parameters__parameter__Function__extra_variable_name_default_value}
        auto const var_default_value = configExtraVar.getConfigParameter<double>("default_value");
        map_extra_variables.emplace(var_name, var_default_value);
    }

    return std::make_unique<FunctionParameter<double>>(
        name, mesh, vec_expressions, map_extra_variables);
}

}  // ProcessLib
