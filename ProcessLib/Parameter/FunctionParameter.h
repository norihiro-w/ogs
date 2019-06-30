/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>
#include <vector>

#include <exprtk.hpp>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

#include "Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{

/// A parameter class evaluating a functon defined by
/// a user-provided maethmatical expression.
///
/// Currently, x, y, and z are supported as variables
/// of the function.
template <typename T>
struct FunctionParameter final : public Parameter<T>
{
    typedef exprtk::symbol_table<T> symbol_table_t;
    typedef exprtk::expression<T> expression_t;
    typedef exprtk::parser<T> parser_t;
    typedef exprtk::parser_error::type error_t;

    /**
     * Constructing from a vector of expressions
     *
     * @param name_       the parameter's name
     * @param mesh_       a mesh object
     * @param vec_expression_str_  a vector of mathematical expressions
     * The vector size specifies the number of components of the parameter.
     */
    FunctionParameter(std::string const& name_,
                         MeshLib::Mesh const& mesh_,
                         std::vector<std::string> const& vec_expression_str_)
        : Parameter<T>(name_), _mesh(mesh_), _vec_expression_str(vec_expression_str_)
    {
        _symbol_table.add_constants();
        _symbol_table.create_variable("x");
        _symbol_table.create_variable("y");
        _symbol_table.create_variable("z");

        _vec_expression.resize(_vec_expression_str.size());
        for (unsigned i=0; i<_vec_expression_str.size(); i++)
        {
            _vec_expression[i].register_symbol_table(_symbol_table);
            parser_t parser;
            if (!parser.compile(_vec_expression_str[i], _vec_expression[i]))
            {
                ERR("Error: %s\tExpression: %s\n", parser.error().c_str(),
                    _vec_expression_str[i].c_str());
            }
        }
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfComponents() const override
    {
        return _vec_expression.size();
    }

    std::vector<T> const operator()(double const /*t*/,
                                     SpatialPosition const& pos) const override
    {
        std::vector<T> cache(getNumberOfComponents());
        auto& x = _symbol_table.get_variable("x")->ref();
        auto& y = _symbol_table.get_variable("y")->ref();
        auto& z = _symbol_table.get_variable("z")->ref();
        if (pos.getCoordinates())
        {
            auto const coords = pos.getCoordinates().get();
            x = coords[0];
            y = coords[1];
            z = coords[2];
        }
        else if (pos.getNodeID())
        {
            auto const& node = *_mesh.getNode(pos.getNodeID().get());
            x = node[0];
            y = node[1];
            z = node[2];
        }

        for (unsigned i=0; i<_vec_expression.size(); i++)
            cache[i] = _vec_expression[i].value();

        return cache;
    }

private:
    MeshLib::Mesh const& _mesh;
    std::vector<std::string> _vec_expression_str;

    symbol_table_t _symbol_table;
    std::vector<expression_t> _vec_expression;
};

std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // ProcessLib
