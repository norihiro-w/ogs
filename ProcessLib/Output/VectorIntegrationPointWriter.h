/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

#include "IntegrationPointWriter.h"

namespace ProcessLib
{

struct VectorIntegrationPointWriter final : public IntegrationPointWriter
{
    explicit VectorIntegrationPointWriter(
        std::string const& name,
        int const n_components,
        int const integration_order,
        std::function<std::vector<std::vector<double>>()>
            callback)
        : _name(name),
          _n_components(n_components),
          _integration_order(integration_order),
          _callback(callback)
    {
    }

    int numberOfComponents() const override { return _n_components; }
    int integrationOrder() const override { return _integration_order; }

    std::string name() const override
    {
        // TODO (naumov) remove ip suffix. Probably needs modification of the
        // mesh properties, s.t. there is no "overlapping" with cell/point data.
        // See getOrCreateMeshProperty.
        return _name;
    }

    std::vector<std::vector<double>> values() const override
    {
        return _callback();
    }

private:
    std::string const _name;
    int const _n_components;
    int const _integration_order;
    std::function<std::vector<std::vector<double>>()> _callback;
};

}  // namespace ProcessLib
