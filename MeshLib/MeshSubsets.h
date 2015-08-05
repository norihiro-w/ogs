/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSUBSETS_H_
#define MESHSUBSETS_H_

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <vector>

#include <logog/include/logog.hpp>

#include "MeshSubset.h"

namespace MeshLib
{

/// Collection of mesh subsets.
class MeshSubsets
{
public:

    /// Single mesh subset constructor.
    MeshSubsets(const MeshSubset* mesh_subset)
    {
        _mesh_subsets.push_back(mesh_subset);
        _n_total_items = mesh_subset->getNTotalItems();
        _n_total_nonghost_items = mesh_subset->getNTotalNonGhostItems();
    }

    /// Construct MeshSubsets from a range of MeshSubset. InputIterator must
    /// dereference to MeshSubset*.
    /// \pre All meshes of each of the MeshSubset objects must unique,
    /// an exception is thrown otherwise.
    template <typename InputIterator>
    MeshSubsets(InputIterator const& first, InputIterator const& last)
        : _mesh_subsets(first, last)
    {
        if (!areMeshSubsetMeshesUnique())
        {
            ERR("Mesh ids of input mesh subsets are not unique.");
            std::abort();
        }

        _n_total_items = std::accumulate(first, last, 0u,
            [](std::size_t const& sum, MeshSubset const* const mesh_subset)
            {
                return sum + mesh_subset->getNTotalItems();
            });
        _n_total_nonghost_items = std::accumulate(first, last, 0u,
            [](std::size_t const& sum, MeshSubset const* const mesh_subset)
            {
                return sum + mesh_subset->getNTotalNonGhostItems();
            });
    }

    /// return the total number of mesh items (in all meshes) where this component is assigned
    std::size_t getNMeshItems() const
    {
        return _n_total_items;
    }

    /// return the total number of mesh items (in all meshes) where this component is assigned
    std::size_t getNNonGhostMeshItems() const
    {
        return _n_total_nonghost_items;
    }

    /// Number of saved mesh subsets.
    std::size_t size() const
    {
        return _mesh_subsets.size();
    }

    /// return MeshSubset
    const MeshSubset& getMeshSubset(std::size_t mesh_index) const
    {
        return *_mesh_subsets[mesh_index];
    }

    std::vector<const MeshSubset*>::const_iterator begin() const
    {
        return _mesh_subsets.begin();
    }

    std::vector<const MeshSubset*>::const_iterator end() const
    {
        return _mesh_subsets.end();
    }


    friend std::ostream& operator<<(std::ostream& os, MeshSubsets const& mss)
    {
        os << "Component size: " << mss.size() << "\n";
        for (std::size_t i=0; i< mss._mesh_subsets.size(); i++)
        {
            os << "# Component " << i  << "\n" << *mss._mesh_subsets[i];
        }
        return os;
    }

private:
    /// returns true if all mesh ids of _mesh_subsets elements are different.
    bool areMeshSubsetMeshesUnique() const
    {
        std::vector<std::size_t> mesh_ids;
        std::transform(_mesh_subsets.cbegin(), _mesh_subsets.cend(),
            std::back_inserter(mesh_ids), std::mem_fun(&MeshSubset::getMeshID));

        std::sort(mesh_ids.begin(), mesh_ids.end());
        auto const it = std::adjacent_find(mesh_ids.cbegin(), mesh_ids.cend());

        return mesh_ids.cend() == it;
    }

private:
    std::vector<const MeshSubset*> _mesh_subsets;

    /// Number of all mesh entities on all subsets.
    std::size_t _n_total_items;

    /// Number of all non-ghost mesh entities on all subsets.
    std::size_t _n_total_nonghost_items;
};

}   // namespace MeshLib

#endif    // MESHSUBSETS_H_
