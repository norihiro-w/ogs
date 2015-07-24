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

#include "MeshComponentMap.h"

#include <iostream>
#include <sstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "BaseLib/MPITools.h"
#include "BaseLib/DebugTools.h"
#include "BaseLib/CodingTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"

namespace AssemblerLib
{

using namespace detail;

std::size_t const MeshComponentMap::nop = std::numeric_limits<std::size_t>::max();

MeshComponentMap::MeshComponentMap(
    const std::vector<MeshLib::MeshSubsets*> &components, ComponentOrder order, bool is_parallel)
: _is_parallel(is_parallel)
{
    // construct dict (and here we number global_index by component type)
    std::size_t global_index_offset = 0;
    std::size_t cell_index = 0;
    _n_nonghost_dofs = 0;
    BaseLib::MPIEnvironment mpi;
#ifdef USE_MPI
    std::vector<std::vector<std::size_t>> vec_local_sizes(components.size(), std::vector<std::size_t>(mpi.size())); // comp id -> dom id -> local size
    std::vector<std::size_t> comp_global_sizes(components.size(), 0);
#endif

    if (is_parallel)
    {
#ifdef DBUG_MCP
        for (auto c : components)
        {
            std::stringstream ss;
            ss << *c;
            DBUG("%s", ss.str().data());
        }
#endif

#ifdef USE_MPI
        for (std::size_t i=0; i<components.size(); i++)
        {
            // send local size for each component
            vec_local_sizes[i][mpi.rank()] = components[i]->getNNonGhostMeshItems();
            std::size_t this_size = components[i]->getNNonGhostMeshItems();
            // get local size for each component
            MPI_Allgather(&this_size, 1, MPI_UNSIGNED_LONG, &vec_local_sizes[i][0], 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        }

//    if (mpi_rank==0) {
//        std::cout << "# vec_local_sizes" << std::endl;
//        for (std::size_t i=0; i<components.size(); i++)
//        {
//            std::cout << i << ": ";
//            for (std::size_t n : vec_local_sizes[i])
//                std::cout << n << " ";
//            std::cout << std::endl;
//        }
//    }

        // compute global size for each component
        for (std::size_t i=0; i<components.size(); i++)
        {
            std::size_t sum = 0;
            for (std::size_t n : vec_local_sizes[i])
                sum += n;
            comp_global_sizes[i] = sum;
        }
#endif
    }

    for (std::size_t comp_id=0; comp_id<components.size(); comp_id++)
    {
        auto c = components[comp_id];
#ifdef USE_MPI
        if (is_parallel) {
        if (comp_id>0)
            global_index_offset = comp_global_sizes[comp_id-1];
        for (int i=0; i<mpi.rank(); i++)
            global_index_offset += vec_local_sizes[comp_id][i];
        }
#endif
        for (unsigned mesh_subset_index = 0; mesh_subset_index < c->size(); mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset = c->getMeshSubset(mesh_subset_index);
            std::size_t const mesh_id = mesh_subset.getMeshID();
            auto &mesh = mesh_subset.getMesh();

            for (std::size_t j=0; j<mesh_subset.getNNodes(); j++)
                if (mesh.isGhostNode(mesh.getNode(j)->getID()))
                    _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Node, j), comp_id, true)); //TODO need to set global id
                else {
                    _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Node, j), comp_id, global_index_offset++));
                    _n_nonghost_dofs++;
                }

            // Elements
            for (std::size_t j=0; j<mesh_subset.getNElements(); j++)
                _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Cell, j), comp_id, cell_index++));
       }
    }

    if (order == ComponentOrder::BY_LOCATION)
    {
        std::size_t offset = 0;
#ifdef USE_MPI
        if (is_parallel) {
            for (int i=0; i<mpi.rank(); i++)
                for (std::size_t j=0; j<components.size(); j++)
                    offset += vec_local_sizes[j][i];
        }
#endif
        renumberByLocation(offset);
    }

#ifdef USE_MPI
    if (is_parallel)
    {
        auto &mlc = _dict.get<ByLocationAndComponent>();
        // for each component
        for (std::size_t comp_id = 0; comp_id < components.size(); comp_id++)
        {
            auto c = components[comp_id];
#ifdef DBUG_MCP
            DBUG("# Component %d", comp_id);
#endif
            for (auto itc = c->begin(); itc != c->end(); ++itc)
            {
                auto* ms = *itc;
                auto mesh_id = ms->getMeshID();
                auto &msh = ms->getMesh();

                // get local ghost ids
                std::vector<std::size_t> myghost_local_ids(ms->getGhostNodeIDs());
                std::vector<std::size_t> myghost_global_ids;
                for (auto local_id : myghost_local_ids)
                    myghost_global_ids.push_back(msh.getGlobalNodeID(local_id));
#ifdef DBUG_MCP
                DBUG("Local ghost nodes (local id): %s", BaseLib::toString(myghost_local_ids).data());
                DBUG("Local ghost nodes (global id): %s", BaseLib::toString(myghost_global_ids).data());
#endif

                // spread them and gather all ghost ids
                int n_myghosts = myghost_global_ids.size();
                std::vector<int> vec_ghost_sizes(mpi.size());
                MPI_Allgather(&n_myghosts, 1, MPI_INT, &vec_ghost_sizes[0], 1, MPI_INT, MPI_COMM_WORLD);
#ifdef DBUG_MCP
                DBUG("Ghost sizes: %s", BaseLib::toString(vec_ghost_sizes).data());
#endif
                std::size_t n_totalghosts = std::accumulate(vec_ghost_sizes.begin(), vec_ghost_sizes.end(), 0);
                std::vector<int> vec_offset(vec_ghost_sizes.size());
                std::partial_sum(vec_ghost_sizes.begin(), --vec_ghost_sizes.end(), ++vec_offset.begin());
#ifdef DBUG_MCP
                DBUG("Ghost offset: %s", BaseLib::toString(vec_offset).data());
#endif
                //
                std::vector<std::size_t> allghost_global_ids(n_totalghosts);
                MPI_Allgatherv(&myghost_global_ids[0], myghost_global_ids.size(), MPI_UNSIGNED_LONG,
                               &allghost_global_ids[0], &vec_ghost_sizes[0], &vec_offset[0],
                               MPI_UNSIGNED_LONG, mpi.communicator());
#ifdef DBUG_MCP
                DBUG("Received all ghost nodes: %s", BaseLib::toString(allghost_global_ids).data());
#endif
                std::sort(allghost_global_ids.begin(), allghost_global_ids.end());
                allghost_global_ids.erase(std::unique(allghost_global_ids.begin(), allghost_global_ids.end()), allghost_global_ids.end());
#ifdef DBUG_MCP
                DBUG("Uniqued all ghost nodes: %s", BaseLib::toString(allghost_global_ids).data());
#endif

                // create a list of dof for the ghost nodes whose master exist in this rank
                std::vector<std::size_t> ghost_master_global_id;
                std::vector<std::size_t> ghost_master_dof;
                for (std::size_t i = 0; i < allghost_global_ids.size(); i++)
                {
                    std::size_t global_item_id = allghost_global_ids[i];
                    //check if it contains given global node id and if it is non-ghost
                    auto local_id = msh.getLocalNodeID(global_item_id);
                    if (local_id == BaseLib::index_npos || msh.isGhostNode(local_id) || !ms->hasNode(msh.getNode(local_id)))
                        continue;
                    auto it = mlc.find(Line(Location(mesh_id, MeshLib::MeshItemType::Node, local_id), comp_id));
                    assert(it != mlc.end());
                    auto dof_id = it->global_index;
                    ghost_master_global_id.push_back(allghost_global_ids[i]);
                    ghost_master_dof.push_back(dof_id);
                }
#ifdef DBUG_MCP
                DBUG("Found master of the ghost nodes: %s", BaseLib::toString(ghost_master_global_id).data());
                DBUG("Global Indexes of those        : %s", BaseLib::toString(ghost_master_dof).data());
#endif
                // send the list to root
                if (!mpi.root())
                {
                    std::size_t n_local_ghosts = ghost_master_global_id.size();
#ifdef DBUG_MCP
                    DBUG("Sending data to root...");
#endif
                    MPI_Send(&n_local_ghosts, 1, MPI_UNSIGNED_LONG, 0, 0, mpi.communicator());
                    MPI_Send(&ghost_master_global_id[0], n_local_ghosts,
                    MPI_UNSIGNED_LONG, 0, 0, mpi.communicator());
                    MPI_Send(&ghost_master_dof[0], n_local_ghosts,
                    MPI_UNSIGNED_LONG, 0, 0, mpi.communicator());
                }
                // root collects dof of all ghosts
                std::vector<std::size_t> allghost_dofs(n_totalghosts);
                if (mpi.root())
                {
                    for (std::size_t i = 0; i < ghost_master_global_id.size(); i++)
                    {
                        auto it = std::find(allghost_global_ids.begin(), allghost_global_ids.end(), ghost_master_global_id[i]);
                        if (it == allghost_global_ids.end())
                            continue;
                        std::size_t pos = std::distance(allghost_global_ids.begin(), it);
                        allghost_dofs[pos] = ghost_master_dof[i];
                    }
                    for (int i = 1; i < mpi.size(); i++)
                    {
#ifdef DBUG_MCP
                        DBUG("Receiving data from rank %d...", i);
#endif
                        std::size_t n_local_ghosts;
                        MPI_Status status;
                        MPI_Recv(&n_local_ghosts, 1, MPI_UNSIGNED_LONG, i, 0, mpi.communicator(), &status);
                        std::vector<std::size_t> rank_ghost_ids(n_local_ghosts);
                        MPI_Recv(&rank_ghost_ids[0], n_local_ghosts,
                        MPI_UNSIGNED_LONG, i, 0, mpi.communicator(), &status);
                        std::vector<std::size_t> rank_ghost_dofs(n_local_ghosts);
                        MPI_Recv(&rank_ghost_dofs[0], n_local_ghosts,
                        MPI_UNSIGNED_LONG, i, 0, mpi.communicator(), &status);
                        for (std::size_t i = 0; i < rank_ghost_ids.size(); i++)
                        {
                            auto it = std::find(allghost_global_ids.begin(), allghost_global_ids.end(), rank_ghost_ids[i]);
                            if (it == allghost_global_ids.end())
                                continue;
                            std::size_t pos = std::distance(allghost_global_ids.begin(), it);
                            allghost_dofs[pos] = rank_ghost_dofs[i];
                        }
                    }
                }
                // spread dof of ghosts to other processes
                MPI_Bcast(&allghost_dofs[0], n_totalghosts, MPI_UNSIGNED_LONG, 0, mpi.communicator());
#ifdef DBUG_MCP
                DBUG("Updated ghost dof: %s", BaseLib::toString(allghost_dofs).data());
#endif
                // now comes back to each process.
                // get dof for their local ghost
//            std::stringstream ss;
//            ss << *this;
//            INFO("DoF table\n%s", ss.str().data());
                for (std::size_t i = 0; i < myghost_global_ids.size(); i++)
                {
                    auto global_id = myghost_global_ids[i];
                    auto local_id = myghost_local_ids[i];
                    auto found = std::find(allghost_global_ids.begin(), allghost_global_ids.end(), global_id);
                    if (found == allghost_global_ids.end())
                    {
                        ERR("DoF of ghost node %d not found", global_id);
                        continue;
                    }
                    std::size_t pos = std::distance(allghost_global_ids.begin(), found);
#ifdef DBUG_MCP
                    DBUG("Check my ghost node %d (global %d), pos=%d", local_id, global_id, pos);
#endif
                    auto it = mlc.find(Line(Location(mesh_id, MeshLib::MeshItemType::Node, local_id), comp_id));
                    assert(it != mlc.end());
                    Line l(*it);
                    l.global_index = allghost_dofs[pos];
                    mlc.replace(it, l);
                }
            }

        }

//    std::stringstream ss;
//    ss << *this;
//    INFO("DoF table\n%s", ss.str().data());
    }
#endif
}

MeshComponentMap
MeshComponentMap::getSubset(std::vector<MeshLib::MeshSubsets*> const& components) const
{
    // New dictionary for the subset.
    ComponentGlobalIndexDict subset_dict;

    std::size_t comp_id = 0;
    for (auto c : components)
    {
        if (c == nullptr)   // Empty component
        {
            comp_id++;
            continue;
        }
        for (std::size_t mesh_subset_index = 0; mesh_subset_index < c->size(); mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset = c->getMeshSubset(mesh_subset_index);
            std::size_t const mesh_id = mesh_subset.getMeshID();
            // Lookup the locations in the current mesh component map and
            // insert the full lines into the subset dictionary.
            for (std::size_t j=0; j<mesh_subset.getNNodes(); j++)
                subset_dict.insert(getLine(Location(mesh_id,
                    MeshLib::MeshItemType::Node, mesh_subset.getNodeID(j)), comp_id));
            for (std::size_t j=0; j<mesh_subset.getNElements(); j++)
                subset_dict.insert(getLine(Location(mesh_id,
                    MeshLib::MeshItemType::Cell, mesh_subset.getElementID(j)), comp_id));
        }
        comp_id++;
    }

    return MeshComponentMap(subset_dict);
}

void MeshComponentMap::renumberByLocation(std::size_t offset)
{
    std::size_t global_index = offset;

    auto &m = _dict.get<ByLocation>(); // view as sorted by mesh item
    for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item)
    {
        Line pos = *itr_mesh_item;
        if (!pos.is_ghost)
            pos.global_index = global_index++;
        m.replace(itr_mesh_item, pos);
    }
}

std::vector<std::size_t> MeshComponentMap::getComponentIDs(const Location &l) const
{
    auto const &m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<std::size_t> vec_compID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_compID.push_back(itr->comp_id);
    return vec_compID;
}

Line MeshComponentMap::getLine(Location const& l,
    std::size_t const comp_id) const
{
    auto const &m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    assert(itr != m.end());     // The line must exist in the current dictionary.
    return *itr;
}

std::size_t MeshComponentMap::getGlobalIndex(Location const& l,
    std::size_t const comp_id) const
{
    auto const &m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    return itr!=m.end() ? itr->global_index : nop;
}

std::vector<std::size_t> MeshComponentMap::getGlobalIndices(const Location &l) const
{
    auto const &m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<std::size_t> global_indices;
    for (auto itr=p.first; itr!=p.second; ++itr)
        global_indices.push_back(itr->global_index);
    return global_indices;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_LOCATION>(
    std::vector<Location> const &ls, bool mask_ghosts) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<std::size_t> global_indices;
    global_indices.reserve(ls.size());

    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
        for (auto itr = p.first; itr != p.second; ++itr)
            if (mask_ghosts && itr->is_ghost)
                global_indices.push_back(-1);
            else
                global_indices.push_back(itr->global_index);
    }

    return global_indices;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_COMPONENT>(
    std::vector<Location> const &ls, bool mask_ghosts) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, std::size_t> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(ls.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
        for (auto itr = p.first; itr != p.second; ++itr)
            if (mask_ghosts && itr->is_ghost)
                pairs.emplace_back(itr->comp_id, -1);
            else
                pairs.emplace_back(itr->comp_id, itr->global_index);
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
        {
            return a.first < b.first;
        };

    // Create vector of global indices from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);

    std::vector<std::size_t> global_indices;
    global_indices.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        global_indices.push_back(p->second);

    return global_indices;
}


template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_LOCATION>(
    std::vector<Location> const &ls, std::size_t const comp_id) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<std::size_t> global_indices;
    global_indices.reserve(ls.size());

    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
        for (auto itr = p.first; itr != p.second; ++itr)
            if (itr->comp_id == comp_id)
                global_indices.push_back(itr->global_index);
    }

    return global_indices;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_COMPONENT>(
    std::vector<Location> const &ls, std::size_t const comp_id) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, std::size_t> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(ls.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
        for (auto itr = p.first; itr != p.second; ++itr)
            if (itr->comp_id == comp_id)
                pairs.emplace_back(itr->comp_id, itr->global_index);
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
        {
            return a.first < b.first;
        };

    // Create vector of global indices from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);

    std::vector<std::size_t> global_indices;
    global_indices.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        global_indices.push_back(p->second);

    return global_indices;
}

std::vector<std::size_t> MeshComponentMap::getGlobalIndicesOfGhosts() const
{
    std::vector<std::size_t> vec;
    for (auto it=_dict.begin(); it!=_dict.end(); ++it)
        if (it->is_ghost)
            vec.push_back(it->global_index);
    std::sort(vec.begin(), vec.end());
    return vec;
}

}   // namespace AssemblerLib
