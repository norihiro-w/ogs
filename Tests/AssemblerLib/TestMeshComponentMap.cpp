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

#include <gtest/gtest.h>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "AssemblerLib/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"
#ifdef USE_MPI
#include "GeoLib/Point.h"
#include "GeoLib/AABB.h"
#include "MeshLib/NodePartitionedMesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#endif

class AssemblerLibMeshComponentMapTest : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef AssemblerLib::MeshComponentMap MeshComponentMap;

    public:
    AssemblerLibMeshComponentMapTest()
        : mesh(nullptr), nodesSubset(nullptr), cmap(nullptr)
    {
#ifndef USE_MPI
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
#else
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, 9u);
        // Global: nnodes = 10, nele = 9
        // Local1: nnodes = 6, ghost node=1, nele = 5; ghost ele = 1
        // Local2: nnodes = 6, ghost node=1, nele = 5; ghost ele = 1
        auto global_mesh = mesh;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        std::vector<std::size_t> glb_node_ids;
        std::size_t n_nghost_elem = 4;
        std::size_t n_global_base_nodes = global_mesh->getNBaseNodes();
        std::size_t n_global_nodes = global_mesh->getNNodes();
        std::size_t n_active_base_nodes;
        std::size_t n_active_nodes;
        MeshLib::Mesh* local_msh = nullptr;
        std::vector<MeshLib::Node*> new_nodes;
        if (myrank == 0) {
            std::vector<std::size_t> rm_node_ids;
            for (unsigned i=6; i<10; i++) rm_node_ids.push_back(i);
            local_msh = MeshLib::removeNodes(*global_mesh, rm_node_ids, global_mesh->getName());
            new_nodes = MeshLib::copyNodeVector(local_msh->getNodes());

            std::size_t start = 0, end = 5;
            for (std::size_t i=start; i<end+1; i++)
                glb_node_ids.push_back(i);
            n_active_base_nodes = local_msh->getNBaseNodes() - 1;
            n_active_nodes = local_msh->getNNodes() - 1;
        } else {
            std::vector<std::size_t> rm_node_ids;
            for (unsigned i=0; i<4; i++) rm_node_ids.push_back(i);
            local_msh = MeshLib::removeNodes(*global_mesh, rm_node_ids, global_mesh->getName());
            new_nodes = MeshLib::copyNodeVector(local_msh->getNodes());
            auto nod = new_nodes[0];
            new_nodes.erase(new_nodes.begin());
            new_nodes.push_back(nod);

            std::size_t start = 5, end = 9;
            for (std::size_t i=start; i<end+1; i++)
                glb_node_ids.push_back(i);
            glb_node_ids.push_back(4);
            n_active_base_nodes = local_msh->getNBaseNodes() - 1;
            n_active_nodes = local_msh->getNNodes() - 1;
        }
        std::size_t n_base_nodes = local_msh->getNBaseNodes();
        auto new_elements = MeshLib::copyElementVector(local_msh->getElements(), new_nodes);
        mesh = new MeshLib::NodePartitionedMesh(
                            local_msh->getName(), new_nodes, glb_node_ids,
                            new_elements, local_msh->getProperties(),
                            n_nghost_elem,
                            n_global_base_nodes,
                            n_global_nodes,
                            n_base_nodes,
                            n_active_base_nodes,
                            n_active_nodes);
        delete local_msh;
        delete global_mesh;
#endif
        nodesSubset = new MeshLib::MeshSubset(*mesh, &mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
    }

    ~AssemblerLibMeshComponentMapTest()
    {
        delete cmap;
        std::remove_if(components.begin(), components.end(),
            [](MeshLib::MeshSubsets* p) { delete p; return true; });
        delete nodesSubset;
        delete mesh;
    }

#ifndef USE_MPI
    static std::size_t const mesh_size = 9;
#else
    static std::size_t const mesh_size = 5;
#endif
    MeshLib::Mesh const* mesh;
    MeshLib::MeshSubset const* nodesSubset;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<MeshLib::MeshSubsets*> components;
    MeshComponentMap const* cmap;
#ifdef USE_MPI
    int myrank;
#endif

    //
    // Functions used for checking.
    //

    // Returns global index of a node location and a component.
    std::size_t giAtNodeForComponent(std::size_t const n, std::size_t const c) const
    {
        return cmap->getGlobalIndex(Location(mesh->getID(), MeshItemType::Node, n), c);
    }

};

#ifdef USE_MPI
TEST_F(AssemblerLibMeshComponentMapTest, CheckMeshes)
{
    if (myrank==0)
    {
        ASSERT_EQ(6u, mesh->getNNodes());
        ASSERT_EQ(5u, mesh->getNActiveNodes());
        ASSERT_FALSE(mesh->isGhostNode(4));
        ASSERT_TRUE(mesh->isGhostNode(5));
        ASSERT_EQ(5u, mesh->getNElements());
        ASSERT_EQ(4u, mesh->getNNonGhostElements());
        ASSERT_EQ(0u, mesh->getGlobalNodeID(0));
        ASSERT_EQ(5u, mesh->getGlobalNodeID(mesh->getNNodes()-1));

        ASSERT_EQ(6u, nodesSubset->getNNodes());
        ASSERT_EQ(5u, nodesSubset->getNNonGhostNodes());
        ASSERT_EQ(6u, nodesSubset->getNTotalItems());
        ASSERT_EQ(5u, nodesSubset->getNTotalNonGhostItems());
    } else {
        ASSERT_EQ(6u, mesh->getNNodes());
        ASSERT_EQ(5u, mesh->getNActiveNodes());
//        ASSERT_FALSE(mesh->isGhostNode(4));
//        ASSERT_TRUE(mesh->isGhostNode(5));
        ASSERT_EQ(5u, mesh->getNElements());
        ASSERT_EQ(4u, mesh->getNNonGhostElements());
        ASSERT_EQ(5u, mesh->getGlobalNodeID(0));
        ASSERT_EQ(4u, mesh->getGlobalNodeID(mesh->getNNodes()-1));

        ASSERT_EQ(6u, nodesSubset->getNNodes());
        ASSERT_EQ(5u, nodesSubset->getNNonGhostNodes());
        ASSERT_EQ(6u, nodesSubset->getNTotalItems());
        ASSERT_EQ(5u, nodesSubset->getNTotalNonGhostItems());
    }
    ASSERT_FALSE(mesh->isGlobal());
    ASSERT_EQ(10u, mesh->getNGlobalNodes());
}

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByComponent)
{
    // - Entries in the vector are arranged in the order of a component type and then node ID
    // - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

//    std::stringstream ss;
//    ss << *cmap;
//    std::cout << "rank=\n" << ss.str() << std::endl;

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
    for (std::size_t i = 0; i < mesh->getNActiveNodes(); i++)
    {
        // Test global indices for the different components of the node.
        if (myrank == 0) {
            ASSERT_EQ(i , giAtNodeForComponent(i, comp0_id));
            ASSERT_EQ(mesh->getNGlobalNodes() + i, giAtNodeForComponent(i, comp1_id));
        } else {
            ASSERT_EQ(5 + i , giAtNodeForComponent(i, comp0_id));
            ASSERT_EQ(mesh->getNGlobalNodes() + 5 + i, giAtNodeForComponent(i, comp1_id));
        }

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }
    for (std::size_t i = mesh->getNActiveNodes(); i < mesh->getNNodes(); i++)
    {
        ASSERT_EQ(static_cast<std::size_t>(-1), giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(static_cast<std::size_t>(-1), giAtNodeForComponent(i, comp1_id));
    }
}

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByLocation)
{
    // - Entries in the vector are arranged in the order of node ID and then a component type
    // - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
    for (std::size_t i = 0; i < mesh->getNActiveNodes(); i++)
    {
        // Test global indices for the different components of the node.
        if (myrank == 0) {
            ASSERT_EQ(2 * i , giAtNodeForComponent(i, comp0_id));
            ASSERT_EQ(2 * i + 1, giAtNodeForComponent(i, comp1_id));
        } else {
            ASSERT_EQ(2 * (5 + i) , giAtNodeForComponent(i, comp0_id));
            ASSERT_EQ(2 * (5 + i) + 1, giAtNodeForComponent(i, comp1_id));
        }

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }

    for (std::size_t i = mesh->getNActiveNodes(); i < mesh->getNNodes(); i++)
    {
        ASSERT_EQ(static_cast<std::size_t>(-1), giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(static_cast<std::size_t>(-1), giAtNodeForComponent(i, comp1_id));
    }
}

#else

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByComponent)
{
    // - Entries in the vector are arranged in the order of a component type and then node ID
    // - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
    for (std::size_t i = 0; i < mesh_size; i++)
    {
        // Test global indices for the different components of the node.
        ASSERT_EQ(i , giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(mesh_size + 1 + i, giAtNodeForComponent(i, comp1_id));

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }
}

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByLocation)
{
    // - Entries in the vector are arranged in the order of node ID and then a component type
    // - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
    for (std::size_t i = 0; i < mesh_size; i++)
    {
        // Test global indices for the different components of the node.
        ASSERT_EQ(2 * i , giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(2 * i + 1, giAtNodeForComponent(i, comp1_id));

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }
}

TEST_F(AssemblerLibMeshComponentMapTest, OutOfRangeAccess)
{
    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, mesh_size + 1), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID() + 1, MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Cell, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, 0), 10));
}

TEST_F(AssemblerLibMeshComponentMapTest, SubsetOfNodesByComponent)
{
    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    // Select some nodes from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 9 }};
    std::vector<MeshLib::Node*> some_nodes;
    for (std::size_t id : ids)
        some_nodes.push_back(const_cast<MeshLib::Node*>(mesh->getNode(id)));

    MeshLib::MeshSubset some_nodes_mesh_subset(*mesh, &some_nodes);

    std::vector<MeshLib::MeshSubsets*> selected_components;
    selected_components.emplace_back(nullptr);  // empty component
    selected_components.emplace_back(new MeshLib::MeshSubsets(&some_nodes_mesh_subset));

    // Subset the original cmap.
    MeshComponentMap cmap_subset = cmap->getSubset(selected_components);

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.size());

    // .. and the content of the subset.
    for (std::size_t id : ids)
    {
        Location const l(mesh->getID(), MeshItemType::Node, id);
        EXPECT_EQ(cmap->getGlobalIndex(l, comp1_id),
            cmap_subset.getGlobalIndex(l, comp1_id));
    }

    for (auto p : selected_components)
        delete p;
}

TEST_F(AssemblerLibMeshComponentMapTest, SubsetOfNodesByLocation)
{
    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    // Select some nodes from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 9 }};
    std::vector<MeshLib::Node*> some_nodes;
    for (std::size_t id : ids)
        some_nodes.push_back(const_cast<MeshLib::Node*>(mesh->getNode(id)));

    MeshLib::MeshSubset some_nodes_mesh_subset(*mesh, &some_nodes);

    std::vector<MeshLib::MeshSubsets*> selected_components;
    selected_components.emplace_back(nullptr);  // empty component
    selected_components.emplace_back(new MeshLib::MeshSubsets(&some_nodes_mesh_subset));

    // Subset the original cmap.
    MeshComponentMap cmap_subset = cmap->getSubset(selected_components);

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.size());

    // .. and the content of the subset.
    for (std::size_t id : ids)
    {
        Location const l(mesh->getID(), MeshItemType::Node, id);
        EXPECT_EQ(cmap->getGlobalIndex(l, comp1_id),
            cmap_subset.getGlobalIndex(l, comp1_id));
    }

    for (auto p : selected_components)
        delete p;
}

#endif
