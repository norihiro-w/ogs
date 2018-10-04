/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include "BaseLib/Algorithm.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
namespace
{
// A class to check whether a node is located on a crack tip with
// the following conditions:
// - the number of connected fracture elements is one
// - the node is not located on a domain boundary
class IsCrackTip
{
public:
    explicit IsCrackTip(MeshLib::Mesh const& mesh)
        : _mesh(mesh), _fracture_element_dim(mesh.getDimension() - 1)
    {
        _is_internal_node.resize(mesh.getNumberOfNodes(), true);

        MeshLib::NodeSearch nodeSearch(mesh);
        nodeSearch.searchBoundaryNodes();
        for (auto i : nodeSearch.getSearchedNodeIDs())
            _is_internal_node[i] = false;
    }

    bool operator()(MeshLib::Node const& node) const
    {
        if (!_is_internal_node[node.getID()] || !_mesh.isBaseNode(node.getID()))
            return false;

        unsigned n_connected_fracture_elements = 0;
        for (MeshLib::Element const* e : node.getElements())
            if (e->getDimension() == _fracture_element_dim)
                n_connected_fracture_elements++;
        assert(n_connected_fracture_elements > 0);

        return (n_connected_fracture_elements == 1);
    }

private:
    MeshLib::Mesh const& _mesh;
    unsigned const _fracture_element_dim;
    std::vector<bool> _is_internal_node;
};


void findFracutreIntersections(
    MeshLib::Mesh const& mesh,
    std::vector<int> const& vec_fracture_mat_IDs,
    std::vector<std::vector<MeshLib::Node*>> const& vec_fracture_nodes,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_branch_nodeID_matIDs,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_junction_nodeID_matIDs
    )
{
    // make a vector all fracture nodes
    std::vector<std::size_t> all_fracture_nodes;
    for (unsigned mat_id = 0; mat_id < vec_fracture_mat_IDs.size(); mat_id++)
        for (auto* node : vec_fracture_nodes[mat_id])
            all_fracture_nodes.push_back(node->getID());

    // create a table of a node id and connected material IDs
    std::map<std::size_t, std::vector<std::size_t>> frac_nodeID_to_matIDs;
    for (unsigned mat_id = 0; mat_id < vec_fracture_mat_IDs.size(); mat_id++)
        for (auto* node : vec_fracture_nodes[mat_id])
            frac_nodeID_to_matIDs[node->getID()].push_back(mat_id);

    auto opt_material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));

    // find branch/junction nodes which connect to multiple fractures
    for (auto entry : frac_nodeID_to_matIDs)
    {
        auto nodeID = entry.first;
        auto const* node = mesh.getNode(entry.first);
        auto const& matIDs = entry.second;
        if (matIDs.size() == 1)
            continue;

        std::map<int,int> vec_matID_counts;
        {
            for (auto matid : matIDs)
                vec_matID_counts[matid] = 0;

            auto const& connected_elements = node->getElements();
            for (auto const* e : connected_elements)
            {
                if (e->getDimension() == mesh.getDimension())
                    continue;
                auto matid = (*opt_material_ids)[e->getID()];
                vec_matID_counts[matid]++;
            }
        }

        bool isBranch = false;
        {
            for (auto entry : vec_matID_counts)
            {
                auto count = entry.second;
                if (count%2==1) {
                    isBranch = true;
                    break;
                }
            }
        }

        if (isBranch)
        {
            std::vector<int> matIDs(2);
            for (auto entry : vec_matID_counts)
            {
                auto matid = entry.first;
                auto count = entry.second;
                if (count%2==0) {
                    matIDs[0] = matid; // master
                } else {
                    matIDs[1] = matid; // slave
                }
            }
            vec_branch_nodeID_matIDs.push_back(std::make_pair(nodeID,matIDs));

        } else {
            std::vector<int> matIDs(2);
            matIDs[0] = std::min(vec_matID_counts.begin()->first, vec_matID_counts.rbegin()->first);
            matIDs[0] = std::max(vec_matID_counts.begin()->first, vec_matID_counts.rbegin()->first);
            vec_junction_nodeID_matIDs.push_back(std::make_pair(nodeID,matIDs));
        }
    }
}

}  // namespace

void getFractureMatrixDataInMesh(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*>& vec_matrix_elements,
    std::vector<int>& vec_fracture_mat_IDs,
    std::vector<std::vector<MeshLib::Element*>>& vec_fracture_elements,
    std::vector<std::vector<MeshLib::Element*>>& vec_fracture_matrix_elements,
    std::vector<std::vector<MeshLib::Node*>>& vec_fracture_nodes,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_branch_nodeID_matIDs,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_junction_nodeID_matIDs)
{
    IsCrackTip isCrackTip(mesh);

    // get vectors of matrix elements and fracture elements
    vec_matrix_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_fracture_elements;
    for (MeshLib::Element* e : mesh.getElements())
    {
        if (e->getDimension() == mesh.getDimension())
            vec_matrix_elements.push_back(e);
        else
            all_fracture_elements.push_back(e);
    }
    DBUG("-> found total %d matrix elements and %d fracture elements",
         vec_matrix_elements.size(), all_fracture_elements.size());

    // get fracture material IDs
    auto opt_material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    for (MeshLib::Element* e : all_fracture_elements)
        vec_fracture_mat_IDs.push_back((*opt_material_ids)[e->getID()]);
    BaseLib::makeVectorUnique(vec_fracture_mat_IDs);
    DBUG("-> found %d fracture material groups", vec_fracture_mat_IDs.size());

    // create a vector of fracture elements for each material
    vec_fracture_elements.resize(vec_fracture_mat_IDs.size());
    for (unsigned frac_id = 0; frac_id < vec_fracture_mat_IDs.size(); frac_id++)
    {
        const auto frac_mat_id = vec_fracture_mat_IDs[frac_id];
        std::vector<MeshLib::Element*>& vec_elements =
            vec_fracture_elements[frac_id];
        std::copy_if(all_fracture_elements.begin(), all_fracture_elements.end(),
                     std::back_inserter(vec_elements),
                     [&](MeshLib::Element* e) {
                         return (*opt_material_ids)[e->getID()] == frac_mat_id;
                     });
        DBUG("-> found %d elements on the fracture %d", vec_elements.size(),
             frac_id);
    }

    // get a vector of fracture nodes for each material
    vec_fracture_nodes.resize(vec_fracture_mat_IDs.size());
    for (unsigned frac_id = 0; frac_id < vec_fracture_mat_IDs.size(); frac_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = vec_fracture_nodes[frac_id];
        for (MeshLib::Element* e : vec_fracture_elements[frac_id])
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                if (isCrackTip(*e->getNode(i)))
                    continue;
                vec_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
            }
        }
        BaseLib::makeVectorUnique(
            vec_nodes, [](MeshLib::Node* node1, MeshLib::Node* node2) {
                return node1->getID() < node2->getID();
            });
        DBUG("-> found %d nodes on the fracture %d", vec_nodes.size(), frac_id);
    }

    // find branch/junction nodes which connect to multiple fractures
    findFracutreIntersections(
        mesh, vec_fracture_mat_IDs, vec_fracture_nodes,
        vec_branch_nodeID_matIDs, vec_junction_nodeID_matIDs);

    // create a vector fracture elements and connected matrix elements,
    // which are passed to a DoF table
    for (auto fracture_elements : vec_fracture_elements)
    {
        std::vector<MeshLib::Element*> vec_ele;
        // first, collect matrix elements
        for (MeshLib::Element* e : fracture_elements)
        {
            // it is sufficient to iterate over base nodes, because they are
            // already connected to all neighbours
            for (unsigned i = 0; i < e->getNumberOfBaseNodes(); i++)
            {
                MeshLib::Node const* node = e->getNode(i);
                if (isCrackTip(*node))
                    continue;
                for (unsigned j = 0; j < node->getNumberOfElements(); j++)
                {
                    // only matrix elements
                    if (node->getElement(j)->getDimension() <
                        mesh.getDimension())
                        continue;
                    vec_ele.push_back(
                        const_cast<MeshLib::Element*>(node->getElement(j)));
                }
            }
        }
        BaseLib::makeVectorUnique(
            vec_ele, [](MeshLib::Element* e1, MeshLib::Element* e2) {
                return e1->getID() < e2->getID();
            });

        // second, append fracture elements
        std::copy(fracture_elements.begin(), fracture_elements.end(),
                  std::back_inserter(vec_ele));

        vec_fracture_matrix_elements.push_back(vec_ele);
    }
}

}  // namespace LIE
}  // namespace ProcessLib
