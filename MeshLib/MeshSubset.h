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

#ifndef MESHSUBSET_H_
#define MESHSUBSET_H_

#include <cassert>
#include <vector>

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"


namespace MeshLib
{

/// A subset of nodes or elements on a single mesh.
class MeshSubset
{
public:
    /// Construct a mesh subset from vector of nodes on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Node pointers.
    /// \param delete_ptr Deletes the vector of Node pointers if true.
    /// \note When delete_ptr is set only the vector is deleted, not the
    /// elements of the vector.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_items,
        bool const delete_ptr = false)
        : _msh(msh), _nodes(vec_items), _eles(nullptr), _delete_ptr(delete_ptr),
          _n_non_ghost_nodes(computeNNonGhostNodes(msh, *vec_items)), _n_non_ghost_elements(0)
    {}

    /// Construct a mesh subset from vector of elements on the given mesh.
    /// \param msh Mesh
    /// \param vec_items Vector of Element pointers.
    /// \param delete_ptr Deletes the vector of Element pointers if true.
    /// \note When delete_ptr is set only the vector is deleted, not the
    /// elements of the vector.
    MeshSubset(const Mesh& msh, std::vector<Element*> const* vec_items,
        bool const delete_ptr = false)
        : _msh(msh), _nodes(nullptr), _eles(vec_items), _delete_ptr(delete_ptr),
          _n_non_ghost_nodes(0), _n_non_ghost_elements(computeNNonGhostElements(msh, *vec_items))
    {}

    /// construct from both nodes and elements
    /// Construct a mesh subset from vector of nodes and a vector of elements on
    /// the given mesh.
    /// \param msh Mesh
    /// \param vec_nodes Vector of Node pointers.
    /// \param vec_eles Vector of Element pointers.
    /// \param delete_ptr Deletes the vector of Node pointers if true.
    /// \note When delete_ptr is set only the vectors are deleted, not the
    /// elements of the vectors.
    MeshSubset(const Mesh& msh, std::vector<Node*> const* vec_nodes,
        std::vector<Element*> const* vec_eles, bool const delete_ptr = false)
        : _msh(msh), _nodes(vec_nodes), _eles(vec_eles), _delete_ptr(delete_ptr),
          _n_non_ghost_nodes(computeNNonGhostNodes(msh, *vec_nodes)),
          _n_non_ghost_elements(computeNNonGhostElements(msh, *vec_eles))
    {}

    ~MeshSubset()
    {
        if (_delete_ptr)
        {
            delete _nodes;
            delete _eles;
        }
    }

    /// return the total number of mesh items
    std::size_t getNTotalItems() const
    {
        return getNNodes() + getNElements();
    }

    /// return the total number of mesh items
    std::size_t getNTotalNonGhostItems() const
    {
        return getNNonGhostNodes() + getNNonGhostElements();
    }

    /// return this mesh ID
    std::size_t getMeshID() const
    {
        return _msh.getID();
    }

    /// return the number of registered nodes
    std::size_t getNNodes() const
    {
        return (_nodes==nullptr) ? 0 : _nodes->size();
    }

    /// return the number of registered non-ghost nodes
    std::size_t getNNonGhostNodes() const
    {
        return _n_non_ghost_nodes;
    }

    /// return the number of registered ghost nodes
    std::size_t getNGhostNodes() const
    {
        return (getNNodes() - getNNonGhostNodes());
    }

    /// Returns the global node id Node::getID() of i-th node in the mesh
    /// subset.
    /// \pre The _nodes must be a valid pointer to a vector of size > i.
    std::size_t getNodeID(std::size_t const i) const
    {
        assert(_nodes && i < _nodes->size());
        return (*_nodes)[i]->getID();
    }

    /// return the number of registered elements
    std::size_t getNElements() const
    {
        return (_eles==nullptr) ? 0 : _eles->size();
    }

    /// return the number of registered non-ghost elements
    std::size_t getNNonGhostElements() const
    {
        return _n_non_ghost_elements;
    }


    /// Returns the global element id Element::getID() of i-th element in the
    /// mesh subset.
    /// \pre The _eles must be a valid pointer to a vector of size > i.
    std::size_t getElementID(std::size_t const i) const
    {
        assert(_eles && i < _eles->size());
        return (*_eles)[i]->getID();
    }

    std::vector<Element*>::const_iterator elementsBegin() const
    {
        return _msh.getElements().cbegin();
    }

    std::vector<Element*>::const_iterator elementsEnd() const
    {
        return _msh.getElements().cend();
    }

    /// Constructs a new mesh subset which is a set intersection of the current
    /// nodes and the provided vector of nodes.
    /// An empty mesh subset may be returned, not a nullptr, in case of empty
    /// intersection or empty input vector.
    MeshSubset*
    getIntersectionByNodes(std::vector<Node*> const& nodes) const
    {
        std::vector<Node*>* active_nodes = new std::vector<Node*>;

        if (_nodes == nullptr || _nodes->empty())
            return new MeshSubset(_msh, active_nodes);   // Empty mesh subset

        for (auto n : nodes)
        {
            auto it = std::find(_nodes->cbegin(), _nodes->cend(), n);
            if (it == _nodes->cend())
                continue;
            active_nodes->push_back(n);
        }

        // Transfer the ownership of active_nodes to the new MeshSubset, which
        // deletes the pointer itself.
        return new MeshSubset(_msh, active_nodes, true);
    }

    const Mesh &getMesh() const
    {
        return _msh;
    }

    std::vector<std::size_t> getGhostNodeIDs() const
    {
        std::vector<std::size_t> ghost_ids;
        for (auto node : *_nodes)
            if (_msh.isGhostNode(node->getID()))
                ghost_ids.push_back(node->getID());
        return ghost_ids;
    }

    bool hasNode(Node const* node) const
    {
        return (std::find(_nodes->begin(), _nodes->end(), node)!=_nodes->end());
    }

    friend std::ostream& operator<<(std::ostream& os, MeshSubset const& ms)
    {
        os << "Mesh ID = " << ms._msh.getID() << ", NNodes = " << (ms._nodes ? ms._nodes->size() : 0) << ", NElements = " << (ms._eles ? ms._eles->size() : 0) << "\n";
        if (ms._nodes) {
            os << "Node local IDs: ";
            for (std::size_t i=0; i< ms._nodes->size(); i++)
                os << (*ms._nodes)[i]->getID() << " ";
            os << "\n";
            os << "Node global IDs: ";
            for (std::size_t i=0; i< ms._nodes->size(); i++)
                os << ms.getMesh().getGlobalNodeID((*ms._nodes)[i]->getID()) << " ";
            os << "\n";
        }
        if (ms._eles) {
            os << "Elements: ";
            for (std::size_t i=0; i< ms._eles->size(); i++)
                os << (*ms._eles)[i]->getID() << " ";
            os << "\n";
        }
        return os;
    }

private:
    Mesh const& _msh;
    std::vector<Node*> const* _nodes;
    std::vector<Element*> const* _eles;
    bool const _delete_ptr = false;
    std::size_t _n_non_ghost_nodes;
    std::size_t _n_non_ghost_elements;

    static std::size_t computeNNonGhostNodes(Mesh const& msh, std::vector<Node*> const& nodes)
    {
        std::size_t n = 0;
        for (auto node : nodes)
            if (!msh.isGhostNode(node->getID()))
                n++;
        return n;
    }

    static std::size_t computeNNonGhostElements(Mesh const& msh, std::vector<Element*> const& elements)
    {
        std::size_t n = 0;
        for (auto e : elements)
            if (!msh.isGhostElement(e->getID()))
                n++;
        return n;
    }
};

}   // namespace MeshLib

#endif  // MESHSUBSET_H_
