/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <string>
#include <vector>

#include "BaseLib/Counter.h"

#include "Node.h"
#include "MeshEnums.h"
#include "Properties.h"
#include "CoordinateSystem.h"

namespace MeshLib
{
	class Element;

/**
 * A basic mesh.
 */
class Mesh : BaseLib::Counter<Mesh>
{
	/* friend functions: */
	friend void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &nodes);

public:
	/// Constructor using a mesh name and an array of nodes and elements
	/// @param name          Mesh name.
	/// @param nodes         A vector of mesh nodes. In case nonlinear nodes are involved, one should
	///                      put them after line ones in the vector and set "n_base_nodes" argument.
	/// @param elements      An array of mesh elements.
	/// @param n_base_nodes  The number of base nodes. This is an optional parameter for nonlinear case.
	///                      If the parameter is set to zero, we consider there are no nonlinear nodes.
	Mesh(const std::string &name,
	     const std::vector<Node*> &nodes,
	     const std::vector<Element*> &elements,
	     Properties const& properties = Properties(),
	     const std::size_t n_base_nodes = 0);

	/// Copy constructor
	Mesh(const Mesh &mesh, std::size_t mesh_id = -1);

	/// Destructor
	virtual ~Mesh();

	/// Add a node to the mesh.
	void addNode(Node* node);

	/// Add an element to the mesh.
	void addElement(Element* elem);

	/// Returns the dimension of the mesh (determined by the maximum dimension over all elements).
	unsigned getDimension() const { return _mesh_dimension; }

	/// Returns the coordinate system of the mesh
	CoordinateSystem getCoordinateSystem() const { return _coord_system; }

	/// Get the node with the given index.
	const Node* getNode(unsigned idx) const { return _nodes[idx]; }

	/// Get the element with the given index.
	const Element* getElement(unsigned idx) const { return _elements[idx]; }

	/// Get the minimum edge length over all elements of the mesh.
	double getMinEdgeLength() const { return _edge_length.first; }

	/// Get the maximum edge length over all elements of the mesh.
	double getMaxEdgeLength() const { return _edge_length.second; }

	/// Get the minimum node distance in the mesh.
    /// The value is calculated from element-wise minimum node distances.
	double getMinNodeDistance() const { return _node_distance.first; }

	/// Get the maximum node distance over all elements of the mesh.
    /// The value is calculated from element-wise maximum node distances.
	double getMaxNodeDistance() const { return _node_distance.second; }

	/// Get the number of elements
	std::size_t getNElements() const { return _elements.size(); }

	/// Get the number of nodes
	std::size_t getNNodes() const { return _nodes.size(); }

	/// Get name of the mesh.
	const std::string getName() const { return _name; }

	/// Get the nodes-vector for the mesh.
	std::vector<Node*> const& getNodes() const { return _nodes; }

	/// Get the element-vector for the mesh.
	std::vector<Element*> const& getElements() const { return _elements; }

	/// Resets the IDs of all mesh-elements to their position in the element vector
	void resetElementIDs();
	
	/// Resets the IDs of all mesh-nodes to their position in the node vector
	void resetNodeIDs();

	/// Changes the name of the mesh.
	void setName(const std::string &name) { this->_name = name; }

	/// Get id of the mesh
	std::size_t getID() const {return _id; }

	/// Get the number of base nodes
	std::size_t getNBaseNodes() const { return _n_base_nodes; }

	/// Return true if the given node is a basic one (i.e. linear order node)
	bool isBaseNode(std::size_t node_idx) const {return node_idx < _n_base_nodes; }

	/// Return true if the mesh has any nonlinear nodes
	bool isNonlinear() const { return (getNNodes() != getNBaseNodes()); }

	MeshLib::Properties & getProperties() { return _properties; }
	MeshLib::Properties const& getProperties() const { return _properties; }

    /// Get the number of nodes of the global mesh for linear elements.
    virtual std::size_t getNGlobalBaseNodes() const
    {
        return getNBaseNodes();
    }

    /// Get the number of all nodes of the global mesh.
    virtual std::size_t getNGlobalNodes() const
    {
        return getNNodes();
    }

    /// Get the global node ID of a node with its local ID.
    virtual std::size_t getGlobalNodeID(const std::size_t node_id) const
    {
        return node_id;
    }

    /// Get the number of the active nodes of the partition for linear elements.
    virtual std::size_t getNActiveBaseNodes() const
    {
        return getNBaseNodes();
    }

    /// Get the number of all active nodes of the partition.
    virtual std::size_t getNActiveNodes() const
    {
        return getNNodes();
    }

    virtual std::size_t getNGhostNodes() const { return getNNodes() - getNActiveNodes(); }

    virtual std::size_t getNGhostBaseNodes() const { return getNBaseNodes() - getNActiveBaseNodes(); }

    virtual const Node* getGhostNode(std::size_t i) const
    {
        return getNode(getNActiveNodes() + i);
    }

    /// Check whether a node with given id is a ghost node.
    virtual bool isGhostNode(const std::size_t /*node_id*/) const
    {
        return false;
    }

    /// Check whether an element with given id is ghost.
    virtual bool isGhostElement(const std::size_t /*element_id*/) const
    {
        return false;
    }

    /// Get the largest ID of active nodes for higher order elements in a partition.
    virtual std::size_t getLargestActiveNodeID() const
    {
        const Node* node = _nodes.back();
        return node->getID();
    }

    /// Get the number of non-ghost elements, or the start entry ID of ghost elements in element vector.
    virtual std::size_t getNNonGhostElements() const
    {
        return getNElements();
    }

    bool isGlobal() const { return getNNodes() == getNGlobalNodes(); }

    virtual bool isPartitioned() const { return false; }

    virtual bool hasGlobalNode(std::size_t global_node_id) const { return global_node_id < getNNodes(); }
    virtual std::size_t getLocalNodeID(std::size_t global_node_id) const { return hasGlobalNode(global_node_id) ? global_node_id : -1; }

protected:
	/// Set the minimum and maximum length over the edges of the mesh.
	void calcEdgeLengthRange();
	/// Set the minimum and maximum node distances within elements.
	void calcNodeDistanceRange();

	/**
	 * Resets the connected elements for the node vector, i.e. removes the old information and
	 * calls setElementsConnectedToNodes to set the new information.
	 * \attention This needs to be called if node neighbourhoods are reset.
	 */
	void resetElementsConnectedToNodes();

	/// Sets the dimension of the mesh.
	void setDimension();

	/// Fills in the neighbor-information for nodes (i.e. which element each node belongs to).
	void setElementsConnectedToNodes();

	/// Fills in the neighbor-information for elements.
	/// Note: Using this implementation, an element e can only have neighbors that have the same dimensionality as e.
	void setElementNeighbors();

	void setNodesConnectedByEdges();

	/// Computes the element-connectivity of nodes. Two nodes i and j are
	/// connected if they are shared by an element.
	void setNodesConnectedByElements();

	std::size_t const _id;
	unsigned _mesh_dimension;
	/// The minimal and maximal edge length over all elements in the mesh
	std::pair<double, double> _edge_length;
	/// The minimal and maximal distance of nodes within an element over all elements in the mesh
	std::pair<double, double> _node_distance;
	std::string _name;
	std::vector<Node*> _nodes;
	std::vector<Element*> _elements;
	std::size_t _n_base_nodes;
	Properties _properties;
	CoordinateSystem _coord_system;
}; /* class */

} /* namespace */

#endif /* MESH_H_ */

