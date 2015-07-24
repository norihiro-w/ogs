/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_MAPTOOLS_H_
#define ASSEMBLERLIB_MAPTOOLS_H_

#include "MathLib/LinAlg/IVector.h"
#include "MathLib/DataType.h"

#include "MeshLib/Location.h"

#include "LocalToGlobalIndexMap.h"

namespace AssemblerLib
{

inline void setLocalVector(
		const LocalToGlobalIndexMap &dofManager,
		size_t comp_id, size_t mesh_id,
		const MathLib::IVector &global_vec, MathLib::IVector &node_value_vec)
{
	auto &meshComponentMap(dofManager.getMeshComponentMap());
	for (size_t i=node_value_vec.getRangeBegin(); i<node_value_vec.getRangeEnd(); i++) {
		MeshLib::Location loc(mesh_id, MeshLib::MeshItemType::Node, i);
		size_t eqs_id = meshComponentMap.getGlobalIndex(loc, comp_id);
		node_value_vec.set(i, global_vec[eqs_id]);
	}
}

inline void getLocalVector(const std::vector<std::size_t> &list_vec_entry_id, const MathLib::LocalVector &global_u, MathLib::LocalVector  &local_u)
{
    size_t valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            valid_entry_cnt++;
    }
    local_u.resize(valid_entry_cnt);
    valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            local_u[valid_entry_cnt++] = global_u[list_vec_entry_id[i]];
    }
}


inline MathLib::LocalVector getLocalVector(const std::vector<std::size_t> &list_vec_entry_id, const MathLib::IVector &global_u)
{
    size_t valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            valid_entry_cnt++;
    }
    MathLib::LocalVector local_u(valid_entry_cnt);
    valid_entry_cnt = 0;
    for (size_t i=0; i<list_vec_entry_id.size(); i++) {
        if (list_vec_entry_id[i] != BaseLib::index_npos)
            local_u[valid_entry_cnt++] = global_u[list_vec_entry_id[i]];
    }
    return local_u;
}

/**
 * copy node values of a component to a solution vector which may consist of multiple components
 *
 * @param dofManager      mapping table
 * @param comp_id         component ID
 * @param mesh_id         mesh ID
 * @param node_value_vec  node value vector
 * @param global_vec      solution vector
 */
inline void setGlobalVector(
		const LocalToGlobalIndexMap &dofManager,
		size_t comp_id, size_t mesh_id,
		const MathLib::IVector &node_value_vec, MathLib::IVector &global_vec)
{
	auto &meshComponentMap(dofManager.getMeshComponentMap());
	for (size_t i=node_value_vec.getRangeBegin(); i<node_value_vec.getRangeEnd(); i++) {
		MeshLib::Location loc(mesh_id, MeshLib::MeshItemType::Node, i);
		size_t eqs_id = meshComponentMap.getGlobalIndex(loc, comp_id);
		if (eqs_id != static_cast<std::size_t>(-1))
			global_vec.set(eqs_id, node_value_vec[i]);
	}
}


/**
 * convert vectors of node id and value to vectors of global id and value
 * @param dofManager        mapping table
 * @param comp_id           component ID
 * @param msh_id            mesh ID
 * @param list_node_id      a vector of node IDs
 * @param list_node_values  a vector of node values
 * @param list_eqs_id       a vector of global IDs
 * @param list_eqs_val      a vector of values
 */
inline void toGlobalIDValues(
		const LocalToGlobalIndexMap &dofManager,
		size_t comp_id, size_t msh_id,
		const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values,
		std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val)
{
	auto &meshComponentMap(dofManager.getMeshComponentMap());
	for (size_t j=0; j<list_node_id.size(); j++) {
		MeshLib::Location loc(msh_id, MeshLib::MeshItemType::Node, list_node_id[j]);
		size_t eqs_id = meshComponentMap.getGlobalIndex(loc, comp_id);
		if (MeshComponentMap::nop == eqs_id)
			continue;
		double bc_val = list_node_values[j];
		list_eqs_id.push_back(eqs_id);
		list_eqs_val.push_back(bc_val);
	}
}

}  // namespace AssemblerLib

#endif  // ASSEMBLERLIB_MAPTOOLS_H_
