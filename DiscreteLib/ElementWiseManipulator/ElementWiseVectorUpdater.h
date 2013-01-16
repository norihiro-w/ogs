/**
 * \file   ElementWiseVectorUpdater.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTWISEVECTORUPDATER_H_
#define ELEMENTWISEVECTORUPDATER_H_

#include <vector>
#include "MeshLib/Elements/Element.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise vector updater
 *
 * \tparam T_VALUE  vector data type
 * \tparam T_LOCAL  local assembler
 */
template <class T_VALUE, class T_LOCAL>
class ElementWiseVectorUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;

    /**
     * constructor
     * @param msh
     * @param a
     */
    ElementWiseVectorUpdater(std::size_t mesh_id, LocalAssemblerType* a)
    : _msh_id(mesh_id), _e_assembler(a)
    {
    }

    /**
     * update a discrete vector for the given element
     * @param e             Mesh element
     * @param dofManager    DoF mapping table
     * @param globalVec     Discrete vector to be updated
     */
    void update(const MeshLib::Element &e, const DofEquationIdTable &dofManager, GlobalVectorType &globalVec)
    {
        LocalVector localVec;
        std::vector<std::size_t> ele_node_ids, ele_node_size_order;
        std::vector<long> local_dofmap_row;//, local_dofmap_column;

        std::vector<T_VALUE> local_u_n;
        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(_msh_id, ele_node_ids, local_dofmap_row); //TODO order
        // local assembly
        localVec.resize(local_dofmap_row.size(), .0);
        _e_assembler->assembly(*e, localVec);
        // update global
        globalVec.addSubvector(local_dofmap_row, &localVec[0]);
    }

private:
    const std::size_t _msh_id;
    LocalAssemblerType* _e_assembler;
};

} //end

#endif //ELEMENTWISEVECTORUPDATER_H_
