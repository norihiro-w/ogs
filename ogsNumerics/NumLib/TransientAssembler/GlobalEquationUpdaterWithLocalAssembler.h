/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GlobalEquationUpdatorWithLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "DiscreteLib/Utils/Tools.h"

namespace NumLib
{

template <class T_LOCAL>
class TransientGlobalEquationUpdaterWithLocalAssembler
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

    TransientGlobalEquationUpdaterWithLocalAssembler(const TimeStep* time, MeshLib::IMesh* msh, DiscreteLib::DofEquationIdTable* dofManager, const GlobalVector* u0, const GlobalVector* u1, LocalAssemblerType* a)
    : _msh(msh), _dofManager(dofManager), _transient_e_assembler(a), _timestep(time), _vec_u0(u0), _vec_u1(u1)
    {

    }

    void update(const MeshLib::IElement &e, MathLib::ILinearEquation &eqs)
    {
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap_row;
        //std::vector<size_t> local_dofmap_column;
        LocalEquation localEQS;
        LocalVector local_u_n1;
        LocalVector local_u_n;

        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        _dofManager->mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_row, DiscreteLib::DofNumberingType::BY_POINT); //TODO order
        // previous and current results
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u1, local_u_n1);
        DiscreteLib::getLocalVector(local_dofmap_row, *_vec_u0, local_u_n);
        // local assembly
        localEQS.create(local_dofmap_row.size());
        _transient_e_assembler->assembly(*_timestep, e, local_u_n1, local_u_n, localEQS);

//        if (i<3) {
//            std::cout << "local A = \n" << *localEQS.getA();
//            std::cout << "local b = \n" << *localEQS.getRHS();
//        }

        // update global
        eqs.addAsub(local_dofmap_row, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }

private:
    MeshLib::IMesh* _msh;
    DiscreteLib::DofEquationIdTable* _dofManager;
    LocalAssemblerType* _transient_e_assembler;
    const TimeStep* _timestep;
    const GlobalVector* _vec_u0;
    const GlobalVector* _vec_u1;
};

} //end
