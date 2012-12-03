/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepKinReduction.h
 *
 * Created on 2012-09-27 by Haibing Shao
 */
 
#include "ReductionKinNodeInfo.h"
 
ReductionKinNodeInfo::ReductionKinNodeInfo(size_t node_id, 
	                                       size_t n_comp, 
	                                       size_t n_eta_mob, 
										   size_t n_eta_immob, 
										   size_t n_xi_mob,
										   size_t n_xi_immob,
										   ogsChem::chemReductionKin* ReductionKin)
    : _node_id(node_id), _n_comp(n_comp), _n_eta_mob(n_eta_mob), 
	  _n_eta_immob(n_eta_immob), _n_xi_mob(n_xi_mob), _n_xi_immob(n_xi_immob), _ReductionKin(ReductionKin)
{
	_Comp_Conc = MathLib::LocalVector::Zero( _n_comp ); 
	_eta_mob   = MathLib::LocalVector::Zero( _n_eta_mob ); ; 
	_eta_immob = MathLib::LocalVector::Zero( _n_eta_immob ); ; 
	_xi_mob    = MathLib::LocalVector::Zero( _n_xi_mob ); 
	_xi_immob  = MathLib::LocalVector::Zero( _n_xi_immob ); 
}

ReductionKinNodeInfo::~ReductionKinNodeInfo()
{

}

void ReductionKinNodeInfo::set_comp_conc( size_t comp_idx, double val )
{
	_Comp_Conc(comp_idx) = val; 
}

void ReductionKinNodeInfo::transform(void)
{
	_ReductionKin->Conc2EtaXi( _Comp_Conc, _eta_mob, _eta_immob, _xi_mob, _xi_immob ); 
}

