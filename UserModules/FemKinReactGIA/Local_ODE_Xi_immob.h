/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Local_ODE_Xi_immob.h
 *
 * Created on 2012-10-08 by Haibing Shao
 */
 
#ifndef LOCAL_ODE_XI_IMMOB_H
#define LOCAL_ODE_XI_IMMOB_H

#include "ChemLib/chemReductionKin.h"

class Local_ODE_Xi_immob
{
public: 

	/**
      * constructor
      */ 
    Local_ODE_Xi_immob(ogsChem::chemReductionKin* reductionKin)
		: _reductionKin(reductionKin)
    {
		_vec_eta_mob      = MathLib::LocalVector::Zero( _reductionKin->get_n_eta_mob() ); 
	    _vec_eta_immob    = MathLib::LocalVector::Zero( _reductionKin->get_n_eta_immob() ); 
	    _vec_xi_mob       = MathLib::LocalVector::Zero( _reductionKin->get_n_xi_mob() ); 
	    _vec_xi_immob     = MathLib::LocalVector::Zero( _reductionKin->get_n_xi_immob() ); 
		_vec_dxi_immob_dt = MathLib::LocalVector::Zero( _reductionKin->get_n_xi_immob() ); 	 
	}
	 
    /**
      * destructor
      */ 
	~Local_ODE_Xi_immob()
    {
        _reductionKin = NULL; 
    }

    /**
      * update the eta and xi values
      */ 
	void update_eta_xi( MathLib::LocalVector & vec_eta_mob, 
		                MathLib::LocalVector & vec_eta_immob, 
						MathLib::LocalVector & vec_xi_mob, 
						MathLib::LocalVector & vec_xi_immob  )
	{
		_vec_eta_mob   = vec_eta_mob; 
		_vec_eta_immob = vec_eta_immob; 
	    _vec_xi_mob    = vec_xi_mob; 
	    _vec_xi_immob  = vec_xi_immob;
	}

    /**
      * evaluate the change of xi_immob over time
      */ 
	MathLib::LocalVector operator() (double time, MathLib::LocalVector vec_xi_immob )
	{
	    _vec_xi_immob = vec_xi_immob; 
		this->_reductionKin->Calc_Xi_immob_Rate( _vec_eta_mob, _vec_eta_immob, _vec_xi_mob, _vec_xi_immob, _vec_dxi_immob_dt); 
		return _vec_dxi_immob_dt; 
	}

private: 

    /**
      * pointer to reduction scheme class
      */ 
	ogsChem::chemReductionKin* _reductionKin; 
	
    /**
      * local vector of eta_mob
      */ 
    MathLib::LocalVector _vec_eta_mob; 

    /**
      * local vector of eta_immob
      */ 
    MathLib::LocalVector _vec_eta_immob; 
    
    /**
      * local vector of xi_mob
      */ 
    MathLib::LocalVector _vec_xi_mob; 

    /**
      * local vector of xi_immob
      */ 
	MathLib::LocalVector _vec_xi_immob;

    /**
      * local vector of dxi_immob_dt
      */ 
	MathLib::LocalVector _vec_dxi_immob_dt; 
 
}; 
 
 #endif