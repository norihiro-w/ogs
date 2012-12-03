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

#ifndef REDUCTION_KIN_NODE_INFO_H
#define REDUCTION_KIN_NODE_INFO_H

#include "ChemLib/chemReductionKin.h"

class ReductionKinNodeInfo
{
public: 
    /**
      * constructor
      */ 
	ReductionKinNodeInfo(size_t node_id, 
		                 size_t n_comp, 
	                     size_t n_eta_mob, 
						 size_t n_eta_immob, 
						 size_t n_xi_mob,
						 size_t n_xi_immob,
						 ogsChem::chemReductionKin* ReductionKin); 
    /**
      * destructor
      */ 
	~ReductionKinNodeInfo(); 

    /**
      * set the value of particular component
      */ 
	void set_comp_conc( size_t comp_idx, double val); 

    /**
      * return the node ID
      */ 
	size_t get_node_id(void) { return _node_id; }

    /**
      * return the value in eta_mob vector
      */ 
	double get_eta_mob_value( size_t eta_mob_idx ) { return _eta_mob(eta_mob_idx); }

    /**
      * return the value in eta_immob vector
      */ 
    double get_eta_immob_value( size_t eta_immob_idx ) { return _eta_immob(eta_immob_idx); }

    /**
      * return the value in xi_mob vector
      */ 
	double get_xi_mob_value( size_t xi_mob_idx ) { return _xi_mob(xi_mob_idx); }

    /**
      * return the value in xi_immob vector
      */ 
	double get_xi_immob_value( size_t xi_immob_idx ) { return _xi_immob(xi_immob_idx); }

    /**
      * transform from concentrations to eta and xi values
      */ 
	void transform(void); 

private: 
    /**
      * node ID
      */ 
	const size_t _node_id; 

    /**
      * number of components
      */ 
	const size_t _n_comp; 

    /**
      * number of eta_mob
      */ 
	const size_t _n_eta_mob; 

    /**
      * number of eta_immob
      */ 
	const size_t _n_eta_immob;

    /**
      * number of xi_mob
      */
	const size_t _n_xi_mob; 

    /**
      * number of xi_immob
      */ 
	const size_t _n_xi_immob; 

    /**
      * component concentrations vector
      */ 
	MathLib::LocalVector _Comp_Conc; 

    /**
      * eta_mob vector
      */ 
	MathLib::LocalVector _eta_mob; 

    /**
      * eta_immob vector
      */ 
	MathLib::LocalVector _eta_immob; 

    /**
      * xi_mob vector
      */ 
	MathLib::LocalVector _xi_mob; 

    /**
      * xi_immob vector
      */ 
	MathLib::LocalVector _xi_immob; 

    /**
      * pointer to the ReductionKin class
      */ 
	ogsChem::chemReductionKin* _ReductionKin;
}; 

#endif  // end of ifndef