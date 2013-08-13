/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.h
 *
 * Created on 2013-03-28 by Haibing Shao
 */
#ifndef CHEM_EQ_REACT_SYS_H
#define CHEM_EQ_REACT_SYS_H

#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionEq.h"
#include "BaseLib/OrderedMap.h"

namespace ogsChem
{

class chemEqReactSys
{
public:
	/**
      * constructor of the class
      */
	chemEqReactSys(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                   std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);
	
	/**
      * destructor of the class
      */
	~chemEqReactSys(void); 

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 
	
	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I; }; 

    /**
      * get the number of mobile components
      */
	size_t get_n_Comp_mob(void) {return _I_mob; }; 

    /**
      * calculate the total mass using 
      * the concentration vector of basis, and
      * the concentration vector of secondary variables
      */
    void calc_tot_mass(LocalVector & vec_conc_basis, 
                       LocalVector & vec_conc_second, 
                       LocalVector & vec_tot_mass); 
    /**
      * calculate the residual of the reaction system using 
      * the concentration vector of all components and
      * the total concentration constrain of the basis species
      */
    void calc_residual(LocalVector & vec_unknowns, 
                       LocalVector & vec_tot_mass_constrain, 
                       LocalVector & vec_residual);
    
    /**
      * calculate the Jacobi of the reaction system analytically using 
      * the concentration vector of all components, 
      * the total concentration constrain of the basis species, 
      * and an already calculated residual vector
      */
    void calc_Jacobi(LocalVector & vec_unknowns,
                     LocalVector & vec_tot_mass_constrain,
                     LocalVector & vec_res_base);

    /**
      * solve the equilibrium reaction system 
      * using the Newton iterations, input values are
      * the vector of concentration values
      * the convergence tolerance, 
      * the maximum number of iterations. 
      * and the status of convergence
      */
    void solve_EqSys_Newton(LocalVector & vec_conc, 
                            size_t & result, 
                            size_t & node_idx, 
                            double iter_tol,
                            double rel_tol,
                            double max_iter); 

    /**
      * solve the system J*dx = -b
      *
      * inputs are: 
      *    - idx_node is the index of node
      *    - J is Jacobi marix
      *    - b is residual vector
      * 
      * if rank of J is equal to number of b and dx, 
      * using direct solver provided by eigen
      * if rank of J is smaller than number of b, 
      * using minimization method
      * delta_x is the returned result. 
      */
    void Min_solv(size_t      & idx_node, 
                  LocalMatrix & J, 
                  LocalVector & b, 
                  LocalVector & delta_x); 
	
private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

	/**
      * a list of all kinetic reactions
      */
	std::vector<ogsChem::chemReactionEq*> & _list_eq_reactions; 

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

    /**
      * stoichiometric matrix from the input reactions
      */	
	LocalMatrix _matStoi_input;

    /**
      * vector of natural log K values
      */	
    LocalVector _vec_lnK; 
	/**
      * _I is the number of components and _J is the number of reactions
      */
	size_t _I, _J; 

	/**
      * number of mobile, sorption and mineral reactions
      */
	size_t _J_mob, _J_sorp, _J_min;

    /**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sec_mob, _I_sec_sorp, _I_sec_min;

    /**
      * number of basis components
      */
	size_t _I_basis;

    /**
      * number of secondary components
      */
	size_t _I_second;

    /**
      * vector of system residual
      * its size is equal to the number of components
      */
    ogsChem::LocalVector _vec_res;

    /**
      * vector of staturation index
      * size of this vector is equal to the amount of minerals
      * if AI = 1, meaning the mineral is saturated and present
      * if AI = 0, meaning the mineral is under-saturated and not present
      */
    ogsChem::LocalVector _AI;

    /**
      * matrix of system Jacobi
      * its size is equal to the n_components by n_components
      */
    ogsChem::LocalMatrix _mat_Jacobi; 

	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);

    /**
      * reading in the logK values oaf each equilibrium reactions
      */
    void read_logK(std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);

	/**
      * count how many mobile, sorption and mineral components
      */
	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp); 

    /**
      * count how many mobile, sorption and mineral reactions
      */
    void countReactions(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, std::vector<ogsChem::chemReactionEq*> & list_eq_reactions);
    
    /**
      * increment the unknown with a damping factor
      * to prevent the concentrations falling into negative values
      */    
    void increment_unknown(LocalVector & x_old, 
                           LocalVector & delta_x, 
                           LocalVector & x_new);
   
    /**
      * update the concentration of minerals
      */
    void update_minerals(LocalVector & vec_unknowns, LocalVector & mass_constrain);
    
    /**
      * calcuate one particular mineral concentration 
      * by the amount of basis concentration and total mass constrain
      */
    double cal_cbarmin_by_total_mass(size_t idx_min, LocalVector & c_basis, LocalVector & tot_mass);

    /**
      * update all concentrations based on p
      */
    void update_concentations(LocalVector & vec_unknowns, LocalVector & vec_concentrations);
};

}

#endif
