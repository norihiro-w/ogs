/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionKin.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#ifndef CHEMREDUCTIONKIN_H
#define CHEMREDUCTIONKIN_H

#include "chemconst.h"
#include "chemReactionKin.h"
#include "chemcomp.h"
#include "BaseLib/OrderedMap.h"

namespace ogsChem
{

class chemReductionKin
{
public:
	/**
      * constructor of the class
      */
	chemReductionKin(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		             std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);
	
	/**
      * destructor of the class
      */
	~chemReductionKin(void); 

    /**
      * convert concentration vector to eta and xi vector
      */
	void Conc2EtaXi(ogsChem::LocalVector &local_conc, 
		            ogsChem::LocalVector &local_eta_mob, 
					ogsChem::LocalVector &local_eta_immob, 
					ogsChem::LocalVector &local_xi_mob, 
					ogsChem::LocalVector &local_xi_immob);

	/**
      * convert eta and xi vector to concentration vector
      */
	void EtaXi2Conc(ogsChem::LocalVector &local_eta_mob, 
		            ogsChem::LocalVector &local_eta_immob, 
					ogsChem::LocalVector &local_xi_mob, 
					ogsChem::LocalVector &local_xi_immob,
					ogsChem::LocalVector &local_conc); 

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 
	
	/**
      * calculate the reaction rates of xi_mob
      */
	void Calc_Xi_mob_Rate(ogsChem::LocalVector &local_eta_mob, 
		                  ogsChem::LocalVector &local_eta_immob, 
					      ogsChem::LocalVector &local_xi_mob,
					      ogsChem::LocalVector &local_xi_immob, 
					      ogsChem::LocalVector &xi_mob_rate ); 

    /**
      * calculate the reaction rates of xi_mob
      */
	void Calc_Xi_immob_Rate(ogsChem::LocalVector &local_eta_mob, 
		                    ogsChem::LocalVector &local_eta_immob, 
					        ogsChem::LocalVector &local_xi_mob,
					        ogsChem::LocalVector &local_xi_immob, 
					        ogsChem::LocalVector &xi_immob_rate ); 

	void Calc_Xi_Rate(ogsChem::LocalVector &local_eta_mob, 
	                  ogsChem::LocalVector &local_eta_immob, 
	                  ogsChem::LocalVector &local_xi_mob,
					  ogsChem::LocalVector &local_xi_immob, 
					  ogsChem::LocalVector &xi_mob_rate, 
					  ogsChem::LocalVector &xi_immob_rate);

	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I; }; 
	
	/**
      * get the length of eta
      */
	size_t get_n_eta(void) {return _n_eta;}; 

	/**
      * get the length of eta_mob
      */
	size_t get_n_eta_mob(void) {return _n_eta_mob;}; 

	/**
      * get the length of eta_immob
      */
	size_t get_n_eta_immob(void) {return _n_eta - _n_eta_mob;}; 

	/**
      * get the length of xi
      */
	size_t get_n_xi(void) {return _n_xi; }; 

	/**
      * get the length of xi
      */
	size_t get_n_xi_mob(void) {return _n_xi_mob; }; 

	/**
      * get the length of xi
      */
	size_t get_n_xi_immob(void) {return _n_xi_immob; }; 

private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

	/**
      * a list of all kinetic reactions
      */
	std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions; 

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * sub-matrixes of S, 
      */	
	LocalMatrix _matS_1, _matS_2, _matS_1_bar, _matS_2_bar; 

	/**
      * complementary orthorgonal matrixes
      */
	LocalMatrix _mat_s_1, _mat_s_2, _matS_1_ast, _matS_2_ast, _matS_1_bar_ast, _matS_2_bar_ast; 

	/**
      * matrix for the calculation of right hand side rate transformation
      */
	LocalMatrix _matA1, _matA2; 

	/**
      * transformation matrix for the calculation between eta, xi and concentrations
      */
	LocalMatrix _mat_c_mob_2_eta_mob,
		        _mat_c_immob_2_eta_immob,   
	            _mat_c_mob_2_xi_mob, 
				_mat_c_immob_2_xi_immob;

	/**
      * transformation matrixes
      */
	// LocalMatrix _c_mob_2_eta_mob, _c_mob_2_xi_mob; 

	/**
      * the size of eta and xi vector
      */
	size_t _n_eta, _n_xi, _n_eta_mob, _n_eta_immob, _n_xi_mob, _n_xi_immob; 

	/**
      * _I is the number of components and _J is the number of reactions
      */
	size_t _I, _J; 

	/**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sorp, _I_min;

	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);

	/**
      * calculate the intemediate parameters for the reduction scheme
      */
	void update_reductionScheme(void); 

	/**
      * return the orthogonal complement of the given matrix
      */
	LocalMatrix orthcomp( LocalMatrix & inMat ); 

	/**
      * count how many mobile and immobile components
      */
	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp); 

};

}

#endif