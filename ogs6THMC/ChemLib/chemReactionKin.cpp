/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReaction.h
 *
 * Created on 2012-09-05 by Haibing Shao
 */
 
#include "chemReactionKin.h"
#include "logog.hpp"

namespace ogsChem
{

chemReactionKin::chemReactionKin()
{
	// set initial values
	_rate = 0.0;
	_rate_constant = 0.0; 
	_rate_constant_order = 1.0; 
    _decay_rate = 0.0; 
}

chemReactionKin::~chemReactionKin(void)
{
}

void chemReactionKin::readReactionStr(std::string & reaction_str)
{
	// TODO
	std::string tmp_str;
}

void chemReactionKin::calcReactionRate(ogsChem::LocalVector & vec_Comp_Conc)
{
	if ( this->_kinReactType == ogsChem::Monod )
		this->_rate = calcReactionRateMonod(vec_Comp_Conc); 
    else
        this->_rate = 0.0; 	
}

double chemReactionKin::calcReactionRateMonod(ogsChem::LocalVector & vec_Comp_Conc)
{
	size_t i, comp_idx; 
	double rate = 1.0; 
	double conc, c_bio;

    // get the concentration of biomass
    c_bio = vec_Comp_Conc(this->_idx_bacteria); 

	// loop over all the monod term
	for (i=0; i<this->_vec_Monod_Comps_Idx.size(); i++ )
	{
		comp_idx = _vec_Monod_Comps_Idx[i]; 
		conc = vec_Comp_Conc( comp_idx ); 
        // non-negative protection --------------------
        // if (conc < 0.0)
        //     conc = 0.0; 
        // end of non-negative protection -------------
		// rate *= [ c / ( c + k_c) ]^order
		rate *= pow( conc / (conc + this->_vec_Monod_Comps_Conc[i]), this->_vec_Monod_Comps_order[i] ); 
	}  // end of for i
	
	// rate constant
	rate *= _rate_constant; // mu_max
	// rate order
	// TODO
	// times c_bio
	rate *= c_bio;
    // decay term
    rate -= _decay_rate * c_bio; 
	// retrun rate
	return rate; 
}


void chemReactionKin::readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                                      ogs5::CKinReact* KRC_reaction)
{
	// get the list of components 
	_vecCompNames = KRC_reaction->reactionpartner; 

	// find the components based on component names
	size_t i; 
	std::string str_comp; 
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator comp_iterator; 
	for (i=0; i<_vecCompNames.size(); i++)
	{
		str_comp = _vecCompNames[i];
		comp_iterator = list_chemComp.find(str_comp);

		if ( comp_iterator != list_chemComp.end() )
		{
			ogsChem::ChemComp* component = comp_iterator->second; 
			this->_vecComponents.push_back(component);		
		}  // end of if	
	}  // end of for i

    // copy the stoichiometric vector 
	_vecStoi = KRC_reaction->stochmet;

	// read the rate parameters
	if ( KRC_reaction->getType() == "monod" )
	{
		this->_kinReactType = ogsChem::Monod; 

        // find the bacteria component
        std::string str_bacteria = KRC_reaction->bacteria_name; 
        this->_idx_bacteria = list_chemComp.find( str_bacteria )->second->getIndex(); 

        if ( this->_idx_bacteria >= list_chemComp.size() )
        {
            ERR("When reading Monod reaction, bacteria name not found! ");
            exit(1);
        }

		// loop over the monod term, 
		for (size_t i=0; i < KRC_reaction->monod.size(); i++ )
		{
			size_t monod_comp_idx; 
			double monod_comp_conc(0.0), monod_term_order(1.0); 
			// read the corresponding components
			std::string comp_str = KRC_reaction->monod[i]->species; 
			// find this component
            // get its index value in all component list
            monod_comp_idx = list_chemComp.find( comp_str )->second->getIndex();
			
            if ( monod_comp_idx >= list_chemComp.size() )
            {
                ERR("When reading Monod terms, component name not found! ");
                exit(1);
            }

			// read the monod component concentration
			monod_comp_conc = KRC_reaction->monod[i]->concentration; 
			// read the monod term order
			monod_term_order = KRC_reaction->monod[i]->order; 

			this->_vec_Monod_Comps_Idx.push_back(   monod_comp_idx   );
			this->_vec_Monod_Comps_Conc.push_back(  monod_comp_conc  ); 
			this->_vec_Monod_Comps_order.push_back( monod_term_order ); 
		}  // end of for i
	}  // end of if KRC_reaction
	else
	{
		this->_kinReactType = ogsChem::NoType; 
	}

	// read the rate constant
	this->_rate_constant = KRC_reaction->rateconstant; 
	// read the rate order
	this->_rate_constant_order = KRC_reaction->rateorder;
    // read the decay rate
    this->_decay_rate = KRC_reaction->decay_rate; 
}


}  // end of namespace