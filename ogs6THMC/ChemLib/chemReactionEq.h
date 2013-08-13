/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReactionEq.h
 *
 * Created on 2013-03-19 by Haibing Shao
 */
#include "chemReaction.h"
#include "ogsFileIO/FemIO/ogs5/rf_kinreact.h"
#include "BaseLib/OrderedMap.h" 

#ifndef CHEM_REACTION_EQ_H
#define CHEM_REACTION_EQ_H

namespace ogsChem
{

class chemReactionEq : public chemReaction
{
public: 
	/**
      * constructor and destructor
      */
	chemReactionEq();
	virtual ~chemReactionEq(void);

	/**
      * return the constant of equilibrium constant
      */
	double get_K(void) { return _eq_const_K; }

    /**
      * return the rate of lnK. 
      */
	double get_ln_K(void) { return _ln_K; }

    /**
      * return the type of equilibrium reaction. 
      */
    EqReactType get_type(void) { return _eqReactType;}

	/**
      * override reading str function
      */
	void readReactionStr(std::string & reaction_str); 

	/**
      * reading one reaction from KRC data structure
      */
	void readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                         ogs5::CKinReact* KRC_reaction); 

    /**
      * evaluate the residual with given log concentrations
      * input is a vector of log concentrations; 
      * the returned value is the residual of the reaction. 
      * if equilibrium is reached, then redisual is zero. 
      */
    virtual double eval(ogsChem::LocalVector & vec_log_Conc)
    {return 0.0;};

protected:
	/**
      * type of kinetic reaction
      */
	EqReactType _eqReactType; 

private:
	/**
      * the equilibrium constant K
      */
	double _eq_const_K; 

    /**
      * natural log of K
      */
	double _ln_K;

};

}  // end of namespace

#endif
