/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReactionEqSorp.h
 *
 * Created on 2013-03-19 by Haibing Shao
 */
#include "chemReactionEq.h"
#include "BaseLib/OrderedMap.h" 

#ifndef CHEM_REACTION_EQ_SORP_H
#define CHEM_REACTION_EQ_SORP_H

namespace ogsChem
{

class chemReactionEqSorp : public chemReactionEq
{
public: 
	/**
      * constructor and destructor
      */
    chemReactionEqSorp(void) : chemReactionEq()
    {
        _eqReactType = ogsChem::SORP_EQ_REACT;
    };
	virtual ~chemReactionEqSorp(void){};

    /**
      * the realization of eval function
      * the returned value is the residual of the reaction. 
      * if equilibrium is reached, then redisual is zero. 
      */
    double eval(ogsChem::LocalVector & vec_log_Conc)
    {
        // TODO
        return 0.0;
    };

};

}  // end of namespace

#endif
