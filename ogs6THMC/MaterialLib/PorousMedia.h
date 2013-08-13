/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "IMedium.h"
#include "NumLib/Function/FunctionLinear.h"
#include "logog.hpp"


namespace NumLib
{
class ITXFunction;
class FunctionLinear1D;
}

namespace MaterialLib
{

struct PorousMedia : public IMedium
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    NumLib::ITXFunction* geo_area;
    NumLib::ITXFunction* dispersivity_long;
    NumLib::ITXFunction* dispersivity_trans; 
	NumLib::ITXFunction* res_saturation; 
    NumLib::ITXFunction* max_saturation;      
    NumLib::ITXFunction* exp_saturation;      
    NumLib::FunctionLinear1D* capp_sat_curve; 
    std::vector<size_t> perm_saturation_model;    
    std::vector<NumLib::FunctionLinear1D*> perm_saturation_curve; 
    size_t capp_sat_model;          
	double minimum_relative_permeability;  

    PorousMedia()
    {
        BaseLib::zeroObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area,
                res_saturation, 
                max_saturation, 
				exp_saturation
                );
    }

    virtual ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area,
                res_saturation,
				max_saturation,
                exp_saturation
                );
        perm_saturation_model.clear();
        BaseLib::releaseObjectsInStdVector(perm_saturation_curve);
        BaseLib::releaseObject(capp_sat_curve);

    }

    virtual MediumType::type getMediumType() const {return MediumType::PorousMedium;};

	virtual double getKrelbySw(const double Sw/*wetting saturation*/, size_t idx_phase)
    {
        double kr = 0.0; 
        bool phase_shift = false; 
        size_t model; 
        model = this->perm_saturation_model[idx_phase];
		minimum_relative_permeability =  1.0e-9; // Could be set as default at a beter place

        switch(model)
        {
	     default: // Error occurs! 
		    ERR("ERROR in PermeabilitySaturationFunction(): Unrecognized relative permeability method.\n");
            exit(0);
	        break;
	    case 0: // CURVE
			this->perm_saturation_curve[idx_phase]->eval( Sw, kr ); 
            if( kr < minimum_relative_permeability )
	        kr = minimum_relative_permeability;
            break;
        }
        return kr; 
    }

    /**
    * return the water saturation value by capillary pressure
    */
    virtual double getSwbyPc(double Pc)
    {	
        double Sw = 0.0;

        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval( Pc, Sw );
            break;
        default: 
            ERR("No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return Sw; 
    }

	virtual double getdSwdPc (double Pc)
    {
        double dSwdPc = 0.0;

        switch ( this->capp_sat_model )
        {
        case 0:  // curve value
            // get the value from curve. 
            capp_sat_curve->eval_slope( Pc, dSwdPc );
            break;
        default: 
            ERR("Error in getSwbyPc: No valid capilary pressure vs water saturation model! "); 
            exit(0); 
            break;
        } // end of switch capp_sat_model
        return dSwdPc; 
    }


};

} //end
 
