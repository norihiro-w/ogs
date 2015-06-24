/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_mmp_new.h
 *
 * Created on 2004-01-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib-Object: MAT-MP
   Task: MediumProperties
   Programing:
   01/2004 OK Implementation
**************************************************************************/

#ifndef rfmat_mp_new_INC
#define rfmat_mp_new_INC
/* Schutz gegen mehrfaches Einfuegen */

// C++ STL
#include <list>
#include <string>
#include <vector>
#include <iostream>

namespace ogs5
{

class CMediumProperties
{
public:
    // Methods
    CMediumProperties(void);              // constructor
    ~CMediumProperties(void);             // destructor
    std::ios::pos_type Read(std::ifstream*);

public:
    // Permeability
    // Permeabilty stress corrector WW
    double* c_coefficient;
    unsigned geo_dimension;
    int permeability_stress_mode;
    //
    // Properties
    // PCS
    std::string pcs_type_name;            //YD
    std::vector<std::string>pcs_name_vector;
    std::vector<std::string> porosity_pcs_name_vector;

    //GEO
    std::string geo_name;
    std::vector<std::string>geo_name_vector; //OK
    double geo_area;
    std::string geo_area_file;            //OK

    double density;
    std::string name;
    int number;
    int porosity_model;                   // porosity
    int porosity_curve;
    double porosity_model_values[15];
    //double porosity;
    double KC_porosity_initial;           // HS 11.2008
    double KC_permeability_initial;       // HS 11.2008
    std::string porosity_file;            //OK/MB
    int tortuosity_model;
    double tortuosity_model_values[10];
    double tortuosity;
    int flowlinearity_model;
    double flowlinearity_model_values[10];
    int storage_model;                    // storativity
    double storage_model_values[10];
    //double storage;
    int conductivity_model;
    double conductivity;
    int unconfined_flow_group;
    int permeability_model;               // permeability
    double permeability;
    double permeability_tensor[9];       
    std::string permeability_tensor_type_name;
    std::string tortuosity_tensor_type_name;
    int permeability_tensor_type;
    int tortuosity_tensor_type;

    int permeability_pressure_model;
    double permeability_pressure_model_values[10];
    double permeability_pressure_rel;
    int permeability_strain_model;        //WX: permeability function strain model. 05.2010
    int permeability_strain_model_value[3]; //WX:permeability fuction strain model value. 05.2010
    //
    // Relative permeability (JT)
    int permeability_saturation_model[2];
    double minimum_relative_permeability;
    double residual_saturation[2];
    double maximum_saturation[2];
    double saturation_exponent[2];
    double perm_saturation_value[2];
    //
    std::string permeability_file;        //SB //OK/MB string permeability_dis_type_file;
    std::string tortuosity_file;          // PCH
    bool entry_pressure_conversion;        //JT
    int capillary_pressure_model;
    int permeability_porosity_model;
    double permeability_porosity_model_values[10];
    double storativity;
    double capillary_pressure_values[5];        //JT2012
    int heat_capacity_model;
    double heat_capacity;                 // thermal properties
    int mass_dispersion_model;
    double mass_dispersion_longitudinal;
    double mass_dispersion_transverse;
    double lgpn;                          //local grid peclet number
    int heat_dispersion_model;
    double heat_dispersion_longitudinal;
    double heat_dispersion_transverse;
    int heat_conductivity_model;
    double heat_conductivity_tensor[9];
    int fct_number;                       // functions
    int heat_diffusion_model;
    int evaporation;                      // if it is 647 then evaporation ON, else OFF: and corresponding heat loss will compensated by heat ST
    double heatflux;
    double vaporfraction;
    //aux
    int m_color[3];
    bool selected;
    int mode;
    // surface water
    // JOD
    double friction_coefficient, friction_exp_slope, friction_exp_depth;
    // JOD
    double overland_width, rill_height, rill_epsilon;
    bool channel;
    double argument;                      //OK
    //Dual Richards transfer coefficient  YD
    double transfer_coefficient;
    //double unsaturated_hydraulic_conductivity;
    double specific_storage;
    int vol_mat_model;                    // CB
    double vol_mat;                       //SB
    int vol_bio_model;                    // CB
    double vol_bio;                       //SB
    double foc;                           // organic carbon content
    bool is_fracture; //NW
};


extern bool MMPRead(const std::string&, std::vector<CMediumProperties*> &mmp_vector);

#define MMP_FILE_EXTENSION ".mmp"

}

#endif
