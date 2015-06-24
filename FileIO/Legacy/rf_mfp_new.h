/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_mfp_new.h
 *
 * Created on 2004-08-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: MAT-FP
   Task: class implementation
   Programing:
   08/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_mfp_new_INC
#define rf_mfp_new_INC

#include <iostream>
#include <string>
#include <vector>

namespace ogs5
{

class CFluidProperties
{
    friend bool MFPRead(const std::string&, std::vector<CFluidProperties*>&);

public:
    CFluidProperties(void);
    ~CFluidProperties(void);
    std::ios::pos_type Read(std::ifstream*);

public:
    bool cal_gravity;                     //YD/WW
    double molar_mass;
    /**
     * density
     */
    double rho_0;
    /**
     * density deviated with respect to the pressure
     */
    double drho_dp;
    /**
     * density deviated with respect to the temperature
     */
    double drho_dT;
    /**
     * density deviated with respect to the concentration
     */
    double drho_dC;

    double diffusion; /*SB:2p */

    // Viscosity
    double my_0;
    double dmy_dp;

    // Thermal properties
    double specific_heat_capacity;
    double heat_conductivity;

    // State variables
    double p_0;
    /**
     * state variable: reference temperature
     */
    double T_0;
    double C_0;

    // Chemical properties
    double T_Latent1, T_Latent2, latent_heat;

    std::string name;
    std::string fluid_name;               //NB4801

    // compressibility
    int compressibility_model_pressure;   //NB
    int compressibility_model_temperature; //NB
    int compressibility_pressure;         //NB
    int compressibility_temperature;      //NB
    int JTC;      //NB

    // Density
    int density_model;

    // Viscosity
    int viscosity_model;
    int heat_conductivity_model;
    int heat_capacity_model;              //YD, shifted to public JOD
    int diffusion_model;                  /* SB:p2 */

    int heat_phase_change_curve;
    std::vector<std::string>density_pcs_name_vector;
    std::vector<std::string>viscosity_pcs_name_vector;
    std::vector<std::string>specific_heat_capacity_pcs_name_vector;
    std::vector<std::string>heat_conductivity_pcs_name_vector;
    std::vector<std::string>enthalpy_pcs_name_vector;
    size_t phase;
};

bool MFPRead(const std::string&, std::vector<CFluidProperties*>&);
#define MFP_FILE_EXTENSION ".mfp"
}

#endif
