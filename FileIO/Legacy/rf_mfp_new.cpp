/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_mfp_new.cpp
 *
 * Created on 2004-08-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: MFP Fluid Properties
   Task:
   Programing:
   08/2004 OK Implementation
   last modified:
**************************************************************************/
#include "rf_mfp_new.h"

#include "makros.h"
// C++ STL
//#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cstdlib>

#include "Ogs5FileTools.h"

//
using namespace std;

namespace ogs5
{

/* Umrechnungen SI - Amerikanisches System */
//WW #include "steam67.h"
//#define PSI2PA 6895.
//#define PA2PSI 1.4503263234227701232777374909355e-4
//#define GAS_CONSTANT    8314.41
//#define COMP_MOL_MASS_AIR    28.96
//#define COMP_MOL_MASS_WATER  18.016
//#define GAS_CONSTANT_V  461.5                     //WW
//#define T_KILVIN_ZERO  273.15                     //AKS


/**************************************************************************
   FEMLib-Method:
   Task: OBJ constructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CFluidProperties::CFluidProperties()
//: name ("WATER")
{
//    phase = 0;
//    // Density
//    density_model = 1;
//    rho_0 = 1000.;
//    drho_dp = 0.;
//    drho_dT = 0.;
//    drho_dC = 0.;
//    // Viscosity
//    viscosity_model = 1;
//    my_0 = 1e-3;
//    dmy_dp = 0.;
//    dmy_dT = 0.;
//    dmy_dC = 0.;
//    // Thermal properties
//    heat_capacity_model = 1;
//    specific_heat_capacity = 4680.;       //CMCD we should give this as SHC not HC GeoSys 4 9/2004
//    heat_conductivity_model = 1;
//    heat_conductivity = 0.6;
//    // Electrical properties
//    // Chemical properties
//    diffusion_model = 1;
//    diffusion = 2.13e-6;
//    // State variables
//    p_0 = 101325.;
//    T_0 = 293.;
//    C_0 = 0.;
//    Z = 1.;
//    cal_gravity = true;
//    // Data
//    mode = 0;                             // Gauss point values
//    Fem_Ele_Std = NULL;
//    // WW
//    molar_mass = COMP_MOL_MASS_AIR;

    cal_gravity = true;
    molar_mass = .0;
    rho_0 = .0;
    drho_dp = .0;
    drho_dT = .0;
    drho_dC = .0;
    diffusion = .0;
    my_0 = .0;
    dmy_dp = .0;
    specific_heat_capacity = .0;
    heat_conductivity = .0;
    p_0 = .0;
    T_0 = .0;
    C_0 = .0;
    T_Latent1 = .0;
    T_Latent2 = .0;
    latent_heat = .0;
    compressibility_model_pressure = 0;
    compressibility_model_temperature = 0;
    compressibility_pressure = .0;
    compressibility_temperature = .0;
    JTC = 0;
    density_model = 0;
    viscosity_model = 0;
    heat_conductivity_model = 0;
    heat_capacity_model = 0;
    diffusion_model = 0;
    heat_phase_change_curve = 0;
    phase = 0;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ deconstructor
   Programing:
   05/2010 OK/AKS Implementation
**************************************************************************/
CFluidProperties::~CFluidProperties(void)
{
//    for (int i = 0; i < (int)component_vector.size(); i++ )
//        component_vector[i] = NULL;
//    component_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   11/2004 SB string streaming
**************************************************************************/
std::ios::pos_type CFluidProperties::Read(std::ifstream* mfp_file)
{


    std::string sub_line;
    std::string line_string;
    std::string delimiter(" ");
    bool new_keyword = false;
    std::string hash("#");
    std::ios::pos_type position;
    std::string sub_string;
    //WW bool new_subkeyword = false;
    std::string dollar("$");
    std::string delimiter_type(":");
    std::stringstream in;
    //========================================================================
    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        //WW new_subkeyword = false;
        position = mfp_file->tellg();
        //SB    mfp_file->getline(buffer,MAX_ZEILE);
        //SB    line_string = buffer;
        line_string = readNonBlankLineFromInputStream(*mfp_file);
        if(line_string.size() < 1)
            break;
        if(line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            break;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$FLUID_TYPE") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> name;   //sub_line
            in.clear();
            continue;
        }
        //....................................................................
        // NB 4.8.01
        if(line_string.find("$FLUID_NAME") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> fluid_name; //sub_line
            in.clear();
            continue;
        }
        //....................................................................
        // NB Oct-2009
        if(line_string.find("$COMPRESSIBILITY") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> compressibility_model_pressure; //sub_line 1 for first phase
            in >> compressibility_pressure; //sub_line 1
            in.clear();
            in >> compressibility_model_temperature; //sub_line 2 for second phase
            in >> compressibility_temperature; //sub_line 2
            in >> JTC;
            in.clear();

            // available models see CFluidProperties::drhodP and CFluidProperties::drhodT
            // 0 incompressible
            // 1 constant slope
            // 2 slope from fct_table
            // 3 difference quotient
            // 4 analytical derivation

            in.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$DAT_TYPE") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> name;   //sub_line
            in.clear();
            continue;
        }
        //YD/WW subkeyword found
        if(line_string.find("$NON_GRAVITY") != string::npos)
        {
            cal_gravity = false;
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$DENSITY") != string::npos)
        {
            //WW new_subkeyword = false;
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> density_model;
            // TF - _rho_fct_name is only used for writing it back to the file
//            if(density_model == 0) // rho = f(x)
//                in >> _rho_fct_name;

            if(density_model == 1) // rho = const
                in >> rho_0;

            if(density_model == 2) // rho(p) = rho_0*(1+beta_p*(p-p_0))
            {
                in >> rho_0;
                in >> p_0;
                in >> drho_dp;
                density_pcs_name_vector.push_back("PRESSURE1");
            }
            if(density_model == 3) // rho(C) = rho_0*(1+beta_C*(C-C_0))
            {
                in >> rho_0;
                in >> C_0;
                in >> drho_dC;
                //        density_pcs_name_vector.push_back("CONCENTRATION1");
                // PCH
                density_pcs_name_vector.push_back("Isochlor");
            }
            if(density_model == 4) // rho(T) = rho_0*(1+beta_T*(T-T_0))
            {
                in >> rho_0;
                in >> T_0;
                in >> drho_dT;
                density_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(density_model == 5) // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
            {
                in >> rho_0;
                in >> C_0;
                in >> drho_dC;
                in >> T_0;
                in >> drho_dT;
                density_pcs_name_vector.push_back("CONCENTRATION1");
                density_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(density_model == 6) // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
            {
                in >> rho_0;
                in >> p_0;
                in >> drho_dp;
                in >> T_0;
                in >> drho_dT;
                density_pcs_name_vector.push_back("PRESSURE1");
                density_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(density_model == 7) // rho(p,p_v,T)
            {
                // no input data required
            }
            if(density_model == 8) // rho(p,T,C)
            {
                in >> C_0;
                density_pcs_name_vector.push_back("PRESSURE1");
                density_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(density_model == 9) //WW
                // Molar mass
                in >> molar_mass;

            if((density_model == 10) //NB 4.8.01  read density from a rho-P-T table
               || (density_model == 11) //NB 4.9.05  Peng-Robinson Equation of State
               || (density_model == 12) //NB 4.9.05  Redlich-Kwong Equation of State
               || (density_model == 13) //NB JUN 09  Fundamental equation
               || (density_model == 14)) //AKS MAY 10  Extended Ideal gas Eq. based on Super compressibility factor
            {
                std::string arg1,arg2,arg3;
                in >> arg1 >> arg2 >> arg3; //get up to three arguments for density model

                if (isdigit(arg1[0]) != 0) // first argument is reference temperature
                {
                    T_0 = atof(arg1.c_str());
                    arg1 = arg2;
                    arg2 = arg3;
                }

                if (arg1.length() == 0) // if no arguments are given use standard
                {
                    arg1 = "PRESSURE1";
                    arg2 = "TEMPERATURE1";
                }
                else if (arg2.length() == 0) // if only PRESSURE argument is given
                    arg2 = "TEMPERATURE1";

                density_pcs_name_vector.push_back(arg1);
            }
            if(density_model == 18) // BG, NB calculated node densities from the phase transition model
            {
            }

            //      mfp_file->ignore(MAX_ZEILE,'\n');
            in.clear();
            continue;
        }
        if(line_string.find("$TEMPERATURE") != string::npos) // subkeyword found 11/2010, BG, NB, DL, SB
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> T_0 >> T_0;
            in.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$VISCOSITY") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> viscosity_model;
            // TF 11/2011 - used only in read- and write-method
//            if(viscosity_model == 0) // my = fct(x)
//                in >> _my_fct_name;
            if(viscosity_model == 1) // my = const

                in >> my_0;
            if(viscosity_model == 2) // my(p) = my_0*(1+gamma_p*(p-p_0))
            {
                in >> my_0;
                in >> p_0;
                in >> dmy_dp;
                viscosity_pcs_name_vector.push_back("PRESSURE1");
            }
            if(viscosity_model == 3) // my(T), Yaws et al. (1976)

                //OK4704
                viscosity_pcs_name_vector.push_back("TEMPERATURE1");
            if(viscosity_model == 4) // my(T), ???
            {
            }
            if(viscosity_model == 5) // my(p,T), Reichenberg (1971)
            {
            }
            if(viscosity_model == 6) // my(C,T),
            {
            }
            if(viscosity_model == 7) // my(p,T,C)
            {
                in >> C_0;
                viscosity_pcs_name_vector.push_back("PRESSURE1");
                viscosity_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(viscosity_model == 9) // my(rho,T)
            {
                std::string arg1,arg2;
                in >> arg1 >> arg2; //get up to three arguments for density model

                if (arg1.length() == 0) // if no arguments are given use standard
                {
                    arg1 = "PRESSURE1";
                    arg2 = "TEMPERATURE1";
                }
                else if (arg2.length() == 0) // if only PRESSURE argument is given

                    arg2 = "TEMPERATURE1";

                viscosity_pcs_name_vector.push_back(arg1);
            }
            if(viscosity_model == 18) // BG, NB calculated node viscosities from the phase transition model
            {
            }

            //    mfp_file->ignore(MAX_ZEILE,'\n');
            in.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$SPECIFIC_HEAT_CAPACITY") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> heat_capacity_model;
            // TF 11/2011 - used only in read- and write-method
//            if(heat_capacity_model == 0) // c = fct(x)
//                in >> heat_capacity_fct_name;
            if(heat_capacity_model == 1) // c = const
            {
                in >> specific_heat_capacity;
                specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
                specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(heat_capacity_model == 2) // my(p,T,C)
            {
                in >> C_0;
                specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
                specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(heat_capacity_model == 3) // YD: improved phase change
            {
                in >> T_Latent1; // Tmin for phase change
                in >> T_Latent2; // Tmax for phase change
                in >> heat_phase_change_curve;
                specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
                specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
                specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
                enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(heat_capacity_model == 4) // YD: improved phase change, function
            {
                in >> T_Latent1; // Tmin for phase change
                in >> T_Latent2; // Tmax for phase change
                in >> specific_heat_capacity; // ^c
                in >> latent_heat; // L
                specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
                specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
                specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
                enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
            }
            in.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$HEAT_CONDUCTIVITY") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> heat_conductivity_model;
            // TF 11/2011 - used only in read- and write-method
//            if(heat_conductivity_model == 0) // my = fct(x)
//                in >> heat_conductivity_fct_name;
            if(heat_conductivity_model == 1) // my = const

                in >> heat_conductivity;
            if(heat_conductivity_model == 2) // my = f(p,T,C)
            {
                in >> C_0;
                heat_conductivity_pcs_name_vector.push_back("PRESSURE1");
                heat_conductivity_pcs_name_vector.push_back("TEMPERATURE1");
            }
            if(heat_conductivity_model == 3) // my = f(p,T) NB
            {
                std::string arg1,arg2;
                in >> arg1 >> arg2; //get up to three arguments for density model

                if (arg1.length() == 0) // if no arguments are given use standard
                {
                    arg1 = "PRESSURE1";
                    arg2 = "TEMPERATURE1";
                }
                else if (arg2.length() == 0) // if only PRESSURE argument is given

                    arg2 = "TEMPERATURE1";

                heat_conductivity_pcs_name_vector.push_back(arg1);
                heat_conductivity_pcs_name_vector.push_back(arg2);
            }
            in.clear();   //OK
            continue;
        }
        // subkeyword found
        if(line_string.find("$PHASE_DIFFUSION") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mfp_file));
            in >> diffusion_model;
            // TF 11/2011 - only read - not used
//            if(diffusion_model == 0) // D = fct(x)
//                in >> dif_fct_name;
            if(diffusion_model == 1) // D = const //MX

                in >> diffusion;
            in.clear();
            continue;
        }
    }
    return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: Master read function
   Programing:
   08/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool MFPRead(const std::string &file_base_name, std::vector<CFluidProperties*> &mfp_vector)
{
    //----------------------------------------------------------------------
    //OK  MFPDelete();
    //----------------------------------------------------------------------
    CFluidProperties* m_mfp = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    std::string mfp_file_name = file_base_name + MFP_FILE_EXTENSION;
    std::ifstream mfp_file (mfp_file_name.data(),std::ios::in);
    if (!mfp_file.good())
        return false;
    mfp_file.seekg(0L,std::ios::beg);
    //========================================================================
    // Keyword loop
    std::cout << "MFPRead ... " << std::flush;
    while (!mfp_file.eof())
    {
        mfp_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos) {
            std::cout << "done, read " << mfp_vector.size() << " fluid properties" <<
            std::endl;
           return true;
        }
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#FLUID_PROPERTIES") != std::string::npos)
        {
            m_mfp = new CFluidProperties();
            position = m_mfp->Read(&mfp_file);
            m_mfp->phase = (int)mfp_vector.size(); //OK4108
            mfp_vector.push_back(m_mfp);
            mfp_file.seekg(position,std::ios::beg);
        }                         // keyword found
    }                                     // eof
    //========================================================================
    // Configuration
    int i;
    int no_fluids = (int)mfp_vector.size();
    if(no_fluids == 1)
    {
        m_mfp = mfp_vector[0];
        m_mfp->phase = 0;
    }
    else if(no_fluids == 2)
        for(i = 0; i < no_fluids; i++)
        {
            m_mfp = mfp_vector[i];
            if(m_mfp->name.find("GAS") != string::npos)
                m_mfp->phase = 0;
            else
                m_mfp->phase = 1;
        }
    //----------------------------------------------------------------------
    // Test
    if(mfp_vector.size() == 0)
    {
        std::cout << "Error in MFPRead: no MFP data" << std::endl;
        abort();
    }
    //----------------------------------------------------------------------
    return true;
}

}

