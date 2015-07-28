/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_msp_new.cpp
 *
 * Created on 2004-08-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: MSP Solid Properties
   Task:
   Programing:
   08/2004 WW Implementation
   last modified:
**************************************************************************/

#include "rf_msp_new.h"

// C++ STL
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cfloat>

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"

// FEM-Makros
#include "makros.h"
#include "Ogs5FileTools.h"

using namespace std;

namespace ogs5
{
using namespace Math_Group;

/**************************************************************************
   FEMLib-Method:
   Task: Master read function
   Programing:
   08/2004 OK Implementation for fluid properties
   08/2004 WW Modification for solid properties
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool MSPRead(const std::string &file_base_name, std::vector<CSolidProperties*> &msp_vector)
{
    //----------------------------------------------------------------------
    //OK  MSPDelete();
    //----------------------------------------------------------------------
    CSolidProperties* m_msp = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    std::string msp_file_name = file_base_name + MSP_FILE_EXTENSION;
    std::ifstream msp_file (msp_file_name.data(),std::ios::in);
    if (!msp_file.good())
        return false;
    msp_file.seekg(0L,std::ios::beg);
    //========================================================================
    // Keyword loop
    //INFO("MSPRead ... ");
    while (!msp_file.eof())
    {
        msp_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos) {
            //INFO("done, read %d solid properties", msp_vector.size());
           return true;
        }
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#SOLID_PROPERTIES") != std::string::npos)
        {
            m_msp = new CSolidProperties();
            m_msp->file_base_name = file_base_name;
            position = m_msp->Read(&msp_file);
            msp_vector.push_back(m_msp);
            msp_file.seekg(position,std::ios::beg);
        }                         // keyword found
    }                                     // eof
    return true;
    //========================================================================
}
/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation for fuild properties
   08/2004 WW Modification for solid properties
   12/2005 WW Creep properties
**************************************************************************/
std::ios::pos_type CSolidProperties::Read(std::ifstream* msp_file)
{
    //char buffer[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::string delimiter(" ");
    bool new_keyword = false;
    std::string hash("#");
    std::ios::pos_type position;
    //  ios::pos_type position0;
    std::string sub_string;
    //WW bool new_subkeyword = false;
    std::string dollar("$");
    std::string delimiter_type(":");
    std::stringstream in_sd;

    int i = 0, Size = 0;

    //========================================================================
    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        //WW new_subkeyword = false;

        position = msp_file->tellg();
        line_string = readNonBlankLineFromInputStream(*msp_file);
        BaseLib::trim(line_string, ':'); //NW
        if(line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            break;
        }
        //....................................................................
        //NAME //OK
        //subkeyword found
        if(line_string.find("$NAME") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> name; //sub_line
            in_sd.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$SWELLING_PRESSURE_TYPE") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> SwellingPressureType;
            if(SwellingPressureType == 1 || SwellingPressureType == 2)
            {
                in_sd >> Max_SwellingPressure;
                in_sd.clear();
            }
            //10.03.2008 WW
            else if (SwellingPressureType == 3 || SwellingPressureType == 4)
            {
                //if(!PCSGet("MULTI_PHASE_FLOW") || !PCSGet("RICHARDS_FLOW"))
                //{
                //    data_Youngs = new Matrix(9);
                //    //! 0: \f$ \kappa_i0     \f$
                //    //! 1: \f$ \alpha_i     \f$
                //    //! 2: \f$ \kappa_{s0}  \f$
                //    //! 3: \f$ \alpha_{sp}  \f$
                //    //! 4: \f$ \alpha_{ss}  \f$
                //    //! 5: \f$ p_ref        \f$
                //    //! 6: \f$ buffer       \f$
                //    if (SwellingPressureType == 3)
                //        in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >>
                //        (*data_Youngs)(2)
                //        >> (*data_Youngs)(3) >> (*data_Youngs)(4) >>
                //        (*data_Youngs)(5);
                //    else if (SwellingPressureType == 4)
                //        in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >>
                //        (*data_Youngs)(2);
                //    in_sd.clear();
                //}
                //else
                //{
                //    std::cout <<
                //    "No multi-phase flow coupled. The thermal elatic model can only be used in H2 coupled proccess."
                //              << std::endl;
                //    std::cout << "Quit the simulation now!" << std::endl;
                //    abort();
                //}
            }
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$DENSITY") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> Density_mode;
            if(Density_mode == 0) // rho = f(x)
            {
                in_sd >> Size;
                in_sd.clear();
                data_Density = new Matrix(Size, 2);
                for(i = 0; i < Size; i++)
                {
                    in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                    in_sd >> (*data_Density)(i,0) >> (*data_Density)(i,1);
                    in_sd.clear();
                }
            }
            else if(Density_mode == 1) // rho = const
            {
                data_Density = new Matrix(1);
                in_sd >> (*data_Density)(0);
                in_sd.clear();
            }
        }
        //....................................................................
        if(line_string.find("$THERMAL") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> sub_line >> ThermalExpansion;
            in_sd.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("CAPACITY") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> Capacity_mode;
            switch(Capacity_mode)
            {
            case 0:       //  = f(x)
                in_sd >> Size;
                in_sd.clear();
                data_Capacity = new Matrix(Size, 2);
                for(i = 0; i < Size; i++)
                {
                    in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                    in_sd >> (*data_Capacity)(i,0) >> (*data_Capacity)(i,1);
                    in_sd.clear();
                }
                break;
            case 1:       //  = const
                data_Capacity = new Matrix(1);
                in_sd >> (*data_Capacity)(0);
                in_sd.clear();
                break;
            case 2:       // boiling model for rock. WW
                // 0. Wet capacity
                // 1. Dry capacity
                // 2. Boiling temperature
                // 3. Boiling temperature range
                // 4. Latent of vaporization
                data_Capacity = new Matrix(5);
                for(i = 0; i < 5; i++)
                    in_sd >> (*data_Capacity)(i);
                in_sd.clear();
                break;
            case 3:       // DECOVALEX THM1, Bentonite
                in_sd.clear();
                break;
            }
        }

        //....................................................................
        // subkeyword found
        if(line_string.compare("CONDUCTIVITY") == 0)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> Conductivity_mode;
            switch(Conductivity_mode)
            {
            case 0:       //  = f(T) //21.12.2009 WW
                in_sd >> Size;
                in_sd.clear();
                data_Conductivity = new Matrix(Size, 2);
                for(i = 0; i < Size; i++)
                {
                    in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                    in_sd >>
                    (*data_Conductivity)(i,
                                         0) >> (*data_Conductivity)(i,1);
                    in_sd.clear();
                }
                //WW
                conductivity_pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 1:       //  = const
                data_Conductivity = new Matrix(1);
                in_sd >> (*data_Conductivity)(0);
                in_sd.clear();
                break;
            case 2:       // boiling model for rock. WW
                // 0. Wet conductivity
                // 1. Dry conductivity
                // 2. Boiling temperature
                // 3. Boiling temperature range
                data_Conductivity = new Matrix(4);
                for(i = 0; i < 4; i++)
                    in_sd >> (*data_Conductivity)(i);
                in_sd.clear();
                capacity_pcs_name_vector.push_back("TEMPERATURE1");
                capacity_pcs_name_vector.push_back("SATURATION1");
                break;
            case 3:       // DECOVALEX THM1, Bentonite
                in_sd.clear();
                capacity_pcs_name_vector.push_back("TEMPERATURE1");
                capacity_pcs_name_vector.push_back("SATURATION1");
                break;
            case 4:       //  = f(S) //21.12.2009 WW
                in_sd >> Size;
                in_sd.clear();
                data_Conductivity = new Matrix(Size, 2);
                for(i = 0; i < Size; i++)
                {
                    in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                    in_sd >>
                    (*data_Conductivity)(i,
                                         0) >> (*data_Conductivity)(i,1);
                    in_sd.clear();
                }
                break;
            }
            in_sd.clear();
        }

        //....................................................................
        if(line_string.find("CONDUCTIVITY_TENSOR") != string::npos) //subkeyword found
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> thermal_conductivity_tensor_type_name;
            thermal_conductivity_tensor_dim = 0; //NW
            switch(thermal_conductivity_tensor_type_name[0])
            {
            case 'I': // isotropic
                thermal_conductivity_tensor_type = 0;
                in_sd >> thermal_conductivity_tensor[0];
                thermal_conductivity_tensor[1] = thermal_conductivity_tensor[2] =
                                                         thermal_conductivity_tensor
                                                         [0];
                break;
            case 'O':      // orthotropic
                thermal_conductivity_tensor_type = 1;
                in_sd >> thermal_conductivity_tensor_dim;
                if(thermal_conductivity_tensor_dim == 0)
                    std::cout <<
                    "Error in CSolidProperties::Read: no tensor dimension"
                              << std::endl;
                if(thermal_conductivity_tensor_dim == 2)
                {
                    in_sd >> thermal_conductivity_tensor[0];
                    in_sd >> thermal_conductivity_tensor[1];
                }
                if(thermal_conductivity_tensor_dim == 3)
                {
                    in_sd >> thermal_conductivity_tensor[0];
                    in_sd >> thermal_conductivity_tensor[1];
                    in_sd >> thermal_conductivity_tensor[2];
                }
                break;
            case 'A':      // anisotropic
                thermal_conductivity_tensor_type = 2;
                in_sd >> thermal_conductivity_tensor_dim;
                if(thermal_conductivity_tensor_dim == 0)
                    std::cout <<
                    "Error in CSolidProperties::Read: no tensor dimension"
                              << std::endl;
                if(thermal_conductivity_tensor_dim == 2)
                {
                    in_sd >> thermal_conductivity_tensor[0];
                    in_sd >> thermal_conductivity_tensor[1];
                    in_sd >> thermal_conductivity_tensor[2];
                    in_sd >> thermal_conductivity_tensor[3];
                }
                if(thermal_conductivity_tensor_dim == 3)
                {
                    in_sd >> thermal_conductivity_tensor[0];
                    in_sd >> thermal_conductivity_tensor[1];
                    in_sd >> thermal_conductivity_tensor[2];
                    in_sd >> thermal_conductivity_tensor[3];
                    in_sd >> thermal_conductivity_tensor[4];
                    in_sd >> thermal_conductivity_tensor[5];
                    in_sd >> thermal_conductivity_tensor[6];
                    in_sd >> thermal_conductivity_tensor[7];
                    in_sd >> thermal_conductivity_tensor[8];
                }
                break;
            default:
                cout <<
                "Error in CSolidProperties::Read: no valid thermal conductivity tensor type"
                     << endl;
                break;
            }
            in_sd.clear();
        }

        //....................................................................
        // subkeyword found
        if(line_string.find("$ELASTICITY") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> sub_line >> PoissonRatio;
            in_sd.clear();
        }
        //....................................................................
        //12.2009. WW
        if(line_string.find("$EXCAVATION") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> excavation;
            in_sd.clear();
        }
        //....................................................................
        if(line_string.find("YOUNGS_MODULUS") != string::npos)
        {                         // subkeyword found
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> Youngs_mode;
            int type = Youngs_mode;
            if(type > 9 && type < 14)
                type = 1000;
            switch(type)  // 15.03.2008 WW
            {
            case 0:       //  = f(x)
                in_sd >> Size;
                in_sd.clear();
                data_Youngs = new Matrix(Size, 2);
                for(i = 0; i < Size; i++)
                {
                    in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                    in_sd >> (*data_Youngs)(i,0) >> (*data_Youngs)(i,1);
                    in_sd.clear();
                }
                break;
            case 1:       //  = const
                data_Youngs = new Matrix(1);
                in_sd >> (*data_Youngs)(0);
                in_sd.clear();
                break; // UJG 24.11.2009
            case 2:       //  = const
                // data_Youngs Lubby1 model
                //  0: E_0
                //  1: a (factor)
                //  2: n (exponent)
                data_Youngs = new Matrix(3);
                in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2);
                in_sd.clear();
                break;
            case 1000:    // case 10-13: transverse isotropic linear elasticity (UJG 24.11.2009)
                // data_Youngs transverse isotropic linear elasticity
                //  0: E_i     (Young's modulus of the plane of isotropy)
                //  1: E_a     (Young's modulus w.r.t. the anisotropy direction)
                //  2: nu_{ia} (Poisson's ratio w.r.t. the anisotropy direction)
                //  3: G_a     (shear modulus w.r.t. the anisotropy direction)
                //  4: n_x     (x-coefficient of the local axis of anisotropy (2D case: -\sin\phi))
                //  5: n_y     (y-coefficient of the local axis of anisotropy (2D case: \cos\phi))
                //  6: n_z     (z-coefficient of the local axis of anisotropy (2D case: 0))
                data_Youngs = new Matrix(7);
                in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >>
                (*data_Youngs)(2) >> (*data_Youngs)(3)
                >> (*data_Youngs)(4) >> (*data_Youngs)(5) >> (*data_Youngs)(6);
                in_sd.clear();
                break;
            }
        }
        //....................................................................
        if(line_string.find("$CREEP") != string::npos)
        {
            if(line_string.find("NORTON") != string::npos)
            {
                Creep_mode = 1;
                /*! \subsection Norton creep model */
                /*! \f$\dot\epsilon_s=A \left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
                // data_Creep:
                //  0: A,   coefficient
                //  1: n,   exponential
                data_Creep = new Matrix(2);
                in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                in_sd >> (*data_Creep)(0);
                in_sd >> (*data_Creep)(1);
                in_sd.clear();
            }
            if(line_string.find("BGRA") != string::npos)
            {
                Creep_mode = 2;
                /*! \subsection Temperature dependent creep model by BGR */
                /*! \f$\dot\epsilon_s=A\exp^{-Q/RT}\left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
                // data_Creep:
                //  0: A,   coefficient
                //  1: n,   exponential
                //  2: Q,   activation energy
                data_Creep = new Matrix(3);
                in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                in_sd >> (*data_Creep)(0);
                in_sd >> (*data_Creep)(1);
                in_sd >> (*data_Creep)(2);
                in_sd.clear();
            }
            //....................................................................
            if(line_string.find("LUBBY2") != string::npos)
            {
                Creep_mode = 1000;
                // data_Creep:
                //  0: eta_m
                //  1: m
                //  2: l
                //  3: eta_k
                //  4: k1
                //  5: k2
                //  6: G_k
                data_Creep = new Matrix(7,2);
                in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                for(i = 0; i < 7; i++)
                    in_sd >> (*data_Creep)(i,0);
                in_sd.clear();
            }
        }
        //....................................................................
        if(line_string.find("$BIOT_CONSTANT") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> biot_const;
            in_sd.clear();
        }
        //....................................................................
        if(line_string.find("BISHOP_COEFFICIENT") != string::npos) //WX
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> bishop_model;
            switch(bishop_model)
            {
            case 1: //constant
                in_sd >> bishop_model_value;
                break;
            case 2: //pow(Se, parameter)
                in_sd >> bishop_model_value;
                break;
            default:
                break;
            }
            in_sd.clear();
        }
        //....................................................................
        if(line_string.find("$STRESS_INTEGRATION_TOLERANCE") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> f_tol >> s_tol;
            in_sd.clear();
        }
        if(line_string.find("$STRESS_UNIT") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> sub_line;
            if(sub_line.compare("MegaPascal") == 0)
                grav_const = 9.81e-6;
            else if(sub_line.find("KiloPascal") == 0)
                grav_const = 9.81e-3;
            in_sd.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$PLASTICITY") != string::npos)
        {
            in_sd.str(readNonBlankLineFromInputStream(*msp_file));
            in_sd >> sub_line;
            in_sd.clear();
            if(sub_line.find("DRUCKER-PRAGER") != string::npos)
            {
                devS = new double[6];
                Plasticity_type = 1;
                // No return mapping
                if(sub_line.find("NORETURNMAPPING") != string::npos)
                {
                    Plasticity_type = 10;
                    if (sub_line.find("TENSIONCUTOFF") != string::npos) //WX: 08.2010
                    {
                        Plasticity_type = 11;
                        dFtds = new double[6]; //WX: 08.2010 Tensile yield function
                        dGtds = new double[6];
                        ConstitutiveMatrix = new Matrix(6,6); //WX: 08.2010
                    }
                    dFds = new double[6];
                    dGds = new double[6];
                    D_dFds = new double[6];
                    D_dGds = new double[6];
                }
                Size = 5;
                if(Plasticity_type == 11)
                    Size = 6;
                /*
                   Material parameters for Cam-Clay model
                   i : parameter
                   0 : The initial yield stress
                   1 : Plastic hardening parameter
                   2 : Internal frection angle
                   3 : Dilatancy angle
                   4 : Localized softening modulus
                   5 : Tensile strength //WX
                 */
            }
            else if(sub_line.find("SINGLE_YIELD_SURFACE") != string::npos)
            {
                Plasticity_type = 2;
                Size = 23;
                /*
                   Material parameters for Single yield surface model
                   i: parameter
                   0: alpha0
                   1: beta0
                   2: delta0
                   3: epsilon0
                   4: kappa0
                   5: gamma0
                   6: m0

                   7: alpha1
                   8: beta1
                   9: delta1
                   10: epsilon1
                   11: kappa1
                   12: gamma1
                   13: m1

                   14: Psi1
                   15: Psi2

                   16: Ch
                   17: Cd
                   18: br
                   19: mr

                   20: Initial stress_xx
                   21: Initial stress_yy
                   22: Initial stress_zz
                 */
            }
            else if(sub_line.find("CAM-CLAY") != string::npos)
            {
                Plasticity_type = 3;
                Size = 10;
                /*
                   Material parameters for Cam-Clay model
                   i: parameter
                   0 : M: slope of the critical line
                   1 : Lamda, the virgin compression index
                   2 : Kappa, swelling index
                   3 : p0, preconsolidation pressure
                   4 : e0, initial void ratio
                   5 : OCR
                   6 : Initial stress_xx
                   7 : Initial stress_yy
                   8 : Initial stress_zz
                   9 : Mimimum p: ( stress_xx + stress_yy + stress_zz )/3
                 */
            }
            else if(sub_line.find("MOHR-COULOMB") != string::npos) //WX
            {
                devS = new double[6];
                ConstitutiveMatrix = new Matrix(6,6);
                Plasticity_type = 4;
                Size = 6;
                /*
                   i    parameter
                   0    cohesion
                   1    friction angle
                   2    dilatance angle
                   3    tension strength
                   4   hardening curve for friction angle
                   5   hardening curve for cohesion
                 */
            }
            else if(sub_line.find("HOEK-BROWN") != string::npos) //WX
            {
                devS = new double[6];
                ConstitutiveMatrix = new Matrix(6,6);
                Plasticity_type = 5;
                Size = 4;
                /*
                   i   parameter
                   0   a
                   1   s
                   2   mb
                   3   sigci
                 */
            }
            data_Plasticity = new Matrix(Size);
            for(i = 0; i < Size; i++)
            {
                in_sd.str(readNonBlankLineFromInputStream(*msp_file));
                in_sd >> (*data_Plasticity)(i);
                in_sd.clear();
            }
        }
        in_sd.clear();
    }
    return position;
}

//==========================================================================

/**************************************************************************
   FEMLib-Method:
   Task: Constructor and destructor
   Programing:
   08/2004 WW Implementation
   04/2005 WW Curve dependency
**************************************************************************/
CSolidProperties::CSolidProperties()
    : data_Youngs(NULL), data_Density(NULL), data_Capacity(NULL),
      data_Conductivity(NULL), data_Plasticity(NULL), data_Creep(NULL)
{
    PoissonRatio = 0.2;
    ThermalExpansion = 0.0;
    biot_const = 1.0;
    // Default, all data are constant
    Density_mode = -1;
    Youngs_mode = -1;
    Capacity_mode = -1;
    Conductivity_mode = -1;
    Creep_mode = -1;
    grav_const = 9.81;                    //WW
    excavation = -1;                      //12.2009. WW
    excavated = false;                    //To be .....  12.2009. WW

    SwellingPressureType = -1;
    Max_SwellingPressure = 0.0;
    // Default, elasticity
    Plasticity_type = -1;

    E = Lambda = G = K = 0.0;
    devS = NULL;
    axisymmetry = false;
    dl2 = 0.0;
    // SYS
    d2G_dSdS = NULL;
    d2G_dSdM = NULL;
    LocalJacobi = NULL;
    inv_Jac = NULL;
    sumA_Matrix = NULL;
    rhs_l = NULL;
    x_l = NULL;
    Li = NULL;
    // Drucker-Prager
    dFds = NULL;
    dGds = NULL;
    D_dFds = NULL;
    D_dGds = NULL;
    dFtds = NULL; //WX
    dGtds = NULL; //WX
    ConstitutiveMatrix = NULL; //WX
    // Curve variable type
    // 0: Time
    // 1: ...
    CurveVariableType_Conductivity = -1;
    mode = 0;                             // Gauss point values //OK
    //
    s_tol = 1e-9;
    f_tol = 1e-9;

    Crotm = NULL;                         // rotation matrix for matrices: UJG 25.11.2009
    D_tran = NULL;                        //UJG/WW

    // Thermal conductivity tensor (default: iso)
    thermal_conductivity_tensor_type = 0;
    thermal_conductivity_tensor_dim = 1;
    thermal_conductivity_tensor[0] = 1.0;

    bishop_model = -1; //15.08.2011. WW

    Al = .0;
    BetaN = .0;
    bishop_model_value = .0;
    csn = .0;
    Hard_Loc = .0;
    Hard = .0;
    HoekB_a = .0;
    HoekB_cohe = .0;
    HoekB_mb = .0;
    HoekB_s = .0;
    HoekB_sigci = .0;
    HoekB_tens = .0;
    Nphi = .0;
    Ntheta = .0;
    tension = .0;
    Xi = .0;
    Y0 = .0;

}
CSolidProperties::~CSolidProperties()
{
    if(data_Density)
        delete data_Density;
    if(data_Density)
        delete data_Youngs;
    if(data_Plasticity)
        delete data_Plasticity;
    if(data_Capacity)
        delete data_Capacity;
    if(data_Conductivity)
        delete data_Conductivity;
    if(data_Creep)
        delete data_Creep;
    if(devS)
        delete [] devS;

    if(Crotm)
        delete Crotm;             // rotation matrix for matrices: UJG 25.11.2009
    if(D_tran)
        delete D_tran;            // rotation matrix for matrices: UJG 25.11.2009
    data_Density = NULL;
    data_Youngs = NULL;
    data_Plasticity = NULL;
    data_Capacity = NULL;
    data_Conductivity = NULL;
    data_Creep = NULL;
    devS = NULL;
    Crotm = NULL;
    D_tran = NULL;

    if(d2G_dSdS)
        delete d2G_dSdS;
    if(d2G_dSdM)
        delete d2G_dSdM;
    if(LocalJacobi)
        delete LocalJacobi;       // To store local Jacobi matrix
    if(inv_Jac)
        delete inv_Jac;           // To store the inverse of the  Jacobi matrix
    if(sumA_Matrix)
        delete sumA_Matrix;
    if(rhs_l)
        delete [] rhs_l;          // To store local unknowns of 15
    if(x_l)
        delete [] x_l;            // To store local unknowns of 15
    if(Li)
        delete [] Li;

    if(dFds)
        delete [] dFds;
    if(dGds)
        delete [] dGds;
    if(D_dFds)
        delete [] D_dFds;
    if(D_dGds)
        delete [] D_dGds;
    if(dFtds)
        delete [] dFtds;  //WX:
    if(dGtds)
        delete [] dGtds;  //WX:
    if(ConstitutiveMatrix)
        delete ConstitutiveMatrix;     //WX:
    dFds = NULL;
    dGds = NULL;
    D_dFds = NULL;
    D_dGds = NULL;
    dFtds = NULL; //WX:
    dGtds = NULL; //WX:
    ConstitutiveMatrix = NULL; //WX:
    d2G_dSdS = NULL;
    d2G_dSdM = NULL;
    LocalJacobi = NULL;
    inv_Jac = NULL;
    sumA_Matrix = NULL;
    rhs_l = NULL;
    x_l = NULL;
    Li = NULL;
}
}
//----------------------------------------------------------------------------

