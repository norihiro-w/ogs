/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_mmp_new.cpp
 *
 * Created on 2004-01-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib-Object: MAT-MP
   Task: MediumProperties
   Programing:
   01/2004 OK Implementation
**************************************************************************/

#include "rf_mmp_new.h"

//#include "makros.h"
// C++ STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <cstdlib>
#include <string>
#include <cstring>

#include "makros.h" // JT
#include "Ogs5FileTools.h"


using namespace std;

namespace ogs5
{

void ScreenMessage(const char* message)
{
#ifdef USE_MPI
    if(myrank > 0)
        return;
#endif
    printf("%s", message);
}

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: constructor
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties::CMediumProperties() :
    geo_dimension(0)
{
    name = "DEFAULT";
    mode = 0;
    selected = false;
    // GEO
    geo_area = 1.0;                       //OK
    porosity_model = 1;
    porosity_model_values[0] = 0.1;
    tortuosity_model = 1;
    tortuosity_model_values[0] = 1.;
    // flow
    storage_model = 1;
    storage_model_values[0] = 0.;
    permeability_model = 1;
    permeability_tensor_type = 0;
    tortuosity_tensor_type = 0;
    permeability_tensor[0] = 1.e-13;
    residual_saturation[0] = 0.0;            // sgr: residual saturation, this phase
    maximum_saturation[0] = 1.0;            // sgm: maximum saturation, this phase
    saturation_exponent[0] = 1.0;            // (set exponent = 1 results in a linear k_rel function)
    conductivity_model = -1;
    flowlinearity_model = 0;
    capillary_pressure_model = -1;
    capillary_pressure_values[4] = 1.0/DBL_EPSILON; // JT: max Pc
    entry_pressure_conversion = false;
    permeability_saturation_model[0] = -1;
    minimum_relative_permeability = 1.0e-9;     // JT: the default value
    unconfined_flow_group = -1;
    permeability_stress_mode = -1;        //WW
    c_coefficient = NULL;                 //WW
    // surface flow
    friction_coefficient = -1;
    //  friction_model = -1;
    // mass transport
    // heat transport
    heat_dispersion_model = -1;           //WW
    heat_dispersion_longitudinal = 0.;
    heat_dispersion_transverse = 0.;
    lgpn = 0.0;
    mass_dispersion_transverse = 0.0;
    mass_dispersion_longitudinal = 0.0;
    heat_diffusion_model = -1;            //WW
    geo_area = 1.0;
//    geo_type_name = "DOMAIN";             //OK
    vol_mat = 0.0;
    vol_bio = 0.0;
    vol_mat_model = 0;
    vol_bio_model = 0;
    foc = 0.0;

    permeability_pressure_model = -1; //01.09.2011. WW
    permeability_strain_model = -1; //01.09.2011. WW

    argument = .0;
    this->channel = false;
    this->conductivity = .0;
    this->density = .0;
    this->evaporation = 0;
    fct_number = 0;
    friction_exp_depth = .0;
    friction_exp_slope = .0;
    heat_capacity = .0;
    heatflux = .0;
    KC_permeability_initial = .0;
    KC_porosity_initial = .0;
    mass_dispersion_model = 0;
    number = 0;
    overland_width = .0;
    permeability_porosity_model = 0;
    permeability_pressure_rel = .0;
    permeability = .0;
    porosity_curve = 0;
    specific_storage = .0;
    storativity = .0;
    tortuosity = .0;
    vaporfraction = .0;
    transfer_coefficient = .0;
    rill_epsilon = .0;
    rill_height = .0;

    is_fracture = false; //NW
    heat_capacity_model = 1; // default
    heat_conductivity_model = 1;
}

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: destructor
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties::~CMediumProperties(void)
{
    if(c_coefficient)
        delete[] c_coefficient;   //WW
    geo_name_vector.clear();
}

////////////////////////////////////////////////////////////////////////////
// IO functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method: MMPRead
   Task: master read functionn
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
bool MMPRead(const std::string& base_file_name, std::vector<CMediumProperties*> &mmp_vector)
{
    //----------------------------------------------------------------------
    //OK  MMPDelete();
    //----------------------------------------------------------------------
    std::cout << "MMPRead ... " << std::flush;
    CMediumProperties* m_mat_mp = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // file handling
    std::string mp_file_name;
    mp_file_name = base_file_name + MMP_FILE_EXTENSION;
    std::ifstream mp_file (mp_file_name.data(),std::ios::in);
    if (!mp_file.good())
    {
        std::cout << "! Error in MMPRead: No material data !" << std::endl;
        return false;
    }
    mp_file.seekg(0L,std::ios::beg);
    //========================================================================
    // keyword loop
    while (!mp_file.eof())
    {
        mp_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != string::npos)
        {
            std::cout << "done, read " << mmp_vector.size() << " medium properties" <<
            std::endl;
            return true;
        }
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#MEDIUM_PROPERTIES") != string::npos)
        {
            m_mat_mp = new CMediumProperties();
            position = m_mat_mp->Read(&mp_file);
            //OK41
            m_mat_mp->number = (int)mmp_vector.size();
            mmp_vector.push_back(m_mat_mp);
            mp_file.seekg(position,std::ios::beg);
        }                         // keyword found
    }                                     // eof
    return true;
    // Tests
}

/**************************************************************************
   FEMLib-Method: CMediumProperties::Read
   Task: read functionn
   Programing:
   02/2004 OK Template
   08/2004 CMCD Implementation
   10/2004 MX/OK Porosity model 3, swelling
   11/2004 CMCD String streaming
   07/2005 MB porosity_file, permeability_file, GEO_TYPE layer
   10/2005 OK GEO_TYPE geo_name_vector
   01/2006 YD PCS_TYPE
   05/2007 PCH Tortuosity tensor
   05/2007 WW Stress permeability coorector. Two models.
   last modification:
**************************************************************************/
// Order of Key Words
/*
         0. $NAME
            (i)    _BORDEN
         1. $GEOTYPE
            (i)        _CLAY
            (ii)    _SILT
            (iii)    _SAND
            (iv)    _GRAVEL
            (v)        _CRYSTALINE
         2. $GEOMETRY
            (i)        _DIMENSION
   (ii)    _AREA
   3. $POROSITY
   4. $TORTUOSITY
   5. $MOBILE_IMOBILE_MODEL
   6. $LITHOLOGY_GRAIN_CLASS
   7. $FLOWLINEARITY
   8. $SORPTION_MODEL
   9. $STORAGE
   11.$PERMEABILITY
   12.$PERMEABILITY_FUNCTION_
   (1)        DEFORMATION
   (2)  PRESSURE
   (3)        SATURATION
   (4)        STRESS
   (5)        VELOCITY
   (6)        POROSITY
   13.$CAPILLARY_PRESSURE
   14.$MASSDISPERSION
   (i)        _LONGITUDINAL
   (ii)    _TRANSVERSE
   15.$HEATDISPERSION
   (i)        _LONGITUDINAL
   (ii)    _TRANSVERSE
   19.$ELECTRIC_CONDUCTIVITY
   20.$UNCONFINED_FLOW_GROUP
   21.$FLUID_EXCHANGE_WITH_OTHER_CONTINUA
 */
std::ios::pos_type CMediumProperties::Read(std::ifstream* mmp_file)
{
    int i, j, k = 0;
    std::string line_string;
    std::stringstream in;
    std::ios::pos_type position;
    std::string dollar("$");
    std::string hash("#");
    //WW bool new_subkeyword = false;
    bool new_keyword = false;
    std::string m_string;
    // WW
    std::stringstream buff;
    std::vector<string> tokens;
    char* pch;
    char seps[] = "+\n";
    char seps1[] = "*";
    double f_buff;
//    size_t indexChWin, indexChLinux;      //JT, DEC 2009
    std::string funfname;                 //JT, DEC 2009

    while (!new_keyword)
    {
        //WW new_subkeyword = false;
        position = mmp_file->tellg();
        line_string = readNonBlankLineFromInputStream(*mmp_file);
        if(line_string.size() < 1)
            break;
        if(line_string.find(hash) != std::string::npos)
        {
            new_keyword = true;
            break;
        }
        //--------------------------------------------------------------------
        //PCS                         //YD
        if(line_string.find("$PCS_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> pcs_type_name;
            in.clear();
            continue;
        }

        //NAME
        //subkeyword found
        if(line_string.find("$NAME") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> name;   //sub_line
            in.clear();
            continue;
        }
        //--------------------------------------------------------------------
        //GEO_TYPE

        //subkeyword found
        if(line_string.find("$GEO_TYPE") != std::string::npos)
        {
            while(!(m_string.find("$") != std::string::npos) &&
                  (!(m_string.find("#") != std::string::npos)))
            {
                position = mmp_file->tellg();
                in.str(readNonBlankLineFromInputStream(*mmp_file));
                in >> m_string >> geo_name;
                in.clear();
            }
            mmp_file->seekg(position,std::ios::beg);
            continue;
        }
        //....................................................................
        // ToDo to GeoLib
        //2i..GEOMETRY_DIMENSION
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$GEOMETRY_DIMENSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> geo_dimension;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // ToDo to GeoLib
        //2ii..GEOMETRY_AREA
        //------------------------------------------------------------------------
//        indexChWin = indexChLinux = 0; //JT, DEC 2009
                                       //subkeyword found
        if(line_string.find("$GEOMETRY_AREA") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> line_string;
            if(line_string.find("FILE") != string::npos)
            {
                in >> geo_area_file;
            }
            else
                geo_area = strtod(line_string.data(),NULL);
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //3..POROSITY
        //------------------------------------------------------------------------
        //CB
        //subkeyword found
        if((line_string.find("$POROSITY") != string::npos) &&
           (!(line_string.find("_DISTRIBUTION") != std::string::npos)))
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> porosity_model;
            switch(porosity_model)
            {
            case 0:       // n=f(x)
                in >> porosity_curve;
                break;
            case 1:       // n=const
                in >> porosity_model_values[0];
                break;
            case 2:       // f(normal effective stress for fracture systems)
                in >> porosity_model_values[0];
                in >> porosity_model_values[1];
                in >> porosity_model_values[2];
                in >> porosity_model_values[3];
                porosity_pcs_name_vector.push_back("PRESSURE1");
                break;
            case 3:       // Chemical swelling model
                in >> porosity_model_values[0]; // Initial porosity
                in >> porosity_model_values[1]; // Specific surface[m^2/g]
                in >> porosity_model_values[2]; // Expansive min. fragtion
                in >> porosity_model_values[3]; // m
                in >> porosity_model_values[4]; // I
                in >> porosity_model_values[5]; // S^l_0
                in >> porosity_model_values[6]; // beta
                porosity_pcs_name_vector.push_back("SATURATION2");
                porosity_pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 4:       // Chemical swelling model (constrained swelling, constant I)
                in >> porosity_model_values[0]; // Initial porosity
                in >> porosity_model_values[1]; // Specific surface[m^2/g]
                in >> porosity_model_values[2]; // Expansive min. fragtion
                in >> porosity_model_values[3]; // m
                in >> porosity_model_values[4]; // I
                in >> porosity_model_values[5]; // S^l_0
                in >> porosity_model_values[6]; // beta
                in >> porosity_model_values[7]; // n_min
                                                //for richard flow only
                porosity_pcs_name_vector.push_back("SATURATION1");
                porosity_pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 5:       // Chemical swelling model (free swelling, constant I)
                in >> porosity_model_values[0]; // Initial porosity
                in >> porosity_model_values[1]; // Specific surface[m^2/g]
                in >> porosity_model_values[2]; // Expansive min. fragtion
                in >> porosity_model_values[3]; // m
                in >> porosity_model_values[4]; // I
                in >> porosity_model_values[5]; // S^l_0
                in >> porosity_model_values[6]; // beta
                porosity_pcs_name_vector.push_back("SATURATION2");
                porosity_pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 6:       // Chemical swelling model (constrained swelling)
                in >> porosity_model_values[0]; // Initial porosity
                in >> porosity_model_values[1]; // Specific surface[m^2/g]
                in >> porosity_model_values[2]; // Expansive min. fragtion
                in >> porosity_model_values[3]; // m
                in >> porosity_model_values[4]; // I
                in >> porosity_model_values[5]; // S^l_0
                in >> porosity_model_values[6]; // beta
                in >> porosity_model_values[7]; // n_min
                porosity_pcs_name_vector.push_back("SATURATION2");
                porosity_pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 7:       // n=f(stress_mean) WW
                in >> porosity_curve;
                break;
            case 10:      // Chemical swelling model (constrained swelling, constant I)
            {
                int m;
                in >> porosity_model_values[0]; // Initial porosity
                in >> m; // m
                if (m > 15)
                    std::cout
                    <<
                    "Maximal number of solid phases is now limited to be 15!!!"
                    << std::endl;
                for (int i = 0; i < m + 1; i++)
                    // molar volume [l/mol]
                    in >> porosity_model_values[i + 1];
                break;
            }
            case 11:      //MB: read from file ToDo
                // in >> porosity_file; // CB
                in >> porosity_model_values[0]; // CB some dummy default value is read
                // CB $POROSITY_DISTRIBUTION should be given as keyword in *.mmp file,
                //     porosities then are to be read in from file by fct.
                //     CMediumProperties::SetDistributedELEProperties
                break;
            case 12:
                in >> porosity_model_values[0]; //WX 03.2011, dependent on strain
                break;
#ifdef GEM_REACT
            case 15:
                in >> porosity_model_values[0]; // set a default value for GEMS calculation
                                                // save this seperately;
                KC_porosity_initial = porosity_model_values[0];

                // KG44: TODO!!!!!!!!!!!!! check the above  ***************

                break;
#endif
#ifdef BRNS
            case 16:
                in >> porosity_model_values[0]; // set a default value for BRNS calculation
                break;
#endif
            default:
                std::cerr << "Error in MMPRead: no valid porosity model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //subkeyword found
        if(line_string.find("$VOL_MAT") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            // in >> idummy >> this->vol_mat; CB
            in >> vol_mat_model >> this->vol_mat;
            switch(vol_mat_model)
            {
            case 1:       // do nothing
                break;
            case 2:       // do nothing
                break;
            default:
                std::cout << "Error in MMPRead: no valid vol_mat_model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //subkeyword found
        if(line_string.find("$VOL_BIO") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            // in >> idummy >> this->vol_bio; CB
            in >> vol_bio_model >> this->vol_bio;
            switch(vol_bio_model)
            {
            case 1:       // do nothing
                break;
            case 2:       // do nothing
                break;
            default:
                std::cout << "Error in MMPRead: no valid vol_bio_model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //4..TORTUOSITY
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$TORTUOSITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> tortuosity_tensor_type_name;
            switch(tortuosity_tensor_type_name[0])
            {
            case '0':     // n=f(x) <- Case zero
                break;
            case '1':     // n=const
                tortuosity_model = 1;
                in >> tortuosity_model_values[0];
                break;
            case 'I':     // isotropic
                tortuosity_tensor_type = 0;
                tortuosity_model = 1;
                in >> tortuosity_model_values[0];
                //CMCD to pick up 2D and 3D Isotropic case.
                tortuosity_model_values[1] = tortuosity_model_values[2] =
                                                     tortuosity_model_values[0];
                break;
            case 'O':     // orthotropic  <- Case Alphabet O
                tortuosity_model = 1;
                tortuosity_tensor_type = 1;
                if(geo_dimension == 0)
                    std::cout <<
                    "Error in CMediumProperties::Read: no geometric dimension"
                              << std::endl;
                if(geo_dimension == 2)
                {
                    in >> tortuosity_model_values[0];
                    in >> tortuosity_model_values[1];
                }
                if(geo_dimension == 3)
                {
                    in >> tortuosity_model_values[0];
                    in >> tortuosity_model_values[1];
                    in >> tortuosity_model_values[2];
                }
                break;
            case 'A':     // anisotropic
                tortuosity_model = 1;
                tortuosity_tensor_type = 2;
                if(geo_dimension == 0)
                    std::cout <<
                    "Error in CMediumProperties::Read: no geometric dimension"
                              << std::endl;
                if(geo_dimension == 2)
                {
                    in >> tortuosity_model_values[0];
                    in >> tortuosity_model_values[1];
                    in >> tortuosity_model_values[2];
                    in >> tortuosity_model_values[3];
                }
                if(geo_dimension == 3)
                {
                    in >> tortuosity_model_values[0];
                    in >> tortuosity_model_values[1];
                    in >> tortuosity_model_values[2];
                    in >> tortuosity_model_values[3];
                    in >> tortuosity_model_values[4];
                    in >> tortuosity_model_values[5];
                    in >> tortuosity_model_values[6];
                    in >> tortuosity_model_values[7];
                    in >> tortuosity_model_values[8];
                }
                break;
            case 'F':     //SB: read from file
                tortuosity_model = 2; //OK
                in >> permeability_file;
                break;
            default:
                std::cout << "Error in MMPRead: no valid tortuosity tensor type" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        //5..MOBILE_IMOBILE_MODEL
        //------------------------------------------------------------------------
        //To do as necessary

        //------------------------------------------------------------------------
        //6..LITHOLOGY_GRAIN_CLASS
        //------------------------------------------------------------------------
        //To do as necessary

        //------------------------------------------------------------------------
        //7..FLOWLINEARITY
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$FLOWLINEARITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> flowlinearity_model;
            switch(flowlinearity_model)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // Read in alpha
                //Alpha
                in >> flowlinearity_model_values[0];
                break;
            case 2:       // For equivalent flow in trianglular elements
                //Alpha
                in >> flowlinearity_model_values[0];
                //Number of Fractures in Equivalent Medium
                in >> flowlinearity_model_values[1];
                //Reynolds Number above which non linear flow occurs.
                in >> flowlinearity_model_values[2];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            default:
                std::cout << "Error in MMPRead: no valid flow linearity model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //8..SORPTION_MODEL
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$ORGANIC_CARBON") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> foc;
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        //9..STORAGE
        //------------------------------------------------------------------------

        //subkeyword found
        if(line_string.find("$STORAGE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> storage_model;
            switch(storage_model)
            {
            case 0:       // S=f(x)
                in >> storage_model_values[0]; //Function of pressure defined by curve
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 1:       // S=const
                in >> storage_model_values[0]; //Constant value in Pa
                break;
            case 2:
                in >> storage_model_values[0]; //S0
                in >> storage_model_values[1]; //increase
                in >> storage_model_values[2]; //sigma(z0)
                in >> storage_model_values[3]; //d_sigma/d_z
                break;
            case 3:
                in >> storage_model_values[0]; //curve number (as real number)
                in >> storage_model_values[1]; //sigma(z0)
                in >> storage_model_values[2]; //_sigma/d_z
                break;
            case 4:
                in >> storage_model_values[0]; //curve number (as real number)
                in >> storage_model_values[1]; //time collation
                in >> storage_model_values[2]; //solid density
                in >> storage_model_values[3]; //curve fitting factor, default 1
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 5:       //Storativity is a function of normal effective stress defined by curve, set up for KTB.
                in >> storage_model_values[0]; //curve number
                in >> storage_model_values[1]; //time collation
                in >> storage_model_values[2]; //Default storage value for material groups > 0
                in >> storage_model_values[3]; //Angular difference between Y direction and Sigma 1
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 6:       //Storativity is a function of normal effective stress defined by curve and distance from borehole, set up for KTB.
                in >> storage_model_values[0]; //curve number
                in >> storage_model_values[1]; //time collation
                in >> storage_model_values[2]; //Default storage value for material groups > 0
                in >> storage_model_values[3]; //Angular difference between Y direction and Sigma 1
                in >> storage_model_values[4]; //Borehole (x) coordinate
                in >> storage_model_values[5]; //Borehole (y) coordinate
                in >> storage_model_values[6]; //Borehole (z) coordinate
                in >> storage_model_values[7]; //Maximum thickness of shear zone
                in >> storage_model_values[8]; //Fracture density
                pcs_name_vector.push_back("PRESSURE1");
                break;
            default:
                cout << "Error in MMPRead: no valid storativity model" << endl;
                break;
            case 7:       //RW/WW
                in >> storage_model_values[0]; //Biot's alpha
                in >> storage_model_values[1]; //Skempton's B coefficient
                in >> storage_model_values[2]; //macroscopic drained bulk modulus
                double val_l = storage_model_values[0] *
                               (1. - storage_model_values[0] *
                                storage_model_values[1])
                               / storage_model_values[1] / storage_model_values[2];
                storage_model_values[1] = val_l;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //10..CONDUCTIVITY_MODEL
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$CONDUCTIVITY_MODEL") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> conductivity_model;
            switch(conductivity_model)
            {
            case 0:       // K=f(x)
                break;
            case 1:       // K=const
                in >> conductivity;
                break;
            case 2:       // Manning
                break;
            case 3:       // Chezy
                break;
            default:
                std::cout << "Error in MMPRead: no valid conductivity model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }

        //subkeyword found
        if(line_string.find("$UNCONFINED") != string::npos)
        {
            unconfined_flow_group = 1;
            continue;
        }

        //------------------------------------------------------------------------
        //11..PERMEABILITY_TENSOR
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$PERMEABILITY_TENSOR") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_tensor_type_name;
            switch(permeability_tensor_type_name[0])
            {
            case 'I':     // isotropic
                permeability_tensor_type = 0;
                permeability_model = 1;
                in >> permeability_tensor[0];
                //CMCD to pick up 2D and 3D Isotropic case.
                permeability_tensor[1] = permeability_tensor[2] =
                                                 permeability_tensor[0];
                break;
            case 'O':     // orthotropic
                permeability_tensor_type = 1;
                if(geo_dimension == 0)
                    std::cout <<
                    "Error in CMediumProperties::Read: no geometric dimension"
                              << std::endl;
                if(geo_dimension == 2)
                {
                    in >> permeability_tensor[0];
                    in >> permeability_tensor[1];
                }
                if(geo_dimension == 3)
                {
                    in >> permeability_tensor[0];
                    in >> permeability_tensor[1];
                    in >> permeability_tensor[2];
                }
                break;
            case 'A':     // anisotropic
                permeability_tensor_type = 2;
                if(geo_dimension == 0)
                    std::cout <<
                    "Error in CMediumProperties::Read: no geometric dimension"
                              << std::endl;
                if(geo_dimension == 2)
                {
                    in >> permeability_tensor[0];
                    in >> permeability_tensor[1];
                    in >> permeability_tensor[2];
                    in >> permeability_tensor[3];
                }
                if(geo_dimension == 3)
                {
                    in >> permeability_tensor[0];
                    in >> permeability_tensor[1];
                    in >> permeability_tensor[2];
                    in >> permeability_tensor[3];
                    in >> permeability_tensor[4];
                    in >> permeability_tensor[5];
                    in >> permeability_tensor[6];
                    in >> permeability_tensor[7];
                    in >> permeability_tensor[8];
                }
                break;
            case 'F':     //SB: read from file
                permeability_model = 2; //OK
                in >> permeability_file;
                break;
            default:
                std::cout <<
                "Error in MMPRead: no valid permeability tensor type" << std::endl;
                break;
            }
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        //12. $PERMEABILITY_FUNCTION
        //                (i)        _DEFORMATION
        //                (ii)    _PRESSURE
        //                (iii)    _SATURATION
        //                (iv)    _STRESS
        //                (v)        _VELOCITY
        //                (vi)    _POROSITY
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //12.1 PERMEABILITY_FUNCTION_DEFORMATION
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$PERMEABILITY_FUNCTION_DEFORMATION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_model;
            switch(permeability_model)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // k=const
                in >> permeability;
                break;
            default:
                std::cout << "Error in MMPRead: no valid permeability model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //WX: 05.2010
        if(line_string.find("$PERMEABILITY_FUNCTION_STRAIN") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_strain_model;
            switch(permeability_strain_model)
            {
            case 0: //strain_volume
                break;
            case 1: //strain_volumeeff plas strain
                in >> permeability_strain_model_value[0];
                break;
            case 2: //eff plas strainif eff plas strain>0, f(strainp). else stain volume
                in >> permeability_strain_model_value[0];
                break;
            case 3: //if eff plas strain>0, f(strainp). else stain volume
                in >> permeability_strain_model_value[0]; //for strain volume
                in >> permeability_strain_model_value[1]; //for eff plas strain
                break;
            case 4: //strain volume + eff plas strain
                in >> permeability_strain_model_value[0];
                in >> permeability_strain_model_value[1];
                break;
            default:
                cout << "Error in MMPRead: no valid permeability strain model" <<
                endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //12.2 PERMEABILITY_FUNCTION_PRESSURE
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$PERMEABILITY_FUNCTION_PRESSURE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_pressure_model;
            switch(permeability_pressure_model)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // k=const
                in >> permeability_pressure_model_values[0];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 2:       // Permebility is a function of effective stress
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                in >> permeability_pressure_model_values[3];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 3:       // Permeability is a function of non linear flow
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                in >> permeability_pressure_model_values[3];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 4:       // Function of effective stress from a curve
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 5:       // Function of overburden converted to effective stress and related to a curve.
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 6:       // Normal effective stress calculated according to fracture orientation, related over a curve : KTB site
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                in >> permeability_pressure_model_values[3];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 7:       // Normal effective stress calculated according to fracture orientation, related over a curve : KTB site, and distance to the borehole.
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                in >> permeability_pressure_model_values[3];
                in >> permeability_pressure_model_values[4];
                in >> permeability_pressure_model_values[5];
                in >> permeability_pressure_model_values[6];
                in >> permeability_pressure_model_values[7];
                in >> permeability_pressure_model_values[8];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 8:       // Effective stress related to depth and curve, Urach
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                pcs_name_vector.push_back("PRESSURE1");
                break;
            case 9:       // Effective stress related to depth and curve, and to temperature, Urach
                in >> permeability_pressure_model_values[0];
                in >> permeability_pressure_model_values[1];
                in >> permeability_pressure_model_values[2];
                in >> permeability_pressure_model_values[3];
                in >> permeability_pressure_model_values[4];
                pcs_name_vector.push_back("PRESSURE1");
                pcs_name_vector.push_back("TEMPERATURE1");
                break;
            case 10: //WX:05.2010 directly from curve
                in >> permeability_pressure_model_values[0]; //WX: curve number 05.2010
                break;
            default:
                std::cout << "Error in MMPRead: no valid permeability model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //12.3 PERMEABILITY_FUNCTION_SATURATION
        //------------------------------------------------------------------------
        if(line_string.find("$PERMEABILITY_SATURATION") != std::string::npos)
        {
            int num_phases = 0;
            //if(H3_Process)
            //    num_phases = 3; // the existence of additional fluids in the .mfp file is not an indicator of a 4-phase system.
            //else if(H2_Process)
            //    num_phases = 2; // the existence of additional fluids in the .mfp file is not an indicator of a 3-phase system.
            //else if(H_Process)
            //    num_phases = 1;
            //
            for(k = 0; k < num_phases; k++)
            {
                in.str(readNonBlankLineFromInputStream(*mmp_file));
                in >> permeability_saturation_model[k];
                //
                switch(permeability_saturation_model[k])
                {
                case 0: // Curve of [Sw,k_rel]
                    in >> perm_saturation_value[k];    // curve number
                    break;
                //
                case 1: // Constant
                    in >> perm_saturation_value[k]; // constant value of k_rel
                    break;
                //
                case 2: // krg = 1.0 - krl (only for gas phase)
                    if(k==0){
                        //ScreenMessage("ERROR in MMPRead: Relative permeability model 2 is only valid for the gas phase.\n");
                        exit(0);
                    }
                    in >> minimum_relative_permeability;
                    break;
                //
                case 3: // JT: Function (power or linear)      WETTING
                    // krw = (m*Se)^exp
                    // Se  = (sl - slr) / (slm - slr)
                    in >> residual_saturation[k];            // slr: residual saturation, this phase
                    in >> maximum_saturation[k];            // slm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // (set exponent = 1 results in a linear k_rel function)
                    in >> perm_saturation_value[k];            // multiplier "m"
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 33: // JT: Function (power or linear)     NON-WETTING 
                    // krg = (m*(1-Se))^exp
                    // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgm --> slm = 1 - sgr
                    in >> residual_saturation[k];            // sgr: residual saturation, this phase
                    in >> maximum_saturation[k];            // sgm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // (set exponent = 1 results in a linear k_rel function)
                    in >> perm_saturation_value[k];            // multiplier "m"
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 4: // 2-phase Van Genuchten/Mualem Model  WETTING    
                    // krw = pow(se,0.5) * pow(1.0-pow(1.0-pow(se,1.0/m),m),2)
                    // Se  = (sl - slr) / (slm - slr)
                    in >> residual_saturation[k];            // slr: residual saturation, this phase
                    in >> maximum_saturation[k];            // slm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always <= 1.0) --> (typical is 0.5) i.e. n = 1 / (1 - exponent) == 2.0
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 44: // 2-phase Van Genuchten/Mualem Model NON-WETTING    
                    // krg = pow(1.0-se,1.0/3.0) * pow(1.0-pow(se,1.0/m),2.0*m)
                    // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgm --> slm = 1 - sgr
                    in >> residual_saturation[k];            // sgr: residual saturation, this phase
                    in >> maximum_saturation[k];            // sgm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always <= 1.0) --> (typical is 0.5) i.e. n = 1 / (1 - exponent) == 2.0
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 6: // Brooks/Corey                           WETTING        
                    // krw = pow(se,3.0+2.0/m)
                    // Se  = (sl - slr) / (slm - slr)
                    in >> residual_saturation[k];            // slr: residual saturation, this phase
                    in >> maximum_saturation[k];            // slm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always >= 1.0) (typical might be 2.0)
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 66: // Brooks/Corey                       NON-WETTING    
                    // krg = pow(1.0-se,2)*(1.0-pow(se,1.0+2.0/m))
                    // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgm --> slm = 1 - sgr
                    in >> residual_saturation[k];            // sgr: residual saturation, this phase
                    in >> maximum_saturation[k];            // sgm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always >= 1.0) (typical might be 2.0)
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 7: // Corey's curves                       WETTING        
                    // krw = pow(se,4)
                    // Se  = (sl - slr) / (slm - slr)
                    in >> residual_saturation[k];            // slr: residual saturation, this phase
                    in >> maximum_saturation[k];            // slm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always >= 1.0)
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                case 77: // Corey's curves                       NON-WETTING
                    // krg = pow((1.0-se),2)*(1.0-pow(se,2))
                    // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgmax --> slm = 1 - sgr
                    in >> residual_saturation[k];            // sgr: residual saturation, this phase
                    in >> maximum_saturation[k];            // sgm: maximum saturation, this phase
                    in >> saturation_exponent[k];            // exponent (always >= 1.0)
                    in >> minimum_relative_permeability;    // minimum relative permeability this phase
                    break;
                //
                default:
                    //ScreenMessage("Error in MMPRead: no valid permeability saturation model.\n");
                    break;
                }
                in.clear();
            }
            continue;
        }

        //------------------------------------------------------------------------
        //12.4 PERMEABILITY_FUNCTION_STRESS
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$PERMEABILITY_FUNCTION_STRESS") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_stress_mode;
            switch(permeability_stress_mode)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // k=const
                in >> permeability;
                break;
            case 2:       // Modified LBNL model. WW
                c_coefficient = new double [18];
                in >> c_coefficient[0] // b_0
                >> c_coefficient[1] // alpha
                >> c_coefficient[2] // br_1
                >> c_coefficient[3] // br_2
                >> c_coefficient[4] // br_3
                >> c_coefficient[5] >> ws; // fracture frequency
                for(i = 6; i < 18; i++)
                    c_coefficient[i] = 0.0;
                for(i = 0; i < 3; i++)
                {
                    in.clear();
                    in.str(readNonBlankLineFromInputStream(*mmp_file));
                    in >> m_string;
                    if(m_string.find("_XX") != string::npos)
                        k = 0;
                    else if(m_string.find("_YY") != string::npos)
                        k = 1;
                    else if(m_string.find("_ZZ") != string::npos)
                        k = 2;
                    m_string.clear();
                    in >> m_string;
                    pch = strtok (const_cast<char*> (m_string.c_str()),seps);
                    buff << pch;
                    buff >> c_coefficient[6 + k * 4];
                    buff.clear();
                    while (pch != NULL)
                    {
                        pch = strtok (NULL, seps);
                        if(pch == NULL)
                            break;
                        string token = pch;
                        tokens.push_back(token);
                    }
                    for(j = 0; j < (int)tokens.size(); j++)
                    {
                        pch = strtok (
                                const_cast<char*> (tokens[j].c_str()),seps1);
                        buff << pch;
                        buff >> f_buff;
                        buff.clear();
                        pch = strtok (NULL,seps1);
                        switch(pch[0])
                        {
                        case 'x':  c_coefficient[k * 4 + 7] = f_buff;
                            break;
                        case 'y':  c_coefficient[k * 4 + 8] = f_buff;
                            break;
                        case 'z':  c_coefficient[k * 4 + 9] = f_buff;
                            break;
                        }
                    }
                    tokens.clear();
                }
                break;
            case 3:       // Barton-Bandis  WW
                c_coefficient = new double [24];
                in >> c_coefficient[0] // JRC
                >> c_coefficient[1] // JCS        //an0
                >> c_coefficient[2] // UCS        //Kn
                >> c_coefficient[7] // sig_h
                >> c_coefficient[8] >> ws; // sig_H
                for(i = 3; i < 7; i++)
                    c_coefficient[i] = 0.0;
                for(i = 9; i < 24; i++)
                    c_coefficient[i] = 0.0;
                for(i = 0; i < 3; i++)
                {
                    in.clear();
                    in.str(readNonBlankLineFromInputStream(*mmp_file));
                    in >> m_string;
                    if(m_string.find("_XX") != string::npos)
                        k = 0;
                    else if(m_string.find("_YY") != string::npos)
                        k = 1;
                    else if(m_string.find("_ZZ") != string::npos)
                        k = 2;
                    m_string.clear();
                    in >> m_string;
                    pch = strtok (const_cast<char*> (m_string.c_str()),seps);
                    buff << pch;
                    buff >> c_coefficient[9 + k * 4];
                    buff.clear();
                    while (pch != NULL)
                    {
                        pch = strtok (NULL, seps);
                        if(pch == NULL)
                            break;
                        string token = pch;
                        tokens.push_back(token);
                    }
                    for(j = 0; j < (int)tokens.size(); j++)
                    {
                        pch = strtok (
                                const_cast<char*> (tokens[j].c_str()),seps1);
                        buff << pch;
                        buff >> f_buff;
                        buff.clear();
                        pch = strtok (NULL,seps1);
                        switch(pch[0])
                        {
                        case 'x':  c_coefficient[k * 4 + 10] = f_buff;
                            break;
                        case 'y':  c_coefficient[k * 4 + 11] = f_buff;
                            break;
                        case 'z':  c_coefficient[k * 4 + 12] = f_buff;
                            break;
                        }
                    }
                    tokens.clear();
                }
                //
                //CalStressPermeabilityFactor3_Coef();
                //
                break;

            default:
                std::cout << "Error in MMPRead: no valid permeability model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //12.5 PERMEABILITY_FUNCTION_VELOCITY
        //------------------------------------------------------------------------
        //WW
        if(line_string.find("$PERMEABILITY_FUNCTION_VELOCITY") != std::string::npos)
        {
            //WW   if(line_string.find("$PERMEABILITY_FUNCTION_STRESS")!=string::npos) { //subkeyword found
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_model;
            switch(permeability_model)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // k=const
                in >> permeability;
                break;
            default:
                std::cout << "Error in MMPRead: no valid permeability model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //12.6 PERMEABILITY_FUNCTION_POROSITY
        //------------------------------------------------------------------------

        //subkeyword found
        if(line_string.find("$PERMEABILITY_FUNCTION_POROSITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_porosity_model;
            switch(permeability_porosity_model)
            {
            case 0:       // k=f(x)
                break;
            case 1:       // k=const
                in >> permeability_porosity_model_values[0];
                mmp_file->ignore(MAX_ZEILE,'\n');
                break;
            case 2:       // Model from Ming Lian
                in >> permeability_porosity_model_values[0];
                in >> permeability_porosity_model_values[1];
                break;
            case 3:       // HS: 11.2008, Kozeny-Carman relationship
                // we set the tensor type first to isotropic and constant as initial value
                permeability_tensor_type = 0;
                permeability_model = 3; // this means permeability depends on K-C relationship
                                        // initial values
                in >> permeability_porosity_model_values[0];
                break;
            case 4:       // HS: 11.2008, Kozeny-Carman_normalized relationship
                // we set the tensor type first to isotropic and constant as initial value
                permeability_tensor_type = 0;
                permeability_model = 4; // this means permeability depends on K-C_normalized relationship
                                        // initial values
                in >> permeability_porosity_model_values[0];
                break;
            case 5:       // HS: 01.2010, Clement 1996 model
                // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
                permeability_tensor_type = 0;
                permeability_model = 5; // Clement original model
                                        // this is initial porosity
                in >> permeability_porosity_model_values[0];
                // this is initial permeability
                in >> permeability_porosity_model_values[1];
                break;
            case 6:       // HS: 01.2010, ,modified Clement, biomass colonies clogging
                // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
                permeability_tensor_type = 0;
                permeability_model = 6; // modified Clement, biomass growing in colonies
                                        // this is initial porosity
                in >> permeability_porosity_model_values[0];
                // this is initial permeability
                in >> permeability_porosity_model_values[1];
                // this is parameter a
                in >> permeability_porosity_model_values[2];
                break;
            case 7:       // HS: 01.2010, ,modified Clement, biofilm clogging
                // M. Thullner et al. 2004, J Contaminant Hydrology 70: 37-62, pp42
                permeability_tensor_type = 0;
                permeability_model = 7; // modified Clement, biomass growing in biofilm
                                        // this is initial porosity
                in >> permeability_porosity_model_values[0];
                // this is initial permeability
                in >> permeability_porosity_model_values[1];
                // this is parameter b
                in >> permeability_porosity_model_values[2];
                // this is prarameter k_fmin
                in >> permeability_porosity_model_values[3];
                break;
            default:
                std::cout << "Error in MMPRead: no valid permeability model" <<
                std::endl;
                break;
            }
            in.clear();
            continue;
        }

        //....................................................................
        //subkeyword found
        if(line_string.find("$CAPILLARY_PRESSURE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> capillary_pressure_model;
            i = 0;
            // Value needed to check for old version formats.
            capillary_pressure_values[2] = -1.0;
            bool old_format = false;
            //
            switch(capillary_pressure_model)
            {
            case 0:       // k=f(Se)
                in >> capillary_pressure_values[0]; // curve
                in >> capillary_pressure_values[1]; // Slr
                in >> capillary_pressure_values[2]; // Slmax
                //
                // JT: Check for old version format.
                if(capillary_pressure_values[2] < 0.0){
                    capillary_pressure_values[1] = residual_saturation[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[2] = maximum_saturation[0];    // old version uses relative permeabilty values for this
                    old_format = true;
                }
                break;
            case 1:       // const
                in >> capillary_pressure_values[0]; // the constant Pc value
                //
                // This next value is NOT REQUIRED. It is here because in the old version this entry was used to represent
                // a constant value of saturation in PP models. SHOULD DELETE THIS CHECK EVENTUALLY!
                in >> capillary_pressure_values[2];
                if(capillary_pressure_values[2] >= 0.0){ // Then a constant saturation value has been entered. This is model #2.
                    //ScreenMessage("WARNING in MMPRead. Capillary pressure model 1 used for a constant saturation. THIS IS NOW MODEL #2. PLEASE SWITCH TO MODEL #2.\n");
                    capillary_pressure_model = 2;
                    capillary_pressure_values[0] = capillary_pressure_values[2];
                }
                //
                // Assign bounds
                capillary_pressure_values[1] = 0.0; // Slr
                capillary_pressure_values[2] = 1.0; // Slmax
                break;
            case 2:          // Constant saturation for pp models (for WX, from JT)
                in >> capillary_pressure_values[0]; // The fixed saturation
                // Assign bounds
                capillary_pressure_values[1] = 0.0; // Slr
                capillary_pressure_values[2] = 1.0; // Slmax
                break;
            case 4:       // van Genuchten
                in >> capillary_pressure_values[0]; // Pb (or "alpha" if [alpha_switch>0])
                in >> capillary_pressure_values[1]; // Slr
                in >> capillary_pressure_values[2]; // Slmax
                in >> capillary_pressure_values[3]; // exponent (always <= 1.0) --> (typical is 0.5) i.e. n = 1 / (1 - exponent) == 2.0
                in >> capillary_pressure_values[4]; // maximum Pc
                in >> i;                            // alpha_switch (default = 0)
                if(i>0) entry_pressure_conversion = true;
                //
                // JT: Check for old version format.
                if(capillary_pressure_values[2] < 0.0){
                    entry_pressure_conversion = true;                        // entry is alpha in old version
                    capillary_pressure_values[1] = residual_saturation[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[2] = maximum_saturation[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[3] = saturation_exponent[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[4] = 1.0e10;
                    old_format = true;
                }
                break;
            case 6:       //  Brook & Corey. 2.05.2008. WW
                in >> capillary_pressure_values[0]; // Pb
                in >> capillary_pressure_values[1]; // Slr
                in >> capillary_pressure_values[2]; // Slmax
                in >> capillary_pressure_values[3]; // exponent (always >= 1.0) (typical might be 2.0)
                in >> capillary_pressure_values[4]; // maximum Pc
                //
                // JT: Check for old version format.
                if(capillary_pressure_values[2] < 0.0){
                    capillary_pressure_values[1] = residual_saturation[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[2] = maximum_saturation[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[3] = saturation_exponent[0];    // old version uses relative permeabilty values for this
                    capillary_pressure_values[4] = 1.0e10;
                    old_format = true;
                }
                break;
            default:
                ScreenMessage("Error in MMPRead: no valid capillary pressure model.\n");
                exit(1);
                break;
            }
            if(old_format){
                ScreenMessage("\n--\n Adopting capillary pressure saturation parameters from the\n");
                ScreenMessage(" relative permeability function for phase 0. Alternatively, you\n");
                ScreenMessage(" may enter capillary pressure specific parameters directly.\n--/n");
            }
            in.clear();
            continue;
        }
        //....................................................................
        //Dual Richards
        if(line_string.find("$TRANSFER_COEFFICIENT") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> transfer_coefficient; //(-)
            //   in >> unsaturated_hydraulic_conductivity;      //(L/T)
            in.clear();
            continue;
        }
        //....................................................................
        //Dual Richards
        if(line_string.find("$SPECIFIC_STORAGE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> specific_storage; //(Pa-1)
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //14 MASSDISPERSION_
        //            (1) LONGITUDINAL
        //            (2) TRANSVERSE
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$MASS_DISPERSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> mass_dispersion_model;
            switch(mass_dispersion_model)
            {
            case 0:       // f(x)
                break;
            case 1:       // Constant value
                in >> mass_dispersion_longitudinal;
                in >> mass_dispersion_transverse;
                in >> lgpn;
                if(lgpn > 0 )
                    cout << "      Limiting Grid Peclet Numbers to " << lgpn <<
                    endl;
                break;
            default:
                std::cout <<
                "Error in CMediumProperties::Read: no valid mass dispersion model"
                          << std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //15 HEATDISPERSION
        //            (1) LONGTIDUINAL
        //            (2) TRANSVERSE
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$HEAT_DISPERSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> heat_dispersion_model;
            switch(heat_dispersion_model)
            {
            case 0:       // f(x)
                break;
            case 1:       // Constant value
                in >> heat_dispersion_longitudinal;
                in >> heat_dispersion_transverse;
                break;
            default:
                std::cout <<
                "Error in CMediumProperties::Read: no valid heat dispersion model"
                          << std::endl;
                break;
            }
            in.clear();
            continue;
        }
        //subkeyword found
        if(line_string.find("$DIFFUSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> heat_diffusion_model;
            in.clear();
            continue;
        }
        //subkeyword found
        if(line_string.find("$EVAPORATION") != string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> evaporation;
            in >> heatflux;
            in >> vaporfraction;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        //16. Surface water
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$SURFACE_FRICTION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> friction_coefficient >> friction_exp_slope >> friction_exp_depth;
            in.clear();
            continue;
        }

        //subkeyword found
        if(line_string.find("$WIDTH") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> overland_width;
            in.clear();
            continue;
        }

        //subkeyword found
        if(line_string.find("$RILL") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> rill_height >> rill_epsilon;
            in.clear();
            continue;
        }

        //subkeyword found
        if(line_string.find("$CHANNEL") != std::string::npos)
        {
            channel = 1;
            continue;
        }

        //------------------------------------------------------------------------
        //19 ELECTRIC_CONDUCTIVITY
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //20 UNCONFINED_FLOW_GROUP
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //21 FLUID_EXCHANGE_WITH_OTHER_CONTINUA
        //------------------------------------------------------------------------

        //------------------------------------------------------------------------
        //11..PERMEABILITY_DISTRIBUTION
        //------------------------------------------------------------------------
//        size_t indexChWin, indexChLinux; //WW
//        indexChWin = indexChLinux = 0;
        std::string funfname;
        //subkeyword found
        if(line_string.find("$PERMEABILITY_DISTRIBUTION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> permeability_file;

            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        //11..POROSITY_DISTRIBUTION
        //------------------------------------------------------------------------
        //subkeyword found
        if(line_string.find("$POROSITY_DISTRIBUTION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*mmp_file));
            in >> porosity_file;
            string file_name = porosity_file;

            in.clear();
            continue;
        }

        //subkeyword found
        if(line_string.find("$FRACTURE") != std::string::npos)
        {
            is_fracture = true; //NW
            continue;
        }
    }
    return position;
}

}
