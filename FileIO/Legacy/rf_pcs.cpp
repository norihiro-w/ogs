/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_pcs.cpp
 *
 * Created on 2004-01-xx by Olaf Kolditz
 */

#include "rf_pcs.h"


// C++
#include <cfloat>
#include <iomanip>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include <logog/include/logog.hpp>

#include "makros.h"
#include "Ogs5FileTools.h"

using namespace std;

namespace ogs5
{

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   last modified:
   10/2010 TF changed process type handling from string to enum
**************************************************************************/
bool PCSRead(const std::string& file_base_name, std::vector<CRFProcess*> &pcs_vector)
{
    //----------------------------------------------------------------------
    char line[MAX_ZEILE];
    int indexCh1a, indexCh2a;
    std::string CommentK("//");
    std::string line_string;
    std::string pcs_file_name;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
    std::ifstream pcs_file(pcs_file_name.data(), ios::in);
    if (!pcs_file.good())
    {
        ERR("Warning: no PCS data *.pcs file is missing");
        return false;
    }

    // rewind the file
    pcs_file.clear();
    pcs_file.seekg(0, std::ios::beg);
    //========================================================================
    // Keyword loop
    //INFO("PCSRead ... ");
    while (!pcs_file.eof())
    {
        pcs_file.getline(line, MAX_ZEILE);
        line_string = line;
        line_string = GetUncommentedLine(line_string);
        if (line_string.find("#STOP") != string::npos)
            break;
        indexCh1a = (int) line_string.find_first_of(CommentK.c_str());
        indexCh2a = (int) line_string.find("#PROCESS");
        //----------------------------------------------------------------------
        // keyword found
        if (indexCh2a > indexCh1a && (indexCh1a == -1))
        {
            CRFProcess* m_pcs = new CRFProcess();
            position = m_pcs->Read(&pcs_file);
            m_pcs->pcs_number = pcs_vector.size();

            pcs_vector.push_back(m_pcs);

            pcs_file.seekg(position, std::ios::beg);
        }                         // keyword found
    }                                     // eof

    //INFO("done, read %d processes", pcs_vector.size());

    return true;
}

CRFProcess::CRFProcess()
: Deactivated_SubDomain(0), mesh_id(0), timegroup_id(0)
{

}

CRFProcess::~CRFProcess()
{
    if (Deactivated_SubDomain!=NULL)
        delete [] Deactivated_SubDomain;
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Read deformation process
   11/2004 OK file streaming
   12/2005 OK MSH_TYPE
   01/2006 OK GEO_TYPE
**************************************************************************/
std::ios::pos_type CRFProcess::Read(std::ifstream* pcs_file)
{
    char line[MAX_ZEILE];
    string line_string;
    string CommentK("//");
    string hash("#");
    bool new_keyword = false;
    bool new_subkeyword = false;
    ios::pos_type position;
    ios::pos_type position_subkeyword;
    std::stringstream line_stream;
    saturation_switch = false;            // JOD for Richards
    //----------------------------------------------------------------------
    while (!new_keyword)
    {
        position = pcs_file->tellg();
        pcs_file->getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find(hash) != string::npos)
        {
            new_keyword = true;
            break;
        }
        //....................................................................
        // WW Comment line
        if (line_string.find_first_of(CommentK.c_str()) != string::npos)
            return position;
        //SB check for comment sign ;
        line_string = GetUncommentedLine(line_string);
        //....................................................................
        // subkeyword found
        if (line_string.find("$PCS_TYPE") != string::npos)
            while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
            {
                position = pcs_file->tellg();
                line_string = readNonBlankLineFromInputStream(*pcs_file);
                if (line_string.find("#") != string::npos)
                    return position;
                if (line_string.find("$") != string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                line_stream.str(line_string);
                std::string pcs_type_name;
                line_stream >> pcs_type_name;
                pcs_type_name_vector.push_back(pcs_type_name);
                this->setProcessType (FiniteElement::convertProcessType(pcs_type_name));
                line_stream.clear();

                if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
                {
                    pcs_no_components++;
                    this->setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
                }
            }
        //....................................................................
        // subkeyword found
        if (line_string.find("$NUM_TYPE") != string::npos)
        {
            *pcs_file >> num_type_name;
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$CPL_TYPE") != string::npos)
        {
            *pcs_file >> cpl_type_name;
            if (cpl_type_name.compare("MONOLITHIC") == 0)
            {
                pcs_monolithic_flow = true;
            }
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$TIM_TYPE") != string::npos)
        {
            *pcs_file >> tim_type_name;
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$APP_TYPE") != string::npos)
        {
            *pcs_file >> rwpt_app;
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
        {
			std::string pri_var_str_tmp;
            *pcs_file >> pri_var_str_tmp;  
			primary_variable_name.push_back(pri_var_str_tmp);
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$ELEMENT_MATRIX_OUTPUT") != string::npos)
        {
            *pcs_file >> Write_Matrix; //WW
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        //WW
        if (line_string.find("$BOUNDARY_CONDITION_OUTPUT") != string::npos)
        {
            write_boundary_condition = true;
            continue;
        }
        //....................................................................
        //BG 05/2012
        if (line_string.find("$OutputMassOfGasInModel") != string::npos)
        {
            OutputMassOfGasInModel = true;
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$ST_RHS") != string::npos)
        {
            *pcs_file >> WriteSourceNBC_RHS; //WW
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        if (line_string.find("$PROCESSED_BC") != string::npos) //25.08.2011. WW
        {
            *pcs_file >> WriteProcessed_BC;
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }

        //....................................................................
        // subkeyword found
        if (line_string.find("$MEMORY_TYPE") != string::npos)
        {
            *pcs_file >> Memory_Type; //WW
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$RELOAD") != string::npos)
        {
            *pcs_file >> reload; //WW
            if (reload == 1 || reload == 3)
                *pcs_file >> nwrite_restart;  //kg44 read number of timesteps between writing restart files
            pcs_file->ignore(MAX_ZEILE, '\n');
            continue;
        }
        // subkeyword found
        if (line_string.find("$DEACTIVATED_SUBDOMAIN") != string::npos)
        {
            //WW
            *pcs_file >> NumDeactivated_SubDomains >> ws;
            Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
            for (int i = 0; i < NumDeactivated_SubDomains; i++)
                *pcs_file >> Deactivated_SubDomain[i] >> ws;
            continue;
        }
        //....................................................................
        // subkeyword found
        if (line_string.find("$MSH_TYPE") != string::npos)
        {
            *pcs_file >> msh_type_name >> ws;
            continue;
        }
        //....................................................................
        //        if (line_string.find("$GEO_TYPE") != string::npos) { //OK
        //            *pcs_file >> geo_type >> geo_type_name >> ws;
        //            continue;
        //        }
        //
        //....................................................................
        // subkeyword found
        if (line_string.find("$MEDIUM_TYPE") != string::npos)
        {
            while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
            {
                position_subkeyword = pcs_file->tellg();
                *pcs_file >> line_string;
                if (line_string.size() == 0)
                    break;
                if (line_string.find("#") != string::npos)
                {
                    new_keyword = true;
                    break;
                }
                if (line_string.find("$") != string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                if (line_string.find("CONTINUUM") != string::npos)
                {
                    *pcs_file >> line_string;
                    //WW
                    double w_m = strtod(line_string.data(), NULL);
                    continuum_vector.push_back(w_m);
                    //WW
                    continuum_vector.push_back(1.0 - w_m);
                    break; //WW
                }
                pcs_file->ignore(MAX_ZEILE, '\n');
            }
            continue;
        }
        //OK
        if (line_string.find("$SATURATION_SWITCH") != string::npos)
        {
            saturation_switch = true;
            continue;
        }
        //SB4900
        if (line_string.find("$USE_VELOCITIES_FOR_TRANSPORT") != string::npos)
        {
            //// Only for fluid momentum process
            //if (this->getProcessType () == FiniteElement::FLUID_MOMENTUM)
            //    use_velocities_for_transport = true;
            continue;
        }
        //Interface to Eclipse and Dumux, BG, 09/2010
        //    if(line_string.find("$SIMULATOR")!=string::npos) { //OK
        if(line_string.compare("$SIMULATOR") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
        {
            *pcs_file >> simulator;
            continue;
        }
        if(line_string.find("$SIMULATOR_PATH") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
        {
            *pcs_file >> simulator_path;
            continue;
        }
        // BG, 09/2010, coupling to Eclipse and DuMux
        if(line_string.find("$SIMULATOR_MODEL_PATH") == 0)
        {
            *pcs_file >> simulator_model_path;
            continue;
        }
        // BG, 09/2010, coupling to Eclipse and DuMux
        if(line_string.find("$USE_PRECALCULATED_FILES") == 0)
        {
            PrecalculatedFiles = true;
            continue;
        }
        // KB, 02/2011, coupling to Eclipse and DuMux
        if(line_string.find("$SIMULATOR_WELL_PATH") == 0)
        {
            *pcs_file >> simulator_well_path;
            continue;
        }
        // BG, NB 11/2010, calculating phase transition for CO2
        if(line_string.find("$PHASE_TRANSITION") == 0)
        {
            string tempstring;
            *pcs_file >> tempstring;
            //if (tempstring == "CO2_H2O_NaCl")
            //    this->Phase_Transition_Model = 1;
            continue;
        }
        //WX:07.2011
        if(line_string.find("$TIME_CONTROLLED_EXCAVATION") == 0)
        {
            *pcs_file >> ExcavMaterialGroup >> ExcavDirection >>
            ExcavBeginCoordinate >> ExcavCurve;
            continue;
        }
        //....................................................................
    }
    //----------------------------------------------------------------------
    return position;
}
}
