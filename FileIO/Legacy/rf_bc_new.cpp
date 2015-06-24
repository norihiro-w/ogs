/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_bc_new.cpp
 *
 * Created on 2004-02-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Class: BC BoundaryConditions
   Task:
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#include "rf_bc_new.h"

// C++ STL
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <list>
#include <sstream>
#include <string>

#include "makros.h"
//// FileIO
#include "ProcessIO.h"
#include "Ogs5FileTools.h"

namespace ogs5
{

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   09/2004 OK POINTS method
   11/2004 MX stream string
**************************************************************************/
std::ios::pos_type CBoundaryCondition::Read(std::ifstream* bc_file,
                                            bool & valid)
{
    std::string line_string;
    bool new_keyword = false;
    std::ios::pos_type position;

    std::string sub_string, strbuff;
    int ibuff;                            //pos,
    double dbuff;                         //WW
    std::stringstream in;

    std::string FilePath;

    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        position = bc_file->tellg();
        line_string = readNonBlankLineFromInputStream (*bc_file);
        if (line_string.size() < 1)
            break;
        if (line_string.find("#") != std::string::npos)
        {
            new_keyword = true;
            break;
        }

        if (line_string.find("$PCS_TYPE") != std::string::npos)
            if (!FileIO::ProcessIO::readProcessInfo (*bc_file, _pcs_type))
                valid = false;

        if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> primaryvariable_name;    // _pcs_pv_name;
            in.clear();
        }

        // HS, this is new. later on we should stick to COMP_NAME, PRIMARY_VARIABLE support will be removed.
        if (line_string.find("$COMP_NAME") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            std::string tmp;
            in >> tmp;    // _pcs_pv_name;
            in.clear();
        }

        if (line_string.find("$GEO_TYPE") != std::string::npos)
            //if (!FileIO::GeoIO::readGeoInfo (this, *bc_file, geo_name, geo_obj,
            //                                 unique_fname))
            //    valid = false;
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> geo_type_name >> geo_name;
            in.clear();
        }

        //PCH
        if (line_string.find("$DIS_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> line_string; //sub_line
            _periodic = false; // JOD

            // Source terms are assign to element nodes directly. 23.02.2009. WW
            if (line_string.find("DIRECT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::DIRECT);
                in >> fname;
                fname = FilePath + fname;
                in.clear();
            }

            if (line_string.find("CONSTANT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> geo_node_value; //sub_line
                in.clear();
            }
            // If a linear function is given. 25.08.2011. WW
            if (line_string.find("FUNCTION") != std::string::npos)
            {
                setProcessDistributionType(FiniteElement::FUNCTION);
                in.clear();
                in.str(readNonBlankLineFromInputStream(*bc_file));
                in >> function_exp; //sub_line
                in.clear();
            }
            if (line_string.find("LINEAR") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::LINEAR);
                // Distribuded. WW
                size_t nLBC;
                in >> nLBC; //sub_line
                in.clear();

                for (size_t i = 0; i < nLBC; i++)
                {
                    in.str(readNonBlankLineFromInputStream(*bc_file));
                    in >> ibuff >> dbuff >> strbuff;
                    in.clear();

                    //           *bc_file>>ibuff>>dbuff;
                    _PointsHaveDistribedBC.push_back(ibuff);
                    _DistribedBC.push_back(dbuff);
                    if (strbuff.size() > 0)
                    {
                        _PointsFCTNames.push_back(strbuff);
                        time_dep_interpol = true;
                    }
                }
                //        bc_file->ignore(MAX_ZEILE,'\n');
            }
        }

        // Time dependent function
        //..Time dependent curve ............................................
        if (line_string.find("$TIM_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> line_string;

            if (line_string.find("CURVE") != std::string::npos)
            {
                //                tim_type_name = "CURVE";
                this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> _curve_index;
                in.clear();

                //        pos1=pos2+1;
                //        sub_string = get_sub_string(buffer,"  ",pos1,&pos2);
                //        _curve_index = atoi(sub_string.c_str());
            }
            continue;
        }

        if (line_string.find("$FCT_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> fct_name; //sub_line
            in.clear();
        }

        if (line_string.find("$MSH_TYPE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> sub_string; //sub_line
            _msh_type_name = "NODE";
            if (sub_string.find("NODE") != std::string::npos)
            {
                in >> _msh_node_number;
                in.clear();
            }
        }

        if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file)); // CONSTANT -21500.0
            in >> line_string;
            if (line_string.find("CONSTANT") != std::string::npos)
            {
                this->setProcessDistributionType(FiniteElement::CONSTANT);
                in >> geo_node_value;
                in.clear();
            }
            in.str(readNonBlankLineFromInputStream(*bc_file)); // 0.0 IF HEAD > 0.04
            std::string pcs_pv_name_cond; // 07/2010 TF temp string
            in >> node_value_cond >> line_string >> pcs_pv_name_cond
            >> line_string >> condition;
            in.clear();
            in.str(readNonBlankLineFromInputStream(*bc_file)); // PCS OVERLAND_FLOW
            std::string pcs_type_name_cond;
            in >> line_string >> pcs_type_name_cond;
            in.clear();
            conditional = true;
        }

        if (line_string.find("$EPSILON") != std::string::npos) // NW
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> epsilon;
            in.clear();
        }
        //....................................................................
        //aktive state of the bc is time controlled  WX
        if (line_string.find("$TIME_CONTROLLED_ACTIVE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> time_contr_curve;
            in.clear();
        }
        //....................................................................
        //bc for excated boundaries WX
        if (line_string.find("$EXCAVATION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*bc_file));
            in >> bcExcav >> MatGr;
            in.clear();
        }
        //....................................................................
    }
    return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: BC constructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::CBoundaryCondition() :
    geo_name (""), _curve_index (-1)
{
    this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
    // FCT
    conditional = false;
    time_dep_interpol = false;
    epsilon = 1e-9;                       //NW
    time_contr_curve = -1;                //WX
    bcExcav = -1;                         //WX
    MatGr = -1;                           //WX
    _msh_node_number = 0;
    _periode_phase_shift = .0;
    _periode_time_length = .0;
    _periodic = false;
    condition = .0;
    geo_node_number = 0;
    geo_node_value = .0;
    node_value_cond = .0;
}



/**************************************************************************
   FEMLib-Method:
   Task: BC deconstructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::~CBoundaryCondition()
{
    // DIS
    node_number_vector.clear();
    geo_node_number = -1;
    geo_node_value = 0.0;

}


/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   05/2010 TF changes due to new GEOLIB integration, some improvements
**************************************************************************/
bool BCRead(std::string const& file_base_name, 
        std::vector<CBoundaryCondition*>& bc_vector)
{
    char line[MAX_ZEILE];
    std::string line_string, bc_file_name;

    // File handling
    bc_file_name = file_base_name + BC_FILE_EXTENSION;

    std::ifstream bc_file(bc_file_name.data(), std::ios::in);
    if (!bc_file.good())
    {
        std::cout << "! Error in BCRead: No boundary conditions !" << std::endl;
        return false;
    }

    // Keyword loop
    std::cout << "BCRead ... " << std::flush;
    while (!bc_file.eof())
    {
        bc_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != std::string::npos)
        {
            std::cout << "done, read " << bc_vector.size()
                      << " boundary conditions" << std::endl;
            return true;
        }
        if (line_string.find("#BOUNDARY_CONDITION") != std::string::npos)
        {
            CBoundaryCondition* bc(new CBoundaryCondition());
            bool valid (true);
            std::ios::pos_type position = bc->Read(&bc_file, valid);
            if (valid)
                bc_vector.push_back(bc);
            else
                delete bc;
            bc_file.seekg(position, std::ios::beg);
        } // keyword found
    } // eof
    return true;
}

}
