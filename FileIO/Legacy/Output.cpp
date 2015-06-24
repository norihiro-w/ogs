/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Output.cpp
 *
 * Created on 2011-04-05 by Lars Bilke
 */

#include "Output.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Ogs5FileTools.h"
#include "makros.h"

using namespace std;

namespace ogs5
{

COutput::COutput() :
    ProcessInfo(), tim_type_name("TIMES"), _id(0), out_amplifier(0.0),
    nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
    VARIABLESHARING = false;    //BG
    _time = .0;
}

COutput::COutput(size_t id) :
    ProcessInfo(), tim_type_name("TIMES"), _id(id), out_amplifier(0.0),
    nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
    VARIABLESHARING = false;    //BG
    _time = .0;
}

COutput::~COutput()
{
    mmp_value_vector.clear();             //OK
}


#define KEYWORD '#'
#define SUBKEYWORD '$'
bool Keyword(const std::string &line)
{
    if(line.find(KEYWORD) != std::string::npos)
        return true;
    else
        return false;
}
/**************************************************************************
   STRLib-Method: SubKeyword
   Task:
   Programing:
   09/2004 OK Implementation
   last modification:
**************************************************************************/
bool SubKeyword(const std::string &line)
{
    if(line.find(SUBKEYWORD) != std::string::npos)
        return true;
    else
        return false;
}

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   07/2004 WW Remove old files
   11/2004 OK string streaming by SB for lines
   03/2005 OK PCS_TYPE
   12/2005 OK DIS_TYPE
   12/2005 OK MSH_TYPE
   08/2008 OK MAT
   06/2010 TF formated, restructured, signature changed, use new GEOLIB data structures
   09/2010 TF signature changed, removed some variables
**************************************************************************/
ios::pos_type COutput::Read(std::ifstream& in_str)
{
    std::string line_string;
    bool new_keyword = false;
    ios::pos_type position;
    bool new_subkeyword = false;
    ios::pos_type position_line;
    bool ok = true;
    string name;
    ios::pos_type position_subkeyword;
    std::stringstream in;

    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        position = in_str.tellg();
        if (new_subkeyword)
            in_str.seekg(position_subkeyword, ios::beg);
        new_subkeyword = false;
        // SB new input        in_str.getline(buffer,MAX_ZEILE);
        // SB new input         line_string = buffer;
        line_string.clear();
        line_string = readNonBlankLineFromInputStream(in_str);
        if (line_string.size() < 1)
            break;

        if (Keyword(line_string))
            return position;

        // subkeyword found
        if (line_string.find("$NOD_VALUES") != string::npos)
        {
            while ((!new_keyword) && (!new_subkeyword))
            {
                position_subkeyword = in_str.tellg();
                //SB input with comments  in_str >> line_string>>ws;
                line_string = readNonBlankLineFromInputStream(in_str);
                if (line_string.find("#") != string::npos)
                    return position;
                if (line_string.find("$") != string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                if (line_string.size() == 0)
                    break;  //SB: empty line
                in.str(line_string);
                in >> name;
                //_alias_nod_value_vector.push_back(name);
                _nod_value_vector.push_back(name);
                in.clear();
            }

            continue;
        }
        //--------------------------------------------------------------------
        // subkeyword found //MX
        if (line_string.find("$PCON_VALUES") != string::npos)
        {
            while ((!new_keyword) && (!new_subkeyword))
            {
                position_subkeyword = in_str.tellg();
                in_str >> line_string >> ws;
                if (line_string.find("#") != string::npos)
                    return position;
                if (line_string.find("$") != string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                if (line_string.size() == 0)
                    break;
                _pcon_value_vector.push_back(line_string);
            }
            continue;
        }

        //--------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$ELE_VALUES") != string::npos)
        {
            ok = true;
            while (ok)
            {
                position_line = in_str.tellg();
                in_str >> line_string;
                if (SubKeyword(line_string))
                {
                    in_str.seekg(position_line, ios::beg);
                    ok = false;
                    continue;
                }
                if (Keyword(line_string))
                    return position;
                _ele_value_vector.push_back(line_string);
                in_str.ignore(MAX_ZEILE, '\n');
            }
            /*
               // Commented by WW
               // Remove files
               tec_file_name = file_base_name + "_domain_ele" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
             */
            continue;
        }
        //-------------------------------------------------------------------- // Added 03.2010 JTARON
        // subkeyword found
        if (line_string.find("$RWPT_VALUES") != string::npos)
        {
            while ((!new_keyword) && (!new_subkeyword))
            {
                position_subkeyword = in_str.tellg();
                line_string = readNonBlankLineFromInputStream(in_str);
                if (line_string.find("#") != string::npos)
                    return position;
                if (line_string.find("$") != string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                if (line_string.size() == 0)
                    break;  //SB: empty line
                in.str(line_string);
                in >> name;
                _rwpt_value_vector.push_back(name);
                in.clear();
            }
            continue;
        }

        //subkeyword found
        if (line_string.find("$GEO_TYPE") != string::npos)
        {
            in_str >> geo_type;
            if (geo_type.find("DOMAIN")==std::string::npos)
                in_str >> geo_name;
            in_str.ignore(MAX_ZEILE, '\n');
            continue;
        }

        // subkeyword found
        if (line_string.find("$TIM_TYPE") != string::npos)
        {
            while ((!new_keyword) && (!new_subkeyword))
            {
                position_subkeyword = in_str.tellg();
                in_str >> line_string;
                if (line_string.size() == 0) //SB
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
                if (line_string.find("STEPS") != string::npos)
                {
                    in_str >> nSteps;
                    tim_type_name = "STEPS"; //OK
                    break; //kg44 I guess that was missing..otherwise it pushes back a time_vector!
                }
                // JT 2010, reconfigured (and added RWPT)... didn't work
                if (line_string.find("STEPPING") != string::npos)
                {
                    double stepping_length, stepping_end, stepping_current;
                    in_str >> stepping_length >> stepping_end;
                    stepping_current = stepping_length;
                    while (stepping_current <= stepping_end)
                    {
                        time_vector.push_back(stepping_current);
                        //                        rwpt_time_vector.push_back(stepping_current);
                        stepping_current += stepping_length;
                    }
                }
                else
                    time_vector.push_back(strtod(line_string.data(), NULL));
                //                    rwpt_time_vector.push_back(strtod(line_string.data(), NULL));
                in_str.ignore(MAX_ZEILE, '\n');
            }
            continue;
        }

        // subkeyword found
        if (line_string.find("$DAT_TYPE") != string::npos)
        {
            in_str >> dat_type_name;
            in_str.ignore(MAX_ZEILE, '\n');
            continue;
        }

        // Coordinates of each node as well as connection list is stored only for the first time step; BG: 05/2011 
        if (line_string.find("$VARIABLESHARING") != string::npos)
        {
           this->VARIABLESHARING = true;
           continue;
        }  

        // subkeyword found
        if (line_string.find("$AMPLIFIER") != string::npos)
        {
            in_str >> out_amplifier;
            in_str.ignore(MAX_ZEILE, '\n');
            continue;
        }

        // subkeyword found
        if (line_string.find("$PCS_TYPE") != string::npos)
        {
            std::string tmp_pcs_type_name;
            in_str >> tmp_pcs_type_name;
            setProcessType(FiniteElement::convertProcessType(tmp_pcs_type_name));
            in_str.ignore(MAX_ZEILE, '\n');
            /* // Comment by WW
               // Remove files
               tec_file_name = pcs_type_name + "_" + "domain" + "_line" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
               tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
               tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
               tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
               tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
               tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
               remove(tec_file_name.c_str());
             */
            continue;
        }

        // subkeyword found
        if (line_string.find("$DIS_TYPE") != string::npos)
        {
            std::string dis_type_name;
            in_str >> dis_type_name;
            setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));
            in_str.ignore(MAX_ZEILE, '\n');
            continue;
        }

        // subkeyword found
        if (line_string.find("$MSH_TYPE") != string::npos)
        {
            in_str >> msh_type_name;
            in_str.ignore(MAX_ZEILE, '\n');
            continue;
        }

        //OK
        if (line_string.find("$MMP_VALUES") != string::npos)
        {
            ok = true;
            while (ok)
            {
                position_line = in_str.tellg();
                in_str >> line_string;
                if (SubKeyword(line_string))
                {
                    in_str.seekg(position_line, ios::beg);
                    ok = false;
                    continue;
                }
                if (Keyword(line_string))
                    return position;
                mmp_value_vector.push_back(line_string);
                in_str.ignore(MAX_ZEILE, '\n');
            }
            continue;
        }

        //OK
        if (line_string.find("$MFP_VALUES") != string::npos)
        {
            ok = true;
            while (ok)
            {
                position_line = in_str.tellg();
                in_str >> line_string;
                if (SubKeyword(line_string))
                {
                    in_str.seekg(position_line, ios::beg);
                    ok = false;
                    continue;
                }
                if (Keyword(line_string))
                    return position;
                mfp_value_vector.push_back(line_string);
                in_str.ignore(MAX_ZEILE, '\n');
            }

            continue;
        }
    }
    return position;
}

}
