/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_out_new.cpp
 *
 * Created on 2004-06-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: OUT
   Task:
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
#include "rf_out_new.h"

// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
using namespace std;

#include "makros.h"
#include "Output.h"

namespace ogs5
{


/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Remove the old files
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   06/2006 WW Remove the old files by new way
   06/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool OUTRead(const std::string& file_base_name,
        std::vector<COutput*> &out_vector)
{
    char line[MAX_ZEILE];
    std::string line_string;
    ios::pos_type position;

    // File handling
    std::string out_file_name = file_base_name + OUT_FILE_EXTENSION;
    std::ifstream out_file(out_file_name.data(), ios::in);
    if (!out_file.good())
        return false;
    out_file.seekg(0L, ios::beg);

    // Keyword loop
    cout << "OUTRead ... " << std::flush;
    while (!out_file.eof())
    {
        out_file.getline(line, MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos) {
            std::cout << "done, read " << out_vector.size() << " output settings" <<
            std::endl;
           return true;
        }

//        //15.01.2008. WW
//        if (line_string.find("#VERSION") != string::npos)
//            output_version = true;  // 02.2011. WW
        //----------------------------------------------------------------------
        // keyword found
        if (line_string.find("#OUTPUT") != string::npos)
        {
            COutput* out(new COutput(out_vector.size()));
            position = out->Read(out_file);
            out_vector.push_back(out);
            out_file.seekg(position, ios::beg);
        }                         // keyword found
    }                                     // eof
    return true;
}

}
