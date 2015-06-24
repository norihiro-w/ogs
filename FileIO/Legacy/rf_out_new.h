/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_out_new.h
 *
 * Created on 2004-06-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: OUT
   Task: class implementation
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_out_new_INC
#define rf_out_new_INC

#include <vector>
#include <string>

namespace ogs5
{

class COutput;

/**
 * read file that stores information about output
 * @param file_base_name base file name (without extension)
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if file reading was successful, else false
 */
bool OUTRead(const std::string& file_base_name,
         std::vector<COutput*> &out_vector);

#define OUT_FILE_EXTENSION ".out"

}

#endif
