/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_ic_new.h
 *
 * Created on 2004-08-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: Initial Conditions IC
   Task: class implementation
   Programing:
   08/2004 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_ic_new_INC
#define rf_ic_new_INC

#define IC_FILE_EXTENSION ".ic"

// C++ STL
#include <fstream>
#include <string>
#include <vector>

// FEM
#include "DistributionInfo.h"                     // TF
#include "ProcessInfo.h"                          // KR

namespace ogs5
{

/**
 * class for handling initial conditions
 */
class CInitialCondition : public ProcessInfo, public DistributionInfo
{
public:
    CInitialCondition();
    ~CInitialCondition();

    /**
     * read initial condition from stream
     * @param in input stream from file
     * @return the new position in the stream after reading
     */
    std::ios::pos_type Read(std::ifstream* in);

public:
    size_t SubNumber;                     //WW
    std::vector<int> subdom_index;        //WW
    std::vector<double> subdom_ic;        //WW
    std::string fname;                    //17.11.2009. PCH

    std::string primaryvariable_name;

    std::string geo_type_name;
    std::string geo_name;                 // TF 05/2010
    double geo_node_value;              //KR

    double gradient_ref_depth;
    double gradient_ref_depth_value;
    double gradient_ref_depth_gradient;
    std::string rfr_file_name;            //OK
};


bool ICRead(const std::string& file_base_name, std::vector<CInitialCondition*> &ic_vector);
}

#endif
