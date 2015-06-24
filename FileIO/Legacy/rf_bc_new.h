/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_bc_new.h
 *
 * Created on 2004-02-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: Boundary Conditions
   Task: class implementation
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_bc_new_INC
#define rf_bc_new_INC

#include <vector>
#include <iostream>


#include "DistributionInfo.h" // TF
#include "ProcessInfo.h"                          // KR


namespace ogs5
{

namespace FileIO
{
class BoundaryConditionIO;
}


class CBoundaryCondition :
    public ProcessInfo,
    public DistributionInfo
{
public:
    friend class CBoundaryConditionsGroup;
    friend class FileIO::BoundaryConditionIO;
    CBoundaryCondition();

    ~CBoundaryCondition();

    /**
     * reads a boundary condition from stream
     * @param in input file stream for reading
     * @param geo_obj pointer to the geometric object manager
     * @param unique_fname the project name
     * @param valid after return the variable valid contains the status of the object,
     * valid is false if there occured an error while reading the data, else true
     * @return the position in the stream after the boundary condition
     */
    // TF
    std::ios::pos_type Read(std::ifstream* in,
                            bool &valid);


public:
    std::vector<std::string> _PointsFCTNames;
    std::vector<size_t> _PointsHaveDistribedBC;
    std::vector<double> _DistribedBC;

    // GEO
    /**
     * the id of the geometric object as string REMOVE CANDIDATE
     */
    std::string geo_name; // TF 05/2010
    std::string geo_type_name;
    std::string primaryvariable_name;

    std::string fname; //27.02.2009. WW
    int _curve_index; // Time function index

    // DIS
    std::vector<long> node_number_vector;
    std::vector<double> node_value_vector;
    long geo_node_number;
    double geo_node_value;

    double _periode_phase_shift; // JOD
    double _periode_time_length; // JOD
    bool _periodic; // JOD

    double node_value_cond; //OK
    double condition; //OK
    double epsilon; //NW. temporally set here for surface interpolation
    bool time_dep_interpol;

    // FCT
    std::string fct_name;
    bool conditional;

    // MSH
    long _msh_node_number;
    std::string _msh_type_name; //OK4105

    // Excavation WX:12.2010
    int bcExcav;
    int MatGr;
    // aktive state is controlled by time curve WX:01.2011
    int time_contr_curve;

    std::string function_exp;
};


//========================================================================
#define BC_FILE_EXTENSION ".bc"

bool BCRead (std::string const& file_base_name,
        std::vector<CBoundaryCondition*> &bc_vector);
}

#endif
