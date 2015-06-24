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

/**
 * \file FEM/Output.h
 * 05/04/2011 LB Refactoring: Moved from rf_out_new.h
 */

#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
#include <vector>

#include "DistributionInfo.h"
#include "ProcessInfo.h"

namespace ogs5
{

class COutput : public ProcessInfo, public DistributionInfo
{
public:
    COutput();
    COutput (size_t id);
    /**
     * method initializes process and mesh attributes
     */
    void init ();
    ~COutput(void);

    /**
     * read from file stream
     * @param in input file stream
     * @param geo_obj object of class GEOObjects that manages the geometric entities
     * @param unique_name the name of the project to access the right geometric entities
     * @return the new position in the stream after reading
     */
    std::ios::pos_type Read(std::ifstream& in);


    bool VARIABLESHARING;                         // Coordinates of each node as well as connection list is stored only for the first time step; BG: 05/2011 

    std::vector<std::string> _nod_value_vector;
    std::vector<std::string> _alias_nod_value_vector;
    // MAT values
    std::vector<std::string> mmp_value_vector; //OK
    std::vector<std::string> mfp_value_vector; //OK

    // MSH
    std::string msh_type_name;            //OK

    // TIM
    std::string tim_type_name;            // STEPS or TIMES ?
    std::vector<double> time_vector;
    double _time;

    /**
     * the position in the global vector out_vector, used only in NODWritePLYDataTEC
     */
    size_t _id;

    std::string file_base_name;
    double out_amplifier;                 //WW to amplify output
                                          //WW/OK

//    MeshLib::CFEMesh* m_msh;
    size_t nSteps;                           // After each nSteps, make output

    // GEO
    /**
     * the id of the geometric object as string REMOVE CANDIDATE
     */
    std::string geo_name;                 // TF 05/2010
    std::string geo_type;

    // File status
    bool _new_file_opened;                //WW

    // DAT
    /**
     * this attribute stores the output format
     */
    std::string dat_type_name;

    // ELE value
    std::vector<std::string> _ele_value_vector;

    // RWPT values
    std::vector<std::string> _rwpt_value_vector;

    // PCON values
    std::vector<std::string> _pcon_value_vector;
};

}
#endif // OUTPUT_H
