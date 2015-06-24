/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_pcs.h
 *
 * Created on 2004-01-xx by Olaf Kolditz
 */

#ifndef rf_pcs_INC
#define rf_pcs_INC

#include <vector>
#include <iostream>
#include <string>

#include "ProcessInfo.h"

#define PCS_FILE_EXTENSION ".pcs"

namespace ogs5
{

/**
 * class manages the physical processes
 */
class CRFProcess : public ProcessInfo
{
    friend bool PCSRead(const std::string&, std::vector<CRFProcess*>&);
public:
    CRFProcess(void);
    virtual ~CRFProcess();

    std::ios::pos_type Read(std::ifstream*);

public:
    size_t pcs_no_components;
    bool pcs_monolithic_flow;
    size_t pcs_number;
    std::vector<double> continuum_vector;
    double PrecalculatedFiles;
    std::vector<int> pcs_number_mass;                        // JT2012
    std::vector<std::string> pcs_type_name_vector;
    int Memory_Type;
    int NumDeactivated_SubDomains;
    int* Deactivated_SubDomain;
    bool write_boundary_condition;        //15.01.2008. WW
    bool OutputMassOfGasInModel;            // BG 05/2012
    bool Write_Matrix;
    int WriteSourceNBC_RHS;
    int WriteProcessed_BC;
    int reload;
    long nwrite_restart;
    std::string geo_type;                 //OK
    std::string geo_type_name;            //OK
    std::vector<long> bc_node_value_in_dom; //WW for domain decomposition
    std::vector<long> bc_local_index_in_dom; //WW for domain decomposition
    std::vector<long> rank_bc_node_value_in_dom; //WW
    std::vector<long> bc_transient_index; //WW/CB
    std::vector<long> st_transient_index;       //WW/CB...BG
    std::vector<long> st_node_value_in_dom; //WW for domain decomposition
    std::vector<long> st_local_index_in_dom; //WW for domain decomposition
    std::vector<long> rank_st_node_value_in_dom; //WW
    std::vector<long> stgem_node_value_in_dom; //KG for domain decomposition
    std::vector<long> stgem_local_index_in_dom; //KG for domain decomposition
    std::vector<long> rank_stgem_node_value_in_dom;
    std::string simulator;                // which solver to use, i.e. GeoSys, ECLIPSE or DuMux
    std::string simulator_path;           // path for executable of external simulator
    std::string simulator_model_path;     // path to exclipse input data file (*.data), with extension
    std::string simulator_well_path;      // path to well schedule ( *.well), with extension
    std::string file_name_base;           //OK
    std::vector<std::string> primary_variable_name;    //OK // HS
    std::string num_type_name;
    int rwpt_app;
    std::string tim_type_name;            //OK
    std::string cpl_type_name;
    bool saturation_switch;               // JOD
    std::string msh_type_name;            //OK
    std::vector<double*> nod_val_vector;  //OK
    std::vector<std::string> nod_val_name_vector;
    std::vector<std::string> ele_val_name_vector;
    std::vector<double*> ele_val_vector;  //PCH
    int ExcavMaterialGroup;               //WX
    int ExcavDirection;                   //WX
    int ExcavCurve;                       //WX
    double ExcavBeginCoordinate;          //WX

    size_t mesh_id;
    size_t timegroup_id;
};

bool PCSRead(const std::string&, std::vector<CRFProcess*> &pcs_vector);


}

#endif
