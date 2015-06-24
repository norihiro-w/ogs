/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_num_new.h
 *
 * Created on 2004-11-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: NUM
   Task: class implementation
   Programing:
   11/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef rf_num_new_INC
#define rf_num_new_INC

#include "makros.h" // JT
//#include "FEMEnums.h"

#define NUM_FILE_EXTENSION ".num"
// C++ STL
//#include "prototyp.h"
#include <fstream>
#include <string>
#include <vector>

namespace ogs5
{

//----------------------------------------------------------------
class CNumerics
{
public:
    CNumerics(std::string);
    ~CNumerics(void);
    std::ios::pos_type Read(std::ifstream*);
public:
    // cf. Computational Geomachanics pp.62 WW
    double* DynamicDamping;
    /// For GMRES solver. 30.06.2010. WW
    long m_cols;
    // method
    std::string method_name;              //OK
    // PCS
    std::string pcs_type_name;
    // RENUMBER
    int renumber_method;
    int renumber_parameter;
    //
    // LS - Linear Solver
    int ls_method;
    int ls_max_iterations;
    int ls_error_method;
    double ls_error_tolerance;
    double ls_theta;
    int ls_precond;
    int ls_storage_method;
    std::string ls_extra_arg; //NW
    //
    // NLS - Non-linear Solver
    std::string nls_method_name;
    int nls_method;                       // Picard or Newton
    int nls_error_method;                 //WW
    int nls_max_iterations;
    double nls_relaxation;
    double nls_error_tolerance[DOF_NUMBER_MAX];        //JT2012: array function of dof
    double nls_plasticity_local_tolerance;
    //
    // CPL WW
    std::string cpl_variable;             // MB
    std::string cpl_process;              // JT
    std::string cpl_variable_JOD;        //JT->JOD. This one defaults to FLUX. I'm not sure what you want to do with it, but cpl_variable must default to "NONE".
    int cpl_max_iterations;
    int cpl_min_iterations;                  // JT2012
    double cpl_error_tolerance[DOF_NUMBER_MAX]; // JT2012: array function of dof
    bool cpl_error_specified;              // JT2012
    bool cpl_master_process;
    //
    // ELE
    int ele_gauss_points;                 // probably element-type-wise
    int ele_mass_lumping;
    int ele_upwind_method;                //CB
    double ele_upwinding;
    int ele_supg_method;                  //NW
    int ele_supg_method_length;           //NW
    int ele_supg_method_diffusivity;      //NW
    //FEM-FCT
    int fct_method;                       //NW
    unsigned int fct_prelimiter_type;     //NW
    double fct_const_alpha;               //NW
    // Deformation
    int GravityProfile;
    // LAGRANGE method //OK
    double lag_quality;
    int lag_max_steps;
    double lag_local_eps;
    int lag_time_weighting;
    double lag_min_weight;
    int lag_use_matrix;
    int lag_vel_method;
};

extern bool NUMRead(const std::string&, std::vector<CNumerics*> &num_vector);

}
#endif
