/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_st_new.h
 *
 * Created on 2004-01-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: Source Terms ST
   Task: class implementation
   Programing:
   01/2003 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_st_new_INC
#define rf_st_new_INC

#include <string>
#include <vector>
#include <iostream>

// FEM
#include "DistributionInfo.h" // TF
#include "ProcessInfo.h"                          // TF

namespace ogs5
{

class CSourceTerm : public ProcessInfo, public DistributionInfo
{
public:
    CSourceTerm();
    ~CSourceTerm();

    std::ios::pos_type Read(std::ifstream* in);
private:                                          // TF, KR
    void ReadDistributionType(std::ifstream* st_file);
    void ReadGeoType(std::ifstream* st_file);

public:
    std::string primaryvariable_name;
    std::string geo_name; // TF 05/2010
    std::string geo_type_name;

    int CurveIndex;
    std::vector<int> element_st_vector;

    double rill_height;
    double sorptivity, constant, rainfall, rainfall_duration, moistureDeficit /*1x*/;
    bool node_averaging, no_surface_water_pressure /*2x*/;

    bool channel;
    double channel_width;
    int geo_node_number;

    double* nodes;
    std::vector<int> node_number_vector;
    std::vector<double> node_value_vector;
    std::vector<int> node_renumber_vector;
    std::string tim_type_name;
    std::string interpolation_method;            //BG
    int TimeInterpolation;                        //BG

    std::string pcs_type_name_cond;
    std::string pcs_pv_name_cond;

    std::string fname;


    double geo_node_value;

    /**
     * is the source term coupled with another source term
     */
    bool _coupled;

    double normaldepth_slope;             // used only once in a global in rf_st_new

    /// Subdomain index for excvation simulation
    // 14.12.2010. WW
    int _sub_dom_idx;

    int fct_method;
    std::string fct_name;

    bool analytical;                      //2x?
    size_t number_of_terms;
    size_t _max_no_terms;                 // used only once in a global in rf_st_new
    size_t _no_an_sol;
    int analytical_material_group;        // used only once in a global in rf_st_new
    int resolution;                       // used only once in a global in rf_st_new
    double analytical_diffusion;          // used only once in a global in rf_st_new
    double analytical_porosity;           // used only once in a global in rf_st_new
    double analytical_tortousity;         // used only once in a global in rf_st_new
    double analytical_linear_sorption_Kd; // used only once in a global in rf_st_new
    double analytical_matrix_density;     // used only once in a global in rf_st_new
    double factor;

    std::string nodes_file;
    int msh_node_number;
    std::string msh_type_name;
    std::vector<int> PointsHaveDistribedBC;
    std::vector<double> DistribedBC;
    std::vector<double> node_value_vectorArea;
    std::vector<double*> normal2surface;
    std::vector<double*> pnt_parameter_vector;
    // 03.2010. WW
    long start_pos_in_st;
    double* GIS_shape_head;               // 07.06.2010. WW
    std::vector<double> precip_times;
    std::vector<std::string> precip_files;

    double _coup_leakance;
};



/**
 * read source term file
 * @param file_base_name base file name (without extension) containing the source terms
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if source terms found in file, else false
 */
bool STRead(const std::string& file_base_name,
         std::vector<CSourceTerm*> &st_vector);

#define ST_FILE_EXTENSION ".st"

}

#endif
