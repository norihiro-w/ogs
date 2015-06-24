/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_tim_new.h
 *
 * Created on 2004-06-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: TIM
   Task: class implementation
   Programing:
   06/2004 OK Implementation
   12/2008 WW PI time step control
   last modified:
**************************************************************************/
#ifndef rf_tim_new_INC
#define rf_tim_new_INC
// C++ STL
//#include <fstream>
#include <iostream>
#include <vector>
//#include <string>
//#include <vector>
#include "makros.h" // JT2012

namespace ogs5
{

//----------------------------------------------------------------
class CTimeDiscretization
{
public:
    CTimeDiscretization(void);
    ~CTimeDiscretization(void);
    std::ios::pos_type Read(std::ifstream*);

public:
    double time_step_length_neumann;      //YD
    double time_step_length;              //YD
    double safty_coe;
    double dt_sum;                        // 17.09.2007 WW
    double this_stepsize;
    double relative_error;
    double absolute_error;
    double reject_factor; //for automatic timestep control: BG
    double h_min;
    double h_max;
    double hacc;
    double erracc;
    double dynamic_control_tolerance[DOF_NUMBER_MAX+1];    //JT2012
    std::string dynamic_control_error_method;            //JT2012
    int dynamic_time_buffer;                            //JT2012
    double dynamic_minimum_suggestion;                    //JT2012
    double dynamic_failure_threshold;                    //JT2012
    int PI_tsize_ctrl_type;
    std::vector<double> critical_time;
    std::string file_base_name;
    // TIM
    std::vector<double>time_step_vector;
    std::vector<int> time_adapt_tim_vector;
    std::vector<double>time_adapt_coe_vector;
    //WW vector<double>fixed_point_vector;

    //WW vector<double> time_step_target_vector; // kg44 for adaptive steps..intermediate time target that need to be reached
    double time_start;
    double time_end;
    double time_current;
    double time_control_manipulate;       //CMCD
    double next_active_time;
    double last_active_time;
    double recommended_time_step;
    double dt_failure_reduction_factor;
    int step_current;
    bool repeat;                          //OK/YD
    bool time_active;                    //JT2012
    bool time_independence;                //JT2012
    bool last_dt_accepted;                //JT2012
    bool minimum_dt_reached;            //JT2012
    long accepted_step_count;            //JT2012
    long rejected_step_count;            //JT2012
    //
    // PCS
    std::string pcs_type_name;            //OK
    // NUM
    std::string time_type_name;           //OK
    std::string time_control_name;
    std::string time_unit;                //WW
    double iter_times;                    //YD
    double multiply_coef;                 //YD
    double max_time_step;                 //YD
    double min_time_step;
    //
    //WW double minish; // JOD
    //WW int sub_steps; // JOD 4.7.10
    bool Write_tim_discrete;              //YD
    std::fstream* tim_discrete;           //YD
    double nonlinear_iteration_error;     //OK/YD
    //WW double max_adaptive_factor; // kg44
    //WW double max_adaptive_concentration_change; // kg44

    int rwpt_numsplits;
};

extern bool TIMRead(const std::string&, std::vector<CTimeDiscretization*> &time_vector);
#define TIM_FILE_EXTENSION ".tim"

}

#endif
