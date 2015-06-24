/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_tim_new.cpp
 *
 * Created on 2004-08-xx by Olaf Kolditz
 */

/**************************************************************************
   FEMLib - Object: TIM
   Task:
   Programing:
   08/2004 OK Implementation
   last modified:
**************************************************************************/
#include "rf_tim_new.h"
// C++ STL
//#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <cstdlib>
// FEM-Makros
#include "makros.h"
#include "Ogs5FileTools.h"

namespace ogs5
{

/**************************************************************************
   FEMLib-Method:
   Task: OBJ constructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::CTimeDiscretization(void)
    : Write_tim_discrete(false),tim_discrete(NULL) //YD
{
    step_current = 0;
    time_start = 0.0;
    time_end = 1.0;
    time_type_name = "CONSTANT";          //OK
    time_control_name = "NONE";           //kg44//JT
    time_unit = "SECOND";
    max_time_step = 1.e10;                //YD
    min_time_step = DBL_EPSILON;          //YD//JT Minimum allowed timestep, this process
    repeat = false;                       //OK/YD
    step_current = 0;                     //WW
    this_stepsize = 0.;                   //WW
    dt_sum = .0;                          //WW
    relative_error = 1.e-4;               //26.08.2008. WW
    absolute_error = 1.e-10;              //26.08.2008. WW
    h_max = 6;                            //27.08.2008. WW
    h_min = 0.2;                          //27.08.2008. WW
    hacc = 0.;                            //27.08.2008. WW
    erracc = 0.;                          //27.08.2008. WW
    PI_tsize_ctrl_type = -1;                 //27.08.2008. WW
    minimum_dt_reached = false;              //JT
    time_active = true;                      //JT
    time_independence = false;              //JT
    dt_failure_reduction_factor = 1.0;    //JT
    accepted_step_count = 0;              //JT
    rejected_step_count = 0;              //JT
    last_active_time = 0.0;                  //JT
    next_active_time = 0.0;                  //JT
    dynamic_time_buffer = 0;              //JT
    for(size_t ii=0; ii<DOF_NUMBER_MAX+1; ii++){
        dynamic_control_tolerance[ii] = -1.0;
    }
    dynamic_failure_threshold = .0;
    dynamic_minimum_suggestion = .0;
    iter_times = 0;
    last_dt_accepted = false;
    multiply_coef = .0;
    nonlinear_iteration_error = .0;
    recommended_time_step = 0;
    reject_factor = .0;
    rwpt_numsplits = 0;
    safty_coe = .0;
    time_control_manipulate = .0;
    time_current = .0;
    time_step_length = .0;
    time_step_length_neumann = .0;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ destructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::~CTimeDiscretization(void)
{
    if(tim_discrete)                      //YD
    {
        tim_discrete->close();
        if(tim_discrete)
            delete tim_discrete;
        time_step_vector.clear();
        time_adapt_tim_vector.clear();
        time_adapt_coe_vector.clear();
    }
}

std::ios::pos_type GetNextSubKeyword(std::ifstream* file,std::string* line, bool* keyword)
{
    char buffer[MAX_ZEILE];
    std::ios::pos_type position;
    position = file->tellg();
    *keyword = false;
    std::string line_complete;
    int i,j;
    // Look for next subkeyword
    while(!(line_complete.find("$") != std::string::npos) && (!file->eof()))
    {
        file->getline(buffer,MAX_ZEILE);
        line_complete = buffer;
        if(line_complete.find("#") != std::string::npos)
        {
            *keyword = true;
            return position;
        }
        //Anf���ngliche Leerzeichen ���berlesen, i=Position des ersten Nichtleerzeichens im string
        i = (int) line_complete.find_first_not_of(" ",0);
        j = (int) line_complete.find(";",i); //Nach Kommentarzeichen ; suchen. j = Position des Kommentarzeichens, j=-1 wenn es keines gibt.
        if(j < 0)
            j = (int)line_complete.length();
        //if(j!=i) break;                         //Wenn das erste nicht-leerzeichen ein Kommentarzeichen ist, zeile ���berlesen. Sonst ist das eine Datenzeile
        if(i != -1)
            *line = line_complete.substr(i,j - i);  //Ab erstem nicht-Leerzeichen bis Kommentarzeichen rauskopieren in neuen substring, falls Zeile nicht leer ist
    }
    return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   11/2004 OK string streaming by SB for lines
   10/2005 YD Time Controls
   08/2008 WW General classic time step size control (PI control)
**************************************************************************/
std::ios::pos_type CTimeDiscretization::Read(std::ifstream* tim_file)
{
    std::string sub_line;
    std::string line_string;
    std::string delimiter(" ");
    bool new_keyword = false;
    std::string hash("#");
    std::ios::pos_type position;
    std::string sub_string;
    bool new_subkeyword = false;
    std::string dollar("$");
    int no_time_steps = 0;
    double time_step_length;
    std::ios::pos_type position_subkeyword;
    std::stringstream line;
    std::string line_complete;
    int iter_times;                       //YD
    double multiply_coef;                 //YD
    int i;

    //========================================================================
    // Schleife ueber alle Phasen bzw. Komponenten
    while(!new_keyword)
    {
        if(new_subkeyword)
            tim_file->seekg(position,std::ios::beg);
        new_subkeyword = false;
        position = GetNextSubKeyword(tim_file,&line_string,&new_keyword);
        if(new_keyword)
            return position;
        /*
            position = tim_file->tellg();
            if(new_subkeyword)
              tim_file->seekg(position_subkeyword,ios::beg);
            new_subkeyword = false;
            tim_file->getline(buffer,MAX_ZEILE);
            line_string = buffer;
           if(line_string.size()<1) // empty line
              continue;
            if(Keyword(line_string))
              return position;
         */
        //....................................................................

        // subkeyword found
        if(line_string.find("$PCS_TYPE") != std::string::npos)
        {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> pcs_type_name;
            line.clear();
            //m_pcs = PCSGet(pcs_type_name); // kg44 inserted to overwrite default Richards_flow
            // this works only of pcs_type is read before adaption
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$TIME_START") != std::string::npos)
        {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> time_start;
            line.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$TIME_END") != std::string::npos)
        {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> time_end;
            line.clear();
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$TIME_UNIT") != std::string::npos)
        {
            *tim_file >> time_unit >> std::ws; //WW unit of time
            continue;
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$INDEPENDENT") != std::string::npos) // JT2012
        {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> time_independence;
            // =0: Process will adopt minimum time step of all processes.
            // =1: Process will have it's own time step (may not execute on every time step). It will be independent.
            line.clear();
            continue;
        }
        //....................................................................
        /* //WW
           if(line_string.find("$TIME_FIXED_POINTS")!=string::npos) { // subkeyword found
           int no_fixed_points;
           double fixed_point;
           line.str(readNonBlankLineFromInputStream(*tim_file));
           line >> no_fixed_points;
           line.clear();
           for(i=0;i<no_fixed_points;i++) {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> fixed_point;
           fixed_point_vector.push_back(fixed_point);
           line.clear();
           }
           continue;
           }
         */
        //....................................................................
        // subkeyword found
        if(line_string.find("$TIME_STEPS") != std::string::npos)
            while((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
            {
                position = tim_file->tellg();
                line_string = readNonBlankLineFromInputStream(*tim_file);
                if(line_string.find("#") != std::string::npos)
                    return position;
                if(line_string.find("$") != std::string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                line.str(line_string);
                line >> no_time_steps;
                line >> time_step_length;
                for(i = 0; i < no_time_steps; i++)
                    time_step_vector.push_back(time_step_length);
                line.clear();
            }
        // subkeyword found
        if(line_string.find("$TIME_SPLITS") != std::string::npos)
        {
            line.str(readNonBlankLineFromInputStream(*tim_file));
            line >> rwpt_numsplits;
            line.clear();
            continue;
        }
        // 25.08.2008. WW
        if(line_string.find("$CRITICAL_TIME") != std::string::npos)
            while((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
            {
                position = tim_file->tellg();
                line_string = readNonBlankLineFromInputStream(*tim_file);
                if(line_string.find("#") != std::string::npos)
                    return position;
                if(line_string.find("$") != std::string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                line.str(line_string);
                double crtime;
                line >> crtime;
                critical_time.push_back(crtime);
                line.clear();
            }
        // subkeyword found
        if(line_string.find("$TIME_CONTROL") != std::string::npos)
            while((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
            {
                position = tim_file->tellg();
                line_string = readNonBlankLineFromInputStream(*tim_file);

                if(line_string.find("#") != std::string::npos)
                    return position;
                if(line_string.find("$") != std::string::npos)
                {
                    new_subkeyword = true;
                    break;
                }
                line.str(line_string);
                line >> time_control_name;
                line.clear();

                if(time_control_name == "PI_AUTO_STEP_SIZE") // 26.08.2008. WW
                {
                    line.str(readNonBlankLineFromInputStream(*tim_file));
                    line >> PI_tsize_ctrl_type >> relative_error >>
                    absolute_error >> this_stepsize;
                    //13.03.2008. WW
                    int real_type = (int)(PI_tsize_ctrl_type / 10);
                    if(real_type < 10 && real_type > 0) //
                    {
                        PI_tsize_ctrl_type = real_type;
                        line >> h_min >> h_max >> max_time_step;
                    }
                    else
                        max_time_step = 0.0;
                    line.clear();
                }
                else if(time_control_name.find("DYNAMIC_VARIABLE") != std::string::npos) // JT2012
                {
                    // DYNAMIC TIME STEP SERIES
                    line.str(readNonBlankLineFromInputStream(*tim_file));
                    int num_tolerances = DOF_NUMBER_MAX;
                    double include_third_variable = -1.0;
                    double third_variable_tolerance = -1.0;
                    //
                    line >> dynamic_control_error_method;                        // Corresponds to FiniteElement::ErrorMethod. Defines how tolerance is applied.
                    line >> time_step_length;                                    // initial_dt
                    line >> min_time_step;                                        // min_dt
                    line >> dynamic_failure_threshold;                            // threshold to force a time failure (recommend 1.1-2.0. If tolerance is exceeeded on first iteration by this factor, dt will be decreased and a failure forced)
                    line >> dynamic_control_tolerance[DOF_NUMBER_MAX];            // max_increase_factor (dt never allowed to increase faster than this factor (i.e. 1.5)
                    //
                    // tolerances[1:dof]: One tolerance for each degree of freedom in process (or only one tolerance for certain methods)
                    //switch(FiniteElement::convertErrorMethod(dynamic_control_error_method))
                    //{
                    //    case FiniteElement::ENORM: // only 1 tolerance required
                    //        num_tolerances = 1;
                    //        line >> dynamic_control_tolerance[0];
                    //        //
                    //        // Next entry is OPTIONAL (may be left blank)!
                    //        line >> include_third_variable; // if >0.0, include a 3rd variable in the global norm (PwSnw model: Pc -- PcPnw model: Snw -- PwPnw model: Snw)
                    //        break;
                    //    //
                    //    case FiniteElement::ERNORM: // only 1 tolerance required
                    //        num_tolerances = 1;
                    //        line >> dynamic_control_tolerance[0];
                    //        //
                    //        // Next entry is OPTIONAL (may be left blank)!
                    //        line >> include_third_variable; // if >0.0, include a 3rd variable in the global norm (PwSnw model: Pc -- PcPnw model: Snw -- PwPnw model: Snw)
                    //        break;
                    //    //
                    //    case FiniteElement::EVNORM: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
                    //        for(int i=0; i<num_tolerances; i++)
                    //            line >> dynamic_control_tolerance[i];
                    //        //
                    //        // Next entry is OPTIONAL (may be left blank)!
                    //        line >> third_variable_tolerance; // tolerance of a third variable (PwSnw model: Pc -- PcPnw model: Snw -- PwPnw model: Snw)
                    //        break;
                    //    //
                    //    case FiniteElement::LMAX: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance required. Applies to x,y,z)
                    //        for(int i=0; i<num_tolerances; i++)
                    //            line >> dynamic_control_tolerance[i];
                    //        //
                    //        // Next entry is OPTIONAL (may be left blank)!
                    //        line >> third_variable_tolerance; // tolerance of a third variable (PwSnw model: Pc -- PcPnw model: Snw -- PwPnw model: Snw)
                    //        break;
                    //    //
                    //    case FiniteElement::BNORM:
                    //        ScreenMessage("ERROR in TIMRead. BNORM not configured for time control.\n");
                    //        ScreenMessage("We suggest ENORM as a valid companion for NEWTON couplings.\n");
                    //        exit(1);
                    //        break;
                    //    //
                    //    default:
                    //        ScreenMessage("ERROR in TIMRead. Invalid error method selected for dynamic time control.\n");
                    //        exit(1);
                    //        break;
                    //}
                    ////
                    if(num_tolerances > 1)
                        dynamic_control_tolerance[DOF_NUMBER_MAX] = third_variable_tolerance;
                    else
                        dynamic_control_tolerance[DOF_NUMBER_MAX] = include_third_variable;
                    line.clear();
                }
                else if(time_control_name.find("DYNAMIC_COURANT") != std::string::npos) // JT2012
                {
                    std::cout<<"Josh apologizes (especially to Marc), but promises DYNAMIC_COURANT will be merged into the next release."<<std::endl;
                    std::cout<<"emoticon:: sad face"<<std::endl;
                    exit(1);
                    //
                    // DYNAMIC TIME STEP SERIES
                    line.str(readNonBlankLineFromInputStream(*tim_file));
                    //
                    line >> time_step_length;                                    // initial_dt
                    line >> min_time_step;                                        // min_dt
                    line >> dynamic_control_tolerance[DOF_NUMBER_MAX];            // max_increase_factor (dt never allowed to increase faster than this factor (i.e. 1.5)
                    line >> dynamic_control_tolerance[0];                        // desired courant number
                    //
                    // ADDITIONAL OPTIONS TO RESTRICT CALCULATION TO CERTAIN ELEMENTS
                    if(time_control_name.find("CONCENTRATION") != std::string::npos){    // DYNAMIC_COURANT_CONCENTRATION
                        line >> dynamic_control_tolerance[1]; // Concentration threshold (elements with concentration beneath this value are not included in Courant restriction)
                    }
                    else if(time_control_name.find("TEMPERATURE") != std::string::npos){// DYNAMIC_COURANT_TEMPERATURE
                        line >> dynamic_control_tolerance[1]; // Temperature threshold (elements with temperature BENEATH this value are not included in Courant restriction)
                        line >> dynamic_control_tolerance[2]; // Temperature threshold (elements with temperature ABOVE   this value are not included in Courant restriction)
                    }
                    //
                    line.clear();
                }
                else if(time_control_name.find("DYNAMIC_PRESSURE") != std::string::npos) // JT2012
                {
                    // DYNAMIC TIME STEP SERIES
                    line.str(readNonBlankLineFromInputStream(*tim_file));
                    //
                    line >> time_step_length;                                    // initial_dt
                    line >> min_time_step;                                        // min_dt
                    line >> dynamic_control_tolerance[DOF_NUMBER_MAX];            // max_increase_factor (dt never allowed to increase faster than this factor (i.e. 1.5)
                    line >> dynamic_control_tolerance[0];                        // pressure tolerance (mean pressure)
                    //
                    line.clear();
                }
                else if(time_control_name == "STEP_SIZE_RESTRICTION") // 26.08.2008. WW
                {
                    line.str(readNonBlankLineFromInputStream(*tim_file));
                    line >> h_min >> h_max;
                    line.clear();
                }
                else if(time_control_name == "NEUMANN"){
                    line.clear();
                }
                else if(time_control_name == "ERROR_CONTROL_ADAPTIVE")
                {
                    //m_pcs->adaption = true;
                    line.clear();
                }
                else if(time_control_name == "SELF_ADAPTIVE")
                {
                    //m_pcs->adaption = true; JOD removed
                    //WW minish = 10;
                    while((!new_keyword) || (!new_subkeyword) ||
                          (!tim_file->eof()))
                    {
                        position = tim_file->tellg();
                        line_string = readNonBlankLineFromInputStream(*tim_file);
                        if(line_string.find("#") != std::string::npos)
                            return position;
                        if(line_string.find("$") != std::string::npos)
                        {
                            new_subkeyword = true;
                            break;
                        }
                        if(line_string.find("MAX_TIME_STEP") !=
                           std::string::npos)
                        {
                            *tim_file >> line_string;
                            max_time_step = strtod(
                                    line_string.data(),NULL);
                            line.clear();
                            // kg44 should not break break;
                        }
                        if(line_string.find("MIN_TIME_STEP") !=
                           std::string::npos)
                        {
                            *tim_file >> line_string;
                            min_time_step = strtod(
                                    line_string.data(),NULL);
                            line.clear();
                            // kg44 should not break break;
                        }
                        /*  //WW
                           if(line_string.find("MINISH")!=string::npos){
                           *tim_file >> line_string;
                           minish = strtod(line_string.data(),NULL);
                           line.clear();
                           }
                         */
                        if(line_string.find("M") == std::string::npos)
                        {
                            line.str(line_string);
                            line >> iter_times;
                            line >> multiply_coef;
                            time_adapt_tim_vector.push_back(iter_times);
                            time_adapt_coe_vector.push_back(
                                    multiply_coef);
                            line.clear();
                        }
                    } // end of while loop adaptive
                // end of if "SELF_ADAPTIVE"
                }
                else{
                    //ScreenMessage("ERROR: Unrecognized time control type.\n");
                    exit(1);
                }
            }             // end of while
        // end of "TIME_CONTROL"
        //....................................................................
        /* //WW
           if(line_string.find("$SUBSTEPS")!=string::npos) { // subkeyword found JOD 4.7.10
           *tim_file>>sub_steps>>ws;
           continue;
           }
         */
        //....................................................................
    }                                     // end of while(!new_keyword)

    return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool TIMRead(const std::string& file_base_name, std::vector<CTimeDiscretization*> &time_vector)
{
    //----------------------------------------------------------------------
    //OK  TIMDelete();
    //----------------------------------------------------------------------
    CTimeDiscretization* m_tim = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    std::string tim_file_name = file_base_name + TIM_FILE_EXTENSION;
    std::ifstream tim_file (tim_file_name.data(),std::ios::in);
    if (!tim_file.good())
        return false;
    tim_file.seekg(0L,std::ios::beg);
    //========================================================================
    // Keyword loop
    std::cout << "TIMRead ... " << std::flush;
    while (!tim_file.eof())
    {
        tim_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos) {
            std::cout << "done, read " << time_vector.size() << " time stepping properties" <<
            std::endl;
           return true;
        }
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#TIME_STEPPING") != std::string::npos)
        {
            m_tim = new CTimeDiscretization();
            position = m_tim->Read(&tim_file);
            m_tim->time_current = m_tim->time_start;
            //----------------------------------------------------------------------
            time_vector.push_back(m_tim);
            tim_file.seekg(position,std::ios::beg);
        }                         // keyword found
    }                                     // eof
    return true;
}

}
