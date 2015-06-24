/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rfmat_cp.cpp
 *
 * Created on 2004-10-xx by Sebastian Bauer
 */

/**************************************************************************/
/* ROCKFLOW - Modul: rfmat_cp.c
 */
/* Task:
   Methods for ComponentProperties
 */
/* Programming:
   10/2004   SB  First Implemented
 */
/**************************************************************************/
#include "rfmat_cp.h"
// C
#include <math.h>
// C++
#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>

#include "makros.h"
#include "Ogs5FileTools.h"

using namespace std;

namespace ogs5
{

/*========================================================================*/
/* Component Properties                                                  */
/*========================================================================*/

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties Constructor
   Programing:
   02/2004 SB Implementation
**************************************************************************/
CompProperties::CompProperties(/* int n // HS we do not need this. */)
    : idx(std::numeric_limits<size_t>::max())
{
    mobil = 1;  // by default, set to mobile species.
    transport_phase = 0; // by default, set to the 1st phase.
    fluid_phase = 0; // by default, set to water

    diffusion_model = -1;
    count_of_diffusion_model_values = 0;

    decay_model = -1;
    count_of_decay_model_values = 0;

    isotherm_model = -1;
    count_of_isotherm_model_values = 0;

    bubble_velocity_model = -1;
    bubble_velocity[0] = bubble_velocity[1] = bubble_velocity[2] = 0.0;

    OutputMassOfComponentInModel = 0;                        // 05/2012 BG
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties Destructor
   Programing:
   02/2004 SB Implementation
**************************************************************************/
CompProperties::~CompProperties(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties read function
   Programing:
   02/2004 SB Implementation - adapted from OK rf_bc_new
   10/2004 SB Adapted to new file structure
   01/2005 OK boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool CPRead(const std::string &file_base_name, std::vector<CompProperties*> &cp_vec)
{
    //----------------------------------------------------------------------
    //OK  MCPDelete();
    //----------------------------------------------------------------------
    CompProperties* m_cp = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    std::string cp_file_name = file_base_name + CP_FILE_EXTENSION;
    ifstream cp_file (cp_file_name.data(),ios::in);
    if (!cp_file.good())
        return false;
    cp_file.seekg(0L,ios::beg);

    //========================================================================
    cp_vec.clear();
    cout << "CPRead" << endl;
    // Schleife ueber alle Phasen bzw. Komponenten
    while (!cp_file.eof())
    {
        cp_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos)
            break;
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#COMPONENT_PROPERTIES") != std::string::npos)
        {
            m_cp = new CompProperties();
            position = m_cp->Read(&cp_file);
            // HS the index of this component is filled
            // one after another by its sequence in mcp file.
            m_cp->idx = cp_vec.size();
            cp_vec.push_back(m_cp);
            m_cp = NULL;
            cp_file.seekg(position,ios::beg);
        }                         // keyword found
    }                                     // eof
//    // immediately check if enough PCS objects are available
//    size_t pcs_mt_count = 0;
//    size_t pcs_rwpt_count = 0;
//    size_t i;
//    for ( i = 0; i < pcs_vector.size(); i++ )
//    {
//        if ( pcs_vector[i]->getProcessType() == FiniteElement::MASS_TRANSPORT )
//            pcs_mt_count++;
//        if ( pcs_vector[i]->getProcessType() == FiniteElement::RANDOM_WALK )
//            pcs_rwpt_count++;
//    }
//    if ( pcs_rwpt_count == 0) // HS, no random walk detected.
//    {
//        if ( pcs_mt_count != cp_vec.size() || pcs_mt_count != cp_name_2_idx.size() )
//        {
//            DisplayMsgLn(
//                    "Mass transport components and Mass transport processes do not fit!");
//            exit(1);
//        }
//        else
//        {
//            // and then link MCP with the PCS.
//            std::map <int, CompProperties*>::iterator cp_iter = cp_vec.begin();
//            for ( i = 0; i < pcs_vector.size(); i++ )
//                if ( pcs_vector[i]->getProcessType() == FiniteElement::MASS_TRANSPORT )
//                {
//                    cp_iter->second->setProcess( pcs_vector[i] );
//                    ++cp_iter;
//                }
//        } // end of else
//    }
    return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties read function
   Programing:
   02/2004 SB Implementation - adapted from OK rf_bc_new
   05/2007 PCH: Anisotropic diffusion coefficient added
**************************************************************************/
ios::pos_type CompProperties::Read(ifstream* rfd_file)
{
    //  char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::string delimiter(" ");
    bool new_keyword = false;
    std::string hash("#");
    std::stringstream in;
    int j;
    //  double *read_help;
    //  long index;
    std::ios::pos_type position;

    //========================================================================
    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        //    new_subkeyword = false;
        position = rfd_file->tellg();

        //    if(!GetLineFromFile(line,rfd_file)) break;
        //    line_string = line;
        line_string = readNonBlankLineFromInputStream(*rfd_file);
        if(line_string.size() < 1)
            break;

        if(line_string.find(hash) != std::string::npos)
        {
            new_keyword = true;
            break;
        }

        /* Keywords nacheinander durchsuchen */
        //....................................................................
        // subkeyword found
        if(line_string.find("$NAME") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> compname; //sub_line
            in.clear();
            //      compname = (char *) sub_line.c_str();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$MOBILE") != std::string::npos)
        {
            //      rfd_file->getline(line,MAX_ZEILE);
            //      line_string = line;
            //      mobil = atoi(line_string.substr(0).data());
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> mobil;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$TRANSPORT_PHASE") != std::string::npos)
        {
            //      rfd_file->getline(line,MAX_ZEILE);
            //      line_string = line;
            //      transport_phase = atoi(line_string.substr(0).data());
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> transport_phase;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$FLUID_PHASE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> fluid_phase;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$MOL_MASS") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> molar_mass;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$CRITICAL_PRESSURE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> pc;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$CRITICAL_TEMPERATURE") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> Tc;
            in.clear();
        }
        //....................................................................
        if(line_string.find("$ACENTRIC_FACTOR") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> omega;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$FLUID_ID") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> fluid_id;
            in.clear();
        }
        if(line_string.find("$MOLAR_VOLUME") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> Vm;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$VOLUME_DIFFUSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> Vd;
            in.clear();
        }
        //....................................................................
        // subkeyword found
        if(line_string.find("$DIFFUSION") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> diffusion_model;
            if(diffusion_model == 0) //read curve

                in >> diffusion_function_name;
            else
            {
                count_of_diffusion_model_values =
                        GetNumberDiffusionValuesCompProperties(diffusion_model);
//                // unknown parameter model
//                if(count_of_diffusion_model_values < 0)
//                {
//                    DisplayMsgLn(" Unknown Diffusion model - program stops !");
//                    exit(1);
//                }

                if (diffusion_model > 0)
                    //        read_help = (double *) Malloc(count_of_diffusion_model_values * sizeof(double));

                    for (j = 0; j < count_of_diffusion_model_values; j++)
                    {
                        if(in.peek() > 0)

                            /*
                                        if(StrTestInv(&sub[p_sub += pos],&pos)) {
                                                StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
                                                sprintf(group_name,"%s%d",name_component_properties,j);
                                                ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
                                                SetInverseVariable(name,group_name,1,(void *)ptr_d);
                                                name=(char *)Free(name);
                                            }
                                            else  */
                            //                    in >> read_help[j];
                            in >> diffusion_model_values[j];
                        else
                            cout <<
                            "Warning: Missing diffusion model values for component "
                                 << compname << endl;
                        if((diffusion_model == 1) && (j == 0))
                            if((diffusion_model_values[j] < 0.0) ||
                               (diffusion_model_values[j] > 100000.0))
                                cout <<
                                "Warning: Funny diffusion model values specified"
                                     << endl;
                    } //end for(j...)

                //        diffusion_model_values = read_help;

//                if (diffusion_model < 0)
//                    DisplayMsgLn(
//                            "Error: Diffusion model must be larger than or 0");
            }
            in.clear();
        }                         // subkeyword found
                                  // subkeyword found
        if(line_string.find("$DECAY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> decay_model;
            if(decay_model == 0) //curve

                in >> decay_function_name;
            else
            {
                count_of_decay_model_values = GetNumberDecayValuesCompProperties(
                        decay_model);

//                if(count_of_decay_model_values < 0) // unknown parameter model
//                {
//                    DisplayMsgLn(
//                            " Unknown Aqueous Decay model - program stops !");
//                    exit(1);
//                }

                if (decay_model > 0)
                    //        read_help = (double *) Malloc(count_of_decay_model_values * sizeof(double));

                    for (j = 0; j < count_of_decay_model_values; j++)
                        /*
                                    if(StrTestInv(&sub[p_sub += pos],&pos)) {
                                            StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
                                            sprintf(group_name,"%s%d",name_component_properties,j);
                                            ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
                                            SetInverseVariable(name,group_name,1,(void *)ptr_d);
                                            name=(char *)Free(name);
                                        }
                                        else  */
                        //                    in >> read_help[j];
                        in >> decay_model_values[j];

                //end for(j...)

                //        decay_model_values = read_help;

//                if (decay_model < -1)
//                    DisplayMsgLn(
//                            "Error: Aqueous decay model must be larger than or 0");
            }
            in.clear();
        }                         // subkeyword found
                                  // subkeyword found
        if(line_string.find("$ISOTHERM") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> isotherm_model;
            if(isotherm_model == 0)
                in >> isotherm_function_name;
            else
            {
                count_of_isotherm_model_values =
                        GetNumberIsothermValuesCompProperties(isotherm_model);

//                if(count_of_isotherm_model_values < 0) // unknown parameter model
//                {
//                    DisplayMsgLn(" Unknown Isotherm model - program stops !");
//                    exit(1);
//                }

                if (isotherm_model > 0)
                    //        read_help = (double *) Malloc(count_of_isotherm_model_values * sizeof(double));

                    for (j = 0; j < count_of_isotherm_model_values; j++)
                        /*
                                    if(StrTestInv(&sub[p_sub += pos],&pos)) {
                                            StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
                                            sprintf(group_name,"%s%d",name_component_properties,j);
                                            ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
                                            SetInverseVariable(name,group_name,1,(void *)ptr_d);
                                            name=(char *)Free(name);
                                        }
                                        else  */
                        //                    in >> read_help[j];
                        in >> isotherm_model_values[j];

                //end for(j...)

                //        isotherm_model_values = read_help;
//                if (isotherm_model < 0)
//                    DisplayMsgLn(
//                            "Error: Isotherm model must be larger than or 0");
            }
            in.clear();
        }                         // subkeyword found

        // subkeyword found
        if(line_string.find("$BUBBLE_VELOCITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> bubble_velocity_model;
            if(bubble_velocity_model == 1)
                in >> bubble_velocity[0] >> bubble_velocity[1] >>
                bubble_velocity[2];
            in.clear();
        }                         // subkeyword found

        //....................................................................

        // parameters for NAPL dissolution CB140708
        // subkeyword found
        if(line_string.find("$MOLAR_DENSITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> molar_density;
            in.clear();
            if(molar_density <= 0) // unphysical entry
            {
//                DisplayMsgLn(
//                        "Error in MOLAR_DENSITY - setting molar_density to 1.0!");
                molar_density = 1.0;
            }
        }
        // subkeyword found
        if(line_string.find("$MOLAR_WEIGHT") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> molar_weight;
            in.clear();
            if(molar_weight <= 0) // unphysical entry
            {
//                DisplayMsgLn("Error in MOLAR_WEIGHT - setting molar_weight to 1.0!");
                molar_weight = 1.0;
            }
        }
        // subkeyword found
        if(line_string.find("$MAXIMUM_AQUEOUS_SOLUBILITY") != std::string::npos)
        {
            in.str(readNonBlankLineFromInputStream(*rfd_file));
            in >> max_solubility;
            in.clear();
            if(max_solubility <= 0) // unphysical entry
            {
//                DisplayMsgLn(
//                        ": Error in MAXIMUM_AQUEOUS_SOLUBILITY - setting max_solubility to 1.0!");
                max_solubility = 1.0;
            }
        }
        //....................................................................
        // subkeyword found
        // Output of the Mass of the component within the model if subkeyword is found, BG 05/2012
        if(line_string.find("$OutputMassOfComponentInModel") != std::string::npos)
        {
            OutputMassOfComponentInModel = 1;
        }
    }                                     //end while
    return position;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberDiffusionValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit des Diffusionsmodells die Anzahl der
   benoetigten Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int diffusion_model : Diffusionsmodell
 */
/* Ergebnis:
   - int - Anzahl der Koeffizienten
 */
/* Programmaenderungen:
   06/2000    AH      Erste Version  (Fall 0 bis 7)
   02/2004    SB      Adapted to CompProperties
 */
/**************************************************************************/
int CompProperties::GetNumberDiffusionValuesCompProperties(int diffusion_model)
{
    int n = -1;

    if (diffusion_model < -1)
        return n;

    switch (diffusion_model)
    {
    case -1:
        n = 0;
        break;                    /* Keine Diffusion */
    case 0:                               /* curve */
        n = 1;
    case 1:
        n = 1;
        break;                    /* Konstanter Diffusionswert */
    case 2:
        n = 1;
        break;                    /* Variabler Diffusionswert (Zeitabhaengig) */

    /* Diffusionskoeffizienten in Wasser Daq [m^2/s] */
    case 3:
        n = 1;
        break;                    /* Worch, 1993 */
    case 4:
        n = 1;
        break;                    /* Hayduk und Laudie, 1974 */
    case 5:
        n = 2;
        break;                    /* Wilke und Chang, 1955 */
    case 6:
        n = 1;
        break;                    /* Stokes-Einstein (fuer Partikel/Makromolekuele) */

    /* Diffusionskoeffizienten in Luft Dg [m^2/s] */
    case 7:
        n = 2;
        break;                    /* FSG-Methode, Lyman et al., 1990 */
    case 8:
        n = 2;
        break;                    /* Archies Law */
    case 9:
        n = 2;
        break;                    /* Archies Law */
    }                                     /* switch */

    /* switch */

    return n;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberDecayValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit des in der sorbierten Phase Zerfallsmodells
   die Anzahl der benoetigten Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int isotherm : Isotherme-Identifikator
 */
/* Ergebnis:
   - int - Anzahl der Isothermen-Koeffizienten
 */
/* Programmaenderungen:
   06/2000    AH      Erste Version (Fall 0 bis 7)
   09/2004      SB      Adapted to CompProperties
 */
/**************************************************************************/
int CompProperties::GetNumberDecayValuesCompProperties(int decay_model)
{
    int n = -1;

    if (decay_model < -1)
        return n;

    switch (decay_model)
    {
    case -1: n = 0;
        break;
    case  0: n = 1;
        break;                    /* No Decay */
    case  1: n = 2;
        break;                    /* Any Order Decay with constant Decay Rate */
    case  2: n = 2;
        break;                    /* Monod or Michaelis-Menten kinetics with constant rate coefficients */
    default:
        //DisplayMsgLn(" Error: Unknown model for decay  ");
        break;
    }                                     /* switch */

    return n;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberIsothermValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit der Isotherme die Anzahl der Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int isotherm : Isotherme-Identifikator
 */
/* Ergebnis:
   - int - Anzahl der Isothermen-Koeffizienten
 */
/* Programmaenderungen:
   10/1999     AH         Erste Version  (Fall 0 bis 11)
   02/2004     SB         Adapted to CompProperties - nur method 0-3+15
   04/2006     CMCD       Included for analytical fracture surface method, Tang et al. 1981, WRR.
 */
/**************************************************************************/
int CompProperties::GetNumberIsothermValuesCompProperties(int isotherm)
{
    int n = 0;

    if (isotherm < -1)
        return -isotherm;

    switch (isotherm)
    {
    case -1: n = 0;                       /* no isotherm */
    case 0:
        n = 1;
        break;                    /* get KD from curve derivative */
    case 1:
        n = 1;
        break;                    /* Henry Isotherm */
    case 2:
        n = 2;
        break;                    /* Freundlich Isotherm */
    case 3:
        n = 2;
        break;                    /* Langmuir Isotherm */
    case 4:
        n = 1;
        break;                    /* Linear Isotherm, for fracture surface CMCD */
    case 5:
        n = 3;
        break;                    /* two-rate model */
    /*
        case 4:
            n = 3;  break;                  // Freundlich Langmuir Isotherm
        case 5:
            n = 4;  break;                  // Double Langmuir Isotherm
        case 6:
            n = 3;  break;                  // Extented Freundlich Isotherm
        case 7:
            n = 3;  break;                  // Gunary Isotherm
        case 8:
            n = 3;  break;                  // Fitter-Sutton Isotherm
       case 9:
       n = 4;  break;                  // Barry Isotherm
       case 10:
       n = 2;  break;                  // Power Isotherm
       case 11:
       n = 2;  break;                  // Modified Kielland Isotherm
     */
    /*    case 15:
          n=1;    break;
     */
    default:
        //DisplayMsgLn(" Error - this ISOTHERM model found ");
        break;
    }                                     /* switch */

    return n;
}
}
