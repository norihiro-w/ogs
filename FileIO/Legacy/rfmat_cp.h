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
/* ROCKFLOW - Modul: rfmat_cp.h
 */
/* Task:
   Methods for ComponentProperties
 */
/* Programming:
   10/2004   SB  First Implemented
 */
/**************************************************************************/

/* Schutz gegen mehrfaches Einfuegen */
#ifndef rfmat_cp_INC
#define rfmat_cp_INC

#include <fstream>
#include <map>
#include <string>
#include <vector>

#define CP_FILE_EXTENSION ".mcp"                  /* File extension for component properties input file */

namespace ogs5
{

/*************************************************************************

   Class ComponentProperties

 **************************************************************************/

class CompProperties// : public ProcessInfo
{
private:
public:
    std::ios::pos_type Read(std::ifstream*);        /* Lesefunktion f? eine Instanz von CompProperties */

    /* constructor */
    CompProperties();
    /* destructor */
    ~CompProperties(void);

    size_t idx;                             /* the unique index of this component. not effective, saved for future*/
    std::string compname;                     /* component name */
    int mobil;                           /* flag mobil */
    int transport_phase;                 /* number of phase, in which component is transported */
    int fluid_phase;
    double molar_mass;
    double pc;                        //Critical pressure [Pa]
    double omega;                // accentric factor [-]
    double Tc;                    // Critical temperature [K]
    int fluid_id;            // Requred to detect a particular fluid from *.mcp files: 0 for CO2; 1 for H2O; 2 for CH4; 3 for N2
    double Vm;                // Vm [m3/kmol]: Molar volume used in interation parameter calculation
    double Vd;                 //[cm3/mol] diffusion volume

    /* Diffusionsmodelle und zugehoerige Beschreibungswerte */
    int diffusion_model;                  /* Zerfallsmodell in geloester Phase */
    int count_of_diffusion_model_values;  /* Anzahl der Parameter zur Spezifikation des Diffusionsmodells */
    double diffusion_model_values[10];    /* Parameter fuer das Diffusionsmodell */
    int diffusion_function_name;
    /* Zugriff auf Number of Parameters */
    /* Zerfallsmodelle und zugehoerige Beschreibungswerte in der geloesten Phase */
    int decay_model;                      /* Zerfallsmodell in geloester Phase */
    int count_of_decay_model_values;      /* Anzahl und Werte zur Spezifikation der */
    double decay_model_values[10];        /* Parameter fuer Zerfallsprozess wie z.B. Zerfallsrate */
    int decay_function_name;
    /* Sorption */
    int isotherm_model;                   /* Isothermen-Typ */
    int count_of_isotherm_model_values;   /* Anzahl der Isothermen-Koeffizienten */
    double isotherm_model_values[10];     /* Isothermen-Koeffizienten */
    int isotherm_function_name;
    /* bubble velocity */
    int bubble_velocity_model;
    double bubble_velocity[3];            /* velocity of rising bubbles */

    /* parameters for NAPL dissolution CB140708 */
    double molar_density;
    double molar_weight;
    double max_solubility;

    int OutputMassOfComponentInModel;


private:
    int GetNumberDiffusionValuesCompProperties(int );
    int GetNumberDecayValuesCompProperties(int); /* Zugriff auf Number of Parameters */
    int GetNumberIsothermValuesCompProperties(int);

};

extern bool CPRead(const std::string &file_base_name, std::vector<CompProperties*> &cp_vec);

}

#endif
