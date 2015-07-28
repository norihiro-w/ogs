/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_cur.cpp
 *
 * Created on 2007-04-xx by Olaf Kolditz
 *
 */

#include "rf_cur.h"

#include <iostream>
#include <sstream>
#include "makros.h"
#include "Ogs5FileTools.h"

namespace ogs5
{
#define RFD_FILE_EXTENSION ".rfd"                 //OK

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void CURRead(const std::string &base_file_name, std::vector<Kurven*> &kurven_vector)
{
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
   //----------------------------------------------------------------------
    // file handling
    std::string cur_file_name = base_file_name + RFD_FILE_EXTENSION;
    std::ifstream cur_file (cur_file_name.c_str(), std::ios::in);
    if (!cur_file.good())
        return;
    cur_file.seekg(0L,std::ios::beg);
    //========================================================================
    // keyword loop
//    std::cout << "CURRead" << std::endl;
    while (!cur_file.eof())
    {
        cur_file.getline(line,MAX_ZEILE);
        line_string = line;
        if(line_string.find("#STOP") != std::string::npos)
            return;
        //----------------------------------------------------------------------
        // keyword found
        if(line_string.find("#CURVE") != std::string::npos)
        {
            Kurven* kurven = new Kurven();
            position = kurven->Read(&cur_file);
            kurven_vector.push_back(kurven);
            cur_file.seekg(position,std::ios::beg);
        }                         // keyword found
    }                                     // eof
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
std::ios::pos_type Kurven::Read(std::ifstream* cur_file)
{
    bool new_keyword = false;
    std::string hash("#");
    std::string line_string;
    std::ios::pos_type position;
    std::stringstream line_stream;
    double d1,d2;
    //----------------------------------------------------------------------
    while (!new_keyword)
    {
        position = cur_file->tellg();
        line_string = readNonBlankLineFromInputStream(*cur_file);
        if(line_string.size() < 1)
            continue;
        //....................................................................
        // Test next keyword
        if(line_string.find(hash) != std::string::npos)
        {
            new_keyword = true;
            continue;
        }
        //--------------------------------------------------------------------
        if(line_string.find(";") != std::string::npos)
            continue;
        //--------------------------------------------------------------------
        //DATA
        line_stream.str(line_string);
        line_stream >> d1 >> d2;
        StuetzStellen* stuetz = new StuetzStellen();
        stuetzstellen.push_back(stuetz);
        stuetz->punkt = d1;
        stuetz->wert = d2;
        line_stream.clear();
        //--------------------------------------------------------------------
    }
    return position;
}

///**************************************************************************
//   FEMLib-Method:
//   04/2007 OK Implementation
//**************************************************************************/
//void CURWrite(std::vector<Kurven*> &kurven_vector)
//{
//    //========================================================================
//    // File handling
//    std::string fct_file_name = "test.cur";
//    std::fstream fct_file (fct_file_name.c_str(), std::ios::trunc | std::ios::out);
//    fct_file.setf(std::ios::scientific, std::ios::floatfield);
//    fct_file.precision(12);
//    if (!fct_file.good())
//        return;
//    fct_file << "GeoSys-CUR: Functions ------------------------------------------------" <<
//    std::endl;
//    //========================================================================
//    for(size_t i = 0; i < kurven_vector.size(); i++)
//    {
//        fct_file << "#CURVES" << std::endl;
//        for(size_t j = 0; j < kurven_vector[i]->stuetzstellen.size(); j++)
//        {
//            StuetzStellen* stuetz = kurven_vector[i]->stuetzstellen[j];
//            fct_file << stuetz->punkt << " " << stuetz->wert <<  std::endl;
//        }
//    }
//    fct_file << "#STOP";
//}

Kurven::~Kurven()
{
    for (size_t i=0; i<stuetzstellen.size(); i++) {
        if (stuetzstellen[i]!=0)
            delete stuetzstellen[i];
    }
}

double Kurven::getCurveValue(int methode, double punkt, int& gueltig)
{
    gueltig = 1;
    int i = 1l;

    if (punkt < stuetzstellen[0]->punkt) {
        gueltig = 0;
        return stuetzstellen[0]->wert;
    } else if (punkt > stuetzstellen.back()->punkt) {
        gueltig = 0;
        return stuetzstellen.back()->wert;
    }

    /* Suchen der Stuetzstelle. Vorraussetzung: Zeitpunkte aufsteigend geordnet */
    while (punkt > stuetzstellen[i]->punkt)
        i++;

    switch (methode)
    {
    default:
    case 0:
        /* Lineare Interpolation */
        return stuetzstellen[i - 1]->wert +
               (stuetzstellen[i]->wert - stuetzstellen[i - 1l]->wert) / (stuetzstellen[i]->punkt - stuetzstellen[i - 1l]->punkt) *
               (punkt - stuetzstellen[i - 1l]->punkt);
    case 1:                               //WW/SF
        // Piece wise constant
        return stuetzstellen[i]->wert;
    }
}

} //ogs5

