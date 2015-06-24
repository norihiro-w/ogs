/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file rf_cur.h
 *
 * Created on 2007-04-xx by Olaf Kolditz
 *
 */

#ifndef RF_CUR_H_
#define RF_CUR_H_

#include <string>
#include <fstream>
#include <vector>

namespace ogs5
{

class Kurven
{
public:
    ~Kurven();

    /**
     *
     * @param cur_file
     * @return
     */
    std::ios::pos_type Read(std::ifstream* cur_file);

    /**
     * Aufgabe:
     * Liefert Wert aus einer Kurve fuer angegebenen Punkt.
     * Liegt der Punkt ausserhalb des durch die Kurve beschriebenen
     * Bereichs, wird der letzte bzw. erste Wert zurueckgeliefert und
     * der Flag fuer die Gueltigkeit auf 0 gesetzt.
     * @param methode   Interpolationsmethode
     * @param punkt     Punkt
     * @param gueltig   Flag fuer die Gueltigkeit des zurueckgelieferten Wertes
     * @return
     */
    double getCurveValue(int methode, double punkt, int& gueltig);

public:
    struct StuetzStellen
    {
        double punkt;
        double wert;
    };
    std::vector<StuetzStellen*> stuetzstellen;
};


void CURRead(const std::string &base_file_name, std::vector<Kurven*> &kurven_vector);

//void CURWrite(std::vector<Kurven*> &kurven_vector);

}
#endif
