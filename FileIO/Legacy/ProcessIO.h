/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessIO.h
 *
 * Created on 2011-04-19 by Thomas Fisher
 */

#ifndef PROCESSIO_H_
#define PROCESSIO_H_

// STL
#include <iostream>

// FEM
#include "FEMEnums.h"

namespace ogs5
{

namespace FileIO
{
/**
 * Small class to read process information.
 */
class ProcessIO
{
public:
    /**
     * read process information from various file type:
     * - boundary condition files
     * - source term files
     *
     * To store the information an object of type ProcessInfo is used (CSourceTerm,
     * CBoundaryCondition, CInitialCondition and COutput inherit from class ProcessInfo).
     * @param in_str the input stream
     * @param pcs_type the process type
     * @return false, if the process is of type INVALID_PROCESS, else true
     * \sa enum ProcessType for valid values
     */
    static bool readProcessInfo (std::istream& in_str, FiniteElement::ProcessType& pcs_type);
};
}
}

#endif /* PROCESSIO_H_ */
