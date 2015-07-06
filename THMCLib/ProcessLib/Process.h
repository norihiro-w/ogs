/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Process.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"

namespace THMCLib
{

/**
 * \brief Interface of Process object
 *
 */
class Process : public NumLib::AbstractTransientMonolithicSystem
{
public:
    /**
     *
     * @param pcs_name          default process name
     * @param n_in_parameters   the number of input parameters
     * @param n_out_parameters  the number of output parameters
     */
    Process(const std::string &pcs_type, size_t n_in_parameters, size_t n_out_parameters)
        : _pcs_type(pcs_type)
    {
        resizeInputParameter(n_in_parameters);
        resizeOutputParameter(n_out_parameters);
    }

    /**
     *
     */
    virtual ~Process() {};

    /**
     * get type name of this process
     * @return
     */
    const std::string& getProcessType() const {return _pcs_type;};

    /**
     * set this process name
     * @param pcs_name
     */
    void setProcessName(const std::string &pcs_name)
    {
        _pcs_name = pcs_name;
    }

    /**
     * set this process name
     * @return
     */
    const std::string& getProcessName() const {return _pcs_name;};

private:
    const std::string _pcs_type;
    std::string _pcs_name;
};


}

