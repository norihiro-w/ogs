/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NamedIO.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "logog/include/logog.hpp"

#include "BaseLib/CodingTools.h"

namespace NumLib
{

/**
 * \brief Input-Output system with named parameters
 *
 * This class just manage parameter names but does not store parameters themselves.
 *
 * \tparam T_BASE A class which stores parameters
 */
template <typename T_BASE>
class NamedIO : public T_BASE
{
public:
    NamedIO() {};
    virtual ~NamedIO() {};

    /**
     *
     * @param i
     * @param key
     */
    virtual void setInputParameterName(size_t i, const std::string& key)
    {
        if (i+1>_in_para_names.size()) _in_para_names.resize(i+1);
        _in_para_names[i] = key;
    }

    /**
     *
     * @param i
     * @param key
     */
    virtual void setOutputParameterName(size_t i, const std::string& key)
    {
        if (i+1>_out_para_names.size()) _out_para_names.resize(i+1);
        _out_para_names[i] = key;
    }

    /**
     *
     * @return
     */
    virtual bool isValid() const
    {
        if (T_BASE::getNumberOfInputParameters()!=_in_para_names.size()) {
            ERR("***Error: the number of input parameter names is not equal to that of values.");
            return false;
        }
        if (T_BASE::getNumberOfOutputParameters()!=_out_para_names.size()) {
            ERR("***Error: the number of output parameter names is not equal to that of values.");
            return false;
        }
        return true;
    }

    /**
     *
     * @param key
     * @return
     */
    virtual bool hasInputParameter(const std::string& key) const
    {
        return _in_para_names.end() != std::find(_in_para_names.begin(), _in_para_names.end(), key);
    }

    virtual bool hasOutputParameter(const std::string& key) const
    {
        return _out_para_names.end() != std::find(_out_para_names.begin(), _out_para_names.end(), key);
    }

    /// get the name of parameter with the given name
    virtual const std::string& getInputParameterName(size_t i) const { return _in_para_names[i];};

    /// get the name of parameter with the given name
    virtual const std::string& getOutputParameterName(size_t i) const { return _out_para_names[i];};

    /// get the index of parameter with the given name
    virtual int getInputParameterID(const std::string&  key) const
    {
        std::vector<std::string>::const_iterator itr = std::find(_in_para_names.begin(), _in_para_names.end(), key);
        if (itr==_in_para_names.end()) {
            return -1;
        } else {
               return itr - _in_para_names.begin();
        }
    };

    /// get the index of parameter with the given name
    virtual int getOutputParameterID(const std::string&  key) const
    {
        std::vector<std::string>::const_iterator itr = std::find(_out_para_names.begin(), _out_para_names.end(), key);
        if (itr==_out_para_names.end()) {
            return -1;
        } else {
               return itr - _out_para_names.begin();
        }
    }

    /// get the number of output parameters
    virtual size_t getNumberOfInputParameterNames() const {return _in_para_names.size();};

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameterNames() const {return _out_para_names.size();};

private:
    std::vector<std::string> _in_para_names;
    std::vector<std::string> _out_para_names;
};

} //end
