/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IOSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cassert>

#include "BaseLib/CodingTools.h"
#include "Parameter.h"

namespace NumLib
{


/**
 * \brief Input-Output system with the fixed number of parameters
 *
 * @tparam T_BASE A base class inherited of ISystemWithInOutParameters
 */
template <typename T_BASE>
class IOSystem : public T_BASE
{
public:
    IOSystem() {};
    virtual ~IOSystem() {};

    /// get the number of input parameters
    virtual size_t getNumberOfInputParameters() const { return _in_parameters.size(); };

    /// set input parameter
    virtual void setInput(size_t i, const Parameter* val)
    {
        _in_parameters[i] = val;
    }

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameters() const {return _out_parameters.size();};

    /// get output parameter
    virtual const Parameter* getOutput(size_t i) const
    {
        return _out_parameters[i];
    }

    // original functions in this class
    /// get output parameter
    template <class T>
    const T* getOutput(size_t i) const
    {
        return static_cast<const T*>(_out_parameters[i]);
    }

protected:
    /// register input parameter
    void resizeInputParameter(size_t n)
    {
        _in_parameters.resize(n, 0);
    }

    /// register output parameter
    void resizeOutputParameter(size_t n)
    {
        _out_parameters.resize(n, 0);
    }

    /// get input parameter
    const Parameter* getInput(size_t i) const
    {
        return _in_parameters[i];
    }

    /// get input parameter
    template <class T>
    const T* getInput(size_t i) const
    {
        return static_cast<const T*>(_in_parameters[i]);
    }

    /// set output parameter
    void setOutput(size_t i, const Parameter* val)
    {
        _out_parameters[i] = val;
    }

private:
    DISALLOW_COPY_AND_ASSIGN(IOSystem);

private:
    std::vector<const Parameter*> _in_parameters;
    std::vector<const Parameter*> _out_parameters;
};



} //end
