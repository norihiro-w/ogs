/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DynamicIOSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include <cassert>

#include "BaseLib/CodingTools.h"
#include "UnnamedParameterSet.h"

namespace NumLib
{

/**
 * \brief Input-Output system with the variable number of parameters
 *
 * - Internal ID
 * - Input ID
 * - Output ID
 *
 * @tparam T_BASE A base class inherited of ISystemWithInOutParameters
 */
template <typename T_BASE>
class DynamicIOSystem : public T_BASE
{
public:
    typedef size_t InternalID;

    DynamicIOSystem() {}
    virtual ~DynamicIOSystem() {}

    /// get the number of input parameters
    virtual size_t getNumberOfInputParameters() const { return _list_input_para_id.size(); }

    /// set input parameter
    virtual void setInput(InternalID id, const Parameter* val) { _shared_parameters.set(id, val); }

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameters() const {return _list_output_para_id.size();}

    /// get output parameter
    virtual const Parameter* getOutput(InternalID parameter_id) const
    {
        return _shared_parameters.get(getInternalIDFromOutputID(parameter_id));
    }

    // original functions in this class
    /// get output parameter
    template <class T>
    const T* getOutput(InternalID parameter_id) const
    {
        return _shared_parameters.get<T>(getInternalIDFromOutputID(parameter_id));
    }

    ///
    void resizeInputParameter(size_t n)
    {
        for (size_t i=0; i<n; i++)
            registerInputParameter(i);
    }

    ///
    void resizeOutputParameter(size_t n)
    {
        for (size_t i=0; i<n; i++)
            registerOutputParameter(i);
    }

    /// register parameter
    InternalID registerInputParameter(size_t key)
    {
        if (_shared_parameters.contain(key)) {
            std::cout << "***Error: the given key already exist." << std::endl;
            return _shared_parameters.find(key);
        }
        InternalID id = 0;
        id = _shared_parameters.add(key, true);
        _list_input_para_id.push_back(id);
        return id;
    }

    ///
    InternalID registerOutputParameter(size_t key)
    {
        size_t out_key = create_output_key(key);
        if (_shared_parameters.contain(out_key)) {
            std::cout << "***Error: the given key already exist." << std::endl;
            return _shared_parameters.find(out_key);
        }
        InternalID id = 0;
        id = _shared_parameters.add(out_key);
        _list_output_para_id.push_back(id);
        return id;
    }

    ///
    inline size_t getInternalIDFromOutputID(size_t i) const {return _list_output_para_id[i];} //TODO

    ///
    inline size_t getInternalIDFromInputID(size_t i) const {return _list_input_para_id[i];} //TODO

protected:
    inline size_t create_output_key(size_t i) const {return i+1000;} //TODO

    /// get input parameter
    const Parameter* getInput(InternalID parameter_id) const
    {
        return _shared_parameters.get(parameter_id);
    }

    /// get input parameter
    template <class T>
    const T* getInput(InternalID parameter_id) const
    {
        return _shared_parameters.get<T>(parameter_id);
    }

    /// set output parameter
    void setOutput(InternalID parameter_id, const Parameter* val)
    {
        _shared_parameters.set(getInternalIDFromOutputID(parameter_id), val);
    }

    /// get parameter set
    UnnamedParameterSet* getParameters() {return &_shared_parameters;}

private:
    UnnamedParameterSet _shared_parameters;
    std::vector<InternalID> _list_input_para_id;
    std::vector<InternalID> _list_output_para_id;

    DISALLOW_COPY_AND_ASSIGN(DynamicIOSystem);
};

#if 0
/**
 * \brief Implementation of ISystemWithInOutParameters
 *
 * @tparam T_BASE A base class inherited of ISystemWithInOutParameters
 */
template <typename T_BASE>
class SystemWithSharedParameters : public T_BASE
{
public:
    typedef size_t InternalID;
    typedef size_t ExternalKey;

    SystemWithSharedParameters() {};
    virtual ~SystemWithSharedParameters() {};

    /// get the number of input parameters
    virtual size_t getNumberOfInputParameters() const { return _list_input_para_id.size(); };

    /// get the name of parameter with the given name
    virtual ExternalKey getInputParameterKey(InternalID i) const { return _list_input_para_id[i];};

    /// get the index of parameter with the given name
    virtual int getInputParameterID(ExternalKey key) const { return _shared_parameters.find(key);};

    /// set input parameter
    virtual void setInput(InternalID id, const Parameter* val) { _shared_parameters.set(id, *val); };

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameters() const {return _list_output_para_id.size();};

    /// get the name of parameter with the given name
    virtual ExternalKey getOutputParameterKey(InternalID i) const { return _list_output_para_id[i];};

    /// get the index of parameter with the given name
    virtual int getOutputParameterID(ExternalKey key) const { return _shared_parameters.find(key);};

    /// get output parameter
    virtual const Parameter* getOutput(InternalID parameter_id) const
    {
        return _shared_parameters.get(parameter_id);
    }

    // original functions in this class
    /// get output parameter
    template <class T>
    const T* getOutput(InternalID parameter_id) const
    {
        return _shared_parameters.get<T>(parameter_id);
    }

    /// register parameter
    InternalID registerInputParameter(ExternalKey key)
    {
        if (_shared_parameters.contain(key)) {
            std::cout << "***Error: the given key already exist." << std::endl;
            return _shared_parameters.find(key);
        }
        InternalID id = 0;
        id = _shared_parameters.add(key, true);
        _list_input_para_id.push_back(id);
        return id;
    }

    InternalID registerOutputParameter(ExternalKey key)
    {
        if (_shared_parameters.contain(key)) {
            std::cout << "***Error: the given key already exist." << std::endl;
            return _shared_parameters.find(key);
        }
        InternalID id = 0;
        id = _shared_parameters.add(key);
        _list_output_para_id.push_back(id);
        return id;
    }

protected:
    /// get input parameter
    const Parameter* getInput(InternalID parameter_id) const
    {
        return _shared_parameters.get(parameter_id);
    }

    /// get input parameter
    template <class T>
    const T* getInput(InternalID parameter_id) const
    {
        return _shared_parameters.get<T>(parameter_id);
    }

    /// set output parameter
    void setOutput(InternalID parameter_id, const Parameter* val)
    {
        _shared_parameters.set(parameter_id, *val);
    }

    /// get parameter set
    UnnamedParameterSet* getParameters() {return &_shared_parameters;};

private:
    UnnamedParameterSet _shared_parameters;
    std::vector<InternalID> _list_input_para_id;
    std::vector<InternalID> _list_output_para_id;

    DISALLOW_COPY_AND_ASSIGN(SystemWithSharedParameters);
};
#endif

} //end
