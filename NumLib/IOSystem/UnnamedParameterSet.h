/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>

#include "BaseLib/CodingTools.h"
#include "Parameter.h"

namespace NumLib
{

/**
 * \brief Container of named variables
 */
class UnnamedParameterSet
{
public:
    ///
    UnnamedParameterSet() {};

    ///
    virtual ~UnnamedParameterSet() { clear(); }

    /// reset data
    void clear()
    {
        _list_var_data.clear();
        _list_attributes.clear();
        for (auto*p : _list_own_data) delete p;
    }

    /// move ownership of memory to another and initialize this
    void move(UnnamedParameterSet &dest)
    {
        dest.clear();
        dest._list_attributes.assign(this->_list_attributes.begin(), this->_list_attributes.end());
        dest._list_var_data.assign(this->_list_var_data.begin(), this->_list_var_data.end());
        // move ownership of memory
        dest._list_own_data.assign(this->_list_own_data.begin(), this->_list_own_data.end());
        _list_own_data.clear();
    }

    /// return the size of the stores variables
    size_t size() const { return _list_var_data.size(); }

    /// return if the container has a variable with the given key
    bool contain(size_t key) const { return (this->find(key)>=0); }

    /// return the name of the variable with the given index
    size_t getKey(size_t var_id) const { return _list_attributes[var_id].key; };

    /// return the index of the variable with the given name
    int find(const size_t key) const
    {
        for (std::vector<ParameterAttribute>::const_iterator itr = _list_attributes.begin(); itr!=_list_attributes.end(); ++itr) {
            if (itr->key == key) {
                return (itr - _list_attributes.begin());
            }
        }

        return -1;
    }

    /// return the variable with the given index
    const Parameter* get(size_t var_id) const
    {
        return _list_var_data[var_id];
    }

    /// return the variable with the given index
    template <class T>
    const T* get(size_t var_id) const
    {
        return static_cast<const T*>(_list_var_data[var_id]);
    }
    /// return the variable with the given index
    template <class T>
    T* get(size_t var_id)
    {
        return const_cast<T*>(static_cast<const T*>(_list_var_data[var_id]));
    }

    /// register the variable
    size_t add(size_t key, bool is_fixed=false)
    {
        size_t new_id = _list_attributes.size();
        _list_attributes.resize(_list_attributes.size()+1);
        _list_attributes[new_id].key = key;
        _list_attributes[new_id].is_fxied = is_fixed;
        _list_var_data.push_back(nullptr);
        return new_id;
    }

    /// set value
    void set(size_t var_id,  const Parameter* v)
    {
        //assert(&v!=0);
        if (v==nullptr) return;
        if (isFixed(var_id)) {
            _list_var_data[var_id] = v;
        } else {
            _list_var_data[var_id] = v->clone();
            _list_own_data.push_back(_list_var_data[var_id]);
        }
        _list_attributes[var_id].is_updated = true;
    }

    bool isFixed(size_t var_id) const
    {
        return _list_attributes[var_id].is_fxied;
    }

    void setFixed(size_t var_id, bool is_fxied)
    {
        _list_attributes[var_id].is_fxied = is_fxied;
    } 

    bool isUpdated(size_t var_id) const
    {
        return _list_attributes[var_id].is_updated;
    }

    void finishUpdate()
    {
        for (size_t i=0; i<_list_attributes.size(); i++)
            _list_attributes[i].is_updated = false;
    } 

    void setName(size_t var_id,std::string &name)
    {
        _list_attributes[var_id].name = name;
    }

    const std::string& getName(size_t var_id) const
    {
        return _list_attributes[var_id].name;
    }

private:
    struct ParameterAttribute
    {
        size_t key;
        std::string name;
        bool is_fxied;
        bool is_updated;
        ParameterAttribute() : key(0), is_fxied(false), is_updated(false) {};
    };

    std::vector<const Parameter*> _list_var_data;
    std::vector<ParameterAttribute> _list_attributes;
    std::vector<const Parameter*> _list_own_data;

    DISALLOW_COPY_AND_ASSIGN(UnnamedParameterSet);
};

}
