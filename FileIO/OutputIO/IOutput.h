/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IOutput.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/OrderedMap.h"
#include "MeshLib/Mesh.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IOutputTiming.h"

/**
 * 
 */
struct OutputVariableInfo
{
    enum DataType
    {
        Char,
        Int,
        Real
    };
    enum ObjectType
    {
        Node,
        Element
    };

    std::string name;
    size_t mesh_id;
    ObjectType object_type;
    DataType data_type;
    size_t nr_of_components;
    NumLib::ITXFunction* value;

    OutputVariableInfo() 
        : name(""), mesh_id(0), object_type(Node), data_type(Char), nr_of_components(1), value(NULL)
    {};

    OutputVariableInfo(const std::string &s, size_t msh_id, ObjectType o, DataType d, size_t n, NumLib::ITXFunction* v)
    : name(s), mesh_id(msh_id), object_type(o), data_type(d), nr_of_components(n), value(v)
    {
    };
};

/**
 * \brief Interface of Output class
 */
class IOutput
{
public:
    ///
    IOutput() : _output_timing(0), _geo_obj(0), _msh(0) {};

    ///
    virtual ~IOutput()
    {
        //BaseLib::releaseObject(_output_timing);
    };

    ///
    void setOutputPath(const std::string &dir, const std::string &base_name)
    {
        _output_dir_path = dir;
        _base_name = base_name;
    }

    ///
    std::string getOutputDir() const {return _output_dir_path;};
    
    ///
    std::string getOutputBaseName() const {return _base_name;};

    ///
    void addNodalVariable(const std::string &var_name)
    {
        _list_nod_var_name.push_back(var_name);
    };

    ///
    void addNodalVariable(const std::vector<std::string> &vec_var_name)
    {
        _list_nod_var_name.insert(_list_nod_var_name.end(), vec_var_name.begin(), vec_var_name.end());
    };
    
    ///
    std::vector<std::string>& getListOfNodalVariables() {return _list_nod_var_name;};
    
    ///
    bool hasNodalVariable(const std::string &var_name) const
    {
        return (std::find(_list_nod_var_name.begin(), _list_nod_var_name.end(), var_name) != _list_nod_var_name.end());
    }

    ///
    void addElementalVariable(const std::string &var_name)
    {
        _list_ele_var_name.push_back(var_name);
    };

    ///
    void addElementalVariable(const std::vector<std::string> &vec_var_name)
    {
        _list_ele_var_name.insert(_list_ele_var_name.end(), vec_var_name.begin(), vec_var_name.end());
    };

    ///
    std::vector<std::string>& getListOfElementalVariables() {return _list_ele_var_name;};

    ///
    bool hasElementalVariable(const std::string &var_name) const
    {
        return (std::find(_list_ele_var_name.begin(), _list_ele_var_name.end(), var_name) != _list_ele_var_name.end());
    }

    ///
    void setGeometry(const GeoLib::GeoObject* geo_obj) {_geo_obj = geo_obj;};
    ///
    const GeoLib::GeoObject* getGeometry() const {return _geo_obj;};

    void setMesh(MeshLib::Mesh* msh) {_msh = msh;};
    MeshLib::Mesh* getMesh() const {return _msh;};

    ///
    void setOutputTiming(IOutputTiming* timing) {_output_timing = timing;};
    ///
    bool isActive(const NumLib::TimeStep &current_time) const {return _output_timing->isActive(current_time);};

    ///
    virtual void write( const NumLib::TimeStep &current_time, 
                        BaseLib::OrderedMap<std::string, OutputVariableInfo> &data) = 0;

private:
    std::string _output_dir_path;
    std::string _base_name;
    IOutputTiming* _output_timing;
    std::vector<std::string> _list_nod_var_name;
    std::vector<std::string> _list_ele_var_name;
    const GeoLib::GeoObject* _geo_obj;
    MeshLib::Mesh* _msh;
};

