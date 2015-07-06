/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <string>

#include <boost/any.hpp>

#include "BaseLib/OrderedMap.h"
#include "MathLib/DataType.h"
#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "NumLib/Fem/Integration/IIntegration.h"
#include "THMCLib/MaterialLib/Models.h"

namespace THMCLib
{
struct FeElementData
{
    typedef MathLib::LocalMatrix ValueArray;
    MathLib::LocalVector p0;
    MathLib::LocalVector p1;
    std::vector<MathLib::LocalVector> v0;
    std::vector<MathLib::LocalVector> v1;
    MathLib::LocalVector T0;
    MathLib::LocalVector T1;

    void reset(const MeshLib::Element &e, unsigned n_gp)
    {
        p0.resize(e.getNNodes());
        p1.resize(e.getNNodes());
        v0.resize(n_gp);
        v1.resize(n_gp);
        T0.resize(e.getNNodes());
        T1.resize(e.getNNodes());
    }

    std::vector<ValueArray> node_values; // varId -> node value vector
    std::vector<ValueArray> ele_values;
    std::vector<ValueArray> gp_values;
    static std::vector<std::string> node_value_names;
    static std::vector<std::string> ele_value_names;
    static std::vector<std::string> gp_value_names;

    template <class T>
    static std::size_t findIndex(const std::vector<T> &vec, const std::string &name)
    {
        auto itr = std::find(vec.begin(), vec.end(), name);
        std::size_t index = std::distance(vec.begin(), itr);
        return (index < vec.size()) ? index : -1;
    }
    static std::size_t getNodeValueIndex(const std::string &name) {return findIndex(node_value_names, name);}
    static std::size_t getElementValueIndex(const std::string &name) {return findIndex(ele_value_names, name);}
    static std::size_t getGpValueIndex(const std::string &name) {return findIndex(gp_value_names, name);}
};

struct FeMeshData
{
    BaseLib::OrderedMap<std::string, boost::any> node_values;
    BaseLib::OrderedMap<std::string, boost::any> element_values;
    BaseLib::OrderedMap<std::string, boost::any> gp_values;
};

inline void updateFeData(FeMeshData &mshData, const NumLib::IFiniteElement &fe, FeElementData &fe_data)
{
    auto& e = *fe.getMeshElement();
    fe_data.reset(e, fe.getIntegrationMethod().getNPoints());
    if (mshData.node_values.count("PRESSURE1")>0) {
        auto &vec0 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("PRESSURE0")->second);
        auto &vec1 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("PRESSURE1")->second);
        for (unsigned i=0; i<e.getNNodes(); i++) {
            fe_data.p0[i] = vec0[e.getNodeIndex(i)];
            fe_data.p1[i] = vec1[e.getNodeIndex(i)];
        }
    }
    if (mshData.gp_values.count("VELOCITY1")>0) {
        typedef std::vector<std::vector<MathLib::LocalVector>> GpVector;
        auto &vec0 = boost::any_cast<GpVector&>(mshData.gp_values.find("VELOCITY0")->second);
        auto &vec1 = boost::any_cast<GpVector&>(mshData.gp_values.find("VELOCITY1")->second);
        const auto n_gp = fe.getIntegrationMethod().getNPoints();
        if (vec0[e.getID()].empty()) {
            vec0[e.getID()].resize(n_gp);
            vec1[e.getID()].resize(n_gp);
        } else {
            for (unsigned i=0; i<n_gp; i++) {
                fe_data.v0[i] = vec0[e.getID()][i];
                fe_data.v1[i] = vec1[e.getID()][i];
            }
        }
    }
    if (mshData.node_values.count("TEMPERATURE1")>0) {
        auto &vec0 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("TEMPERATURE0")->second);
        auto &vec1 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("TEMPERATURE1")->second);
        for (unsigned i=0; i<e.getNNodes(); i++) {
            fe_data.T0[i] = vec0[e.getNodeIndex(i)];
            fe_data.T1[i] = vec1[e.getNodeIndex(i)];
        }
    }

}

inline MaterialLib::StateVariables getStateVariables(const FeElementData &fe_data, const NumLib::DynamicShapeMatrices& sh)
{
    MaterialLib::StateVariables var;
    if (fe_data.p1.size()>0) var.p = sh.N.transpose()*fe_data.p1;
    if (fe_data.T1.size()>0) var.T = sh.N.transpose()*fe_data.T1;
    return var;
}
} //THMCLib
