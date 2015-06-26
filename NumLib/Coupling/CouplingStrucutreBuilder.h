/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "ICoupledProblem.h"
#include "MonolithicProblem.h"
#include "PartitionedProblem.h"
#include "Algorithm/PartitionedAlgorithmFactory.h"

//template<typename Ptree>
//void write_xml
//(
//	std::basic_ostream< typename Ptree::key_type::value_type > & stream,
//	const Ptree & pt,
//	const boost::property_tree::xml_writer_settings< typename Ptree::key_type::value_type > & settings =
//			boost::property_tree::xml_writer_settings< typename Ptree::key_type::value_type >()
//);

namespace NumLib
{

template <class T_I, class T_M, class T_P, class T_ALGORITHM>
class TemplateCouplingStrucutreBuilder
{
    std::vector<T_M*> _vec_m;
    std::vector<std::string> _vec_m_name;
public:
    TemplateCouplingStrucutreBuilder() {};

    template <class T_EQS_FACTORY>
    T_I* build(const boost::property_tree::ptree *option, T_EQS_FACTORY &eqs_fac)
    {
        auto op_cpl = option->get_child_optional("coupling");
//        boost::property_tree::xml_writer_settings<char> settings('\t', 1);
//        boost::property_tree::write_xml(std::cout, *op_cpl, settings);
//        std::cout << std::endl;
        if (!op_cpl) {
            INFO("tag<coupling> not found.");
            return NULL;
        }
        if (auto op_sub = op_cpl.get().get_child_optional("M")) {
            T_M *sys = buildMonolithicSystem(&op_sub.get(), eqs_fac);
            return sys;
        } else if (auto op_sub = op_cpl.get().get_child_optional("P")) {
            T_P *sys = buildPartitionedSystem(&op_sub.get(), eqs_fac);
            return sys;
        }
        return NULL;
    }

    std::vector<T_M*>& getListOfMonolithicSystem()
    {
        return _vec_m;
    }

    std::vector<std::string>& getListOfMonolithicSystemName()
    {
        return _vec_m_name;
    }

private:
    template <class T_EQS_FACTORY>
    T_M* buildMonolithicSystem(const boost::property_tree::ptree *option, T_EQS_FACTORY &eqs_fac)
    {
//		boost::property_tree::xml_writer_settings<char> settings('\t', 1);
//		boost::property_tree::write_xml(std::cout, *option, settings);
//		std::cout << std::endl;
        std::string eqs_type = option->get<std::string>("type");
        if (eqs_type.empty()) {
            ERR("***Error in TemplateCouplingStrucutreBuilder: attribute \"type\" is not found.");
            return NULL;
        }
        T_M* eqs = eqs_fac.create(eqs_type);
        if (eqs==NULL) {
            ERR("***Error in TemplateCouplingStrucutreBuilder: Monolithic system type <%s> is not found.", eqs_type.c_str());
            return NULL;
        }
        {
            auto it_range = option->equal_range("in");
            unsigned i = 0;
            for (auto v = it_range.first; v!=it_range.second; ++v)
                eqs->setInputParameterName(i++, (*v).second.data());
        }
        {
            auto it_range = option->equal_range("out");
            unsigned i = 0;
            for (auto v = it_range.first; v!=it_range.second; ++v)
                eqs->setOutputParameterName(i++, (*v).second.data());
        }
        _vec_m.push_back(eqs);
        std::string eqs_name = option->get_value<std::string>("name");
        if (eqs_name.empty())
            eqs_name = eqs_type;
        _vec_m_name.push_back(eqs_name);
        return eqs;
    }

    template <class T_EQS_FACTORY>
    T_P* buildPartitionedSystem(const boost::property_tree::ptree *option, T_EQS_FACTORY &eqs_fac)
    {
//        boost::property_tree::xml_writer_settings<char> settings('\t', 1);
//        boost::property_tree::write_xml(std::cout, *option, settings);
//        std::cout << std::endl;
        T_P* part = new T_P();
        {
            auto range = option->equal_range("in");
            unsigned i = 0;
            for (auto it=range.first; it!=range.second; ++it)
                part->setInputParameterName(i++, it->second.data());
            part->resizeInputParameter(i);
        }
        {
            auto range = option->equal_range("out");
            unsigned i = 0;
            for (auto it=range.first; it!=range.second; ++it)
                part->setOutputParameterName(i++, it->second.data());
            part->resizeOutputParameter(i);
        }
        //alg
        size_t max_itr = option->get<size_t>("max_itr");
        double epsilon = option->get<double>("epsilon");
        std::string alg = option->get<std::string>("algorithm");
        part->setAlgorithm(*T_ALGORITHM::create(alg, max_itr, epsilon));
//        IPartitionedAlgorithm* alg = T_ALGORITHM::create(option->getOption("algorithm"), checker, max_itr, epsilon);
//        if (alg!=0) {
//        }
        //problems
        if (auto op = option->get_child_optional("problems")) {
            for (auto& v : *op) {
                std::string str = v.first;
                auto &op_sub = v.second;
                T_I* sys = 0;
                if (str.compare("M")==0) {
                    sys = buildMonolithicSystem(&op_sub, eqs_fac);
                    if (sys!=0) part->addProblem(*sys);
                } else if (str.compare("P")==0) {
                    sys = buildPartitionedSystem(&op_sub, eqs_fac);
                    if (sys!=0) part->addProblem(*sys, true);
                }
            }
        }
        part->connectParameters();
        return part;
    }
};

typedef class TemplateCouplingStrucutreBuilder<ICoupledSystem,TemplateSteadyMonolithicSystem,PartitionedProblem,PartitionedAlgorithmFactory> CouplingStrucutreBuilder;


} //end


