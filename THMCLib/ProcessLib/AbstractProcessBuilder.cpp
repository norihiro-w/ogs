/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractProcessBuilder.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#include "AbstractProcessBuilder.h"

#include <iostream>
#include <logog/include/logog.hpp>

#include "Process.h"
#include "ProcessInfo.h"
#include "ProcessFactoryBase.h"

namespace THMCLib
{

ProcessInfo* AbstractProcessBuilder::registerProcess(const std::string &pcs_type, ProcessFactoryBase* pcs_buid)
{
    this->_map_pcs_name2new[pcs_type] = pcs_buid;
    //ProcessInfo* pcs_info = new ProcessInfo();
    //return pcs_info;

    return NULL;
}

bool AbstractProcessBuilder::hasRegisterd(const std::string &pcs_type) const
{
    return _map_pcs_name2new.count(pcs_type)>0;
}

Process* AbstractProcessBuilder::create(const std::string &pcs_type) const
{
    if (!hasRegisterd(pcs_type)) return 0;

    std::map<std::string, ProcessFactoryBase*>::const_iterator itr = this->_map_pcs_name2new.find(pcs_type);
    return itr->second->createProcess();
}

void AbstractProcessBuilder::output() const
{
    std::map<std::string, ProcessFactoryBase*>::const_iterator itr;
    INFO("------------------------------------------------------------------");
    INFO("List of available modules");
    for (itr=this->_map_pcs_name2new.begin(); itr!=_map_pcs_name2new.end(); ++itr) {
        INFO("* %s", itr->first.c_str());
    }
    INFO("------------------------------------------------------------------");
    INFO("");
}

} //end
