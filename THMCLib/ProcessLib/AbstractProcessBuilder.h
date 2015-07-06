/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractProcessBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <map>

#include "ProcessFactoryImpl.h"


namespace THMCLib
{
class Process;
class ProcessInfo;
class ProcessFactoryBase;

/**
 * \brief Abstract class for any process builder.
 *
 * Any subclass should follow singleton pattern.
 */
class AbstractProcessBuilder
{
public:
    virtual ~AbstractProcessBuilder() {}

    /// register a process to this builder
    ProcessInfo* registerProcess(const std::string &pcs_type, ProcessFactoryBase* pcs_buid);

    /// check if a process has been already registered
    bool hasRegisterd(const std::string &pcs_type) const;

    /// create the new instance of a specified process
    /// @return return 0 if a process is not found
    Process* create(const std::string &pcs_type) const;

    /// list available processes to standard IO
    void output() const;

protected:
    AbstractProcessBuilder() {}; // only called by sub classes

private:
    std::map<std::string, ProcessFactoryBase*> _map_pcs_name2new;
};

} //end

//#define PROCRSS_REGISTER
#ifdef PROCRSS_REGISTER

#define OGS_PROCESS(pcs_name, pcs_classname)\
ProcessLib::ProcessInfo* const pcs_classname::_pcs_info =\
        ProcessLib::ProcessBuilder::getInstance()->registerProcess(\
        #pcs_name,\
        new ProcessLib::ProcessFactoryImpl<pcs_classname >);

#define OGS_DEF_PROCESS(pcs_name, pcs_classname)\
int class_##pcs_name = 0;

#define OGS_LINK_PROCESS(pcs_name, pcs_classname)\
extern int class_##pcs_name; \
int link_##pcs_name = class_##pcs_name;

#else

#define OGS_PROCESS(pcs_name, pcs_classname)
#define OGS_DEF_PROCESS(pcs_name, pcs_classname)
#define OGS_LINK_PROCESS(pcs_name, pcs_classname)
#define OGS_ADD_PROCESS(pcs_name, pcs_classname)\
        registerProcess(#pcs_name, new THMCLib::ProcessFactoryImpl<pcs_classname>);
#define OGS_ADD_PROCESS_SYS(pcs_name, pcs_classname, discrete_system)\
        registerProcess(#pcs_name, new THMCLib::ProcessFactoryImpl<pcs_classname<discrete_system> >);
#define OGS_ADD_PROCESS_SYS_SOLVER(pcs_name, pcs_classname, discrete_system, linear_solver)\
        registerProcess(#pcs_name, new THMCLib::ProcessFactoryImpl<pcs_classname<discrete_system,linear_solver> >);

#endif

