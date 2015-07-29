/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <cstring>
#include <string>
#include <vector>

extern std::vector<std::string> g_args; // defined in main.cpp

inline char *convert(const std::string & s)
{
   char *pc = new char[s.size()+1];
   std::strcpy(pc, s.c_str());
   return pc;
}

inline std::vector<char*> createCmdLineArguments(const std::string &str_in_project, const std::string &str_out_dir)
{
    std::vector<std::string> vec_args;
    vec_args.push_back("ogs");
//    vec_args.push_back("-i");
    vec_args.push_back(str_in_project);
    vec_args.push_back("-o");
    vec_args.push_back(str_out_dir);
    for (std::size_t i=1; i<g_args.size(); i++)
        vec_args.push_back(g_args[i]);
    std::vector<char*> argv;
    std::transform(vec_args.begin(), vec_args.end(), std::back_inserter(argv), convert);
    return argv;
}
