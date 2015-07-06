/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ComparisonTools.h
 *
 * Created on 2012-12-10 by Norihiro Watanabe
 */

#pragma once

#include <string>

bool compareAciiVTUFiles(const std::string &str_ref_file, const std::string &str_test_file, const double epsilon);
