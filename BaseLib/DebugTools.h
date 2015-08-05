/**
 * \author Norihiro Watanabe
 * \date   2014-03-07
 * \brief  Helper tools for debugging
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DEBUGTOOLS_H
#define DEBUGTOOLS_H

#include <ostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <vector>

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
//    os << "\n";
    return os;
}

namespace BaseLib
{

template<typename T>
std::string toString(const std::vector<T> &v) {
    std::stringstream ss;
    ss << v;
    return ss.str();
}

} // BaseLib

#endif //DEBUGTOOLS_H

