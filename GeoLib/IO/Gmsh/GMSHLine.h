/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHLINE_H_
#define GMSHLINE_H_

#include <iosfwd>

namespace GeoLib
{
namespace IO
{
namespace GMSH
{

class GMSHLine {
public:
    GMSHLine(std::size_t start_point_id, std::size_t end_point_id);
    virtual ~GMSHLine();
    void write(std::ostream &os, std::size_t id) const;
    void resetLineData(std::size_t start_point_id, std::size_t end_point_id);

private:
    std::size_t _start_pnt_id;
    std::size_t _end_pnt_id;
};

} // end namespace GMSH
} // end namespace IO
} // end namespace GeoLib

#endif /* GMSHLINE_H_ */
