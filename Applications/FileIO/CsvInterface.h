/**
 * @file CsvInterface.h
 * @author Karsten Rink
 * @date 2015-03-25
 * @brief Definition of the CsvInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CSVINTERFACE_H_
#define CSVINTERFACE_H_

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "BaseLib/IO/CsvInterface.h"

namespace GeoLib {
    class Point;
}

namespace FileIO {

/**
 * Interface for reading CSV file formats.
 */
class CsvInterface  : public BaseLib::IO::CsvInterface
{

public:
    CsvInterface() {}

    /**
     * Reads 3D points from a CSV file. It is assumed that the file has a header
     * specifying a name for each of the columns. The first three columns will be
     * interpreted as x-, y- and z-coordinate, respectively.
     * \param fname    Name of the file to be read
     * \param delim    Deliminator, default is ','
     * \param points   A vector containing the 3D points read from the file
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points);

    /**
     * Reads 3D points from a CSV file. It is assumed that the file has a header
     * specifying a name for each of the columns. The columns specified in the
     * function call will be used for reading x-, y- and z-coordinates,
     * respectively If z_column_name is an empty string or not given at all, all
     * z-coordinates will be set to zero.
     * \param fname           Name of the file to be read
     * \param delim           Deliminator, default is ','
     * \param points          A vector containing the 3D points read from the file
     * \param x_column_name   Name of the column to be interpreted as x-coordinate
     * \param y_column_name   Name of the column to be interpreted as y-coordinate
     * \param z_column_name   Name of the column to be interpreted as z-coordinate
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::string const& x_column_name,
                          std::string const& y_column_name,
                          std::string const& z_column_name = "");

    /**
     * Reads 3D points from a headerless CSV file, so columns for x-, y- and
     * z-coordinates have to be specified using indices (starting with 0).
     * If z_column_idx is not given (or set to numeric_limits::max()), all
     * z-coordinates will be set to zero.
     * \param fname          Name of the file to be read
     * \param delim          Deliminator, default is ','
     * \param points         A vector containing the 3D points read from the file
     * \param x_column_idx   Index of the column to be interpreted as x-coordinate
     * \param y_column_idx   Index of the column to be interpreted as y-coordinate
     * \param z_column_idx   Index of the column to be interpreted as z-coordinate
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::size_t x_column_idx,
                          std::size_t y_column_idx,
                          std::size_t z_column_idx = std::numeric_limits<std::size_t>::max());

private:
    /// Actual point reader for public methods
    static int readPoints(std::ifstream &in, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::array<std::size_t, 3> const& column_idx);

};

} // FileIO

#endif /* CSVINTERFACE_H_ */
