/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CsvInterface.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace BaseLib
{
namespace IO
{

CsvInterface::CsvInterface()
: _writeCsvHeader(true)
{
}

std::size_t CsvInterface::findColumn(std::string const& line, char delim, std::string const& column_name)
{
    std::list<std::string> const fields = BaseLib::splitString(line, delim);
    if (fields.empty())
        return std::numeric_limits<std::size_t>::max();

    std::size_t count(0);
    for (auto it = fields.cbegin(); it != fields.cend(); ++it)
    {
        if ((*it).compare(column_name) == 0)
            break;
        else
            count++;
    }

    if (count == fields.size())
        return std::numeric_limits<std::size_t>::max();

    return count;
}

void CsvInterface::addIndexVectorForWriting(std::size_t s)
{
    std::vector<int> idx_vec(s);
    std::iota(idx_vec.begin(), idx_vec.end(), 0);
    addVectorForWriting("Index", idx_vec);
}

bool CsvInterface::write()
{
    if (_data.empty())
    {
        ERR ("CsvInterface::write() - No data to write.");
        return false;
    }

    std::size_t const n_vecs (_data.size());
    std::size_t const vec_size (getVectorSize(0));

    if (_writeCsvHeader)
    {
        _out << _vec_names[0];
        for (std::size_t i=1; i<n_vecs; ++i)
            _out << "\t" << _vec_names[i];
        _out << "\n";
    }

    for (std::size_t j=0; j<vec_size; ++j)
    {
        writeValue(0,j);
        for (std::size_t i=1; i<n_vecs; ++i)
        {
            _out << "\t";
            writeValue(i,j);
        }
        _out << "\n";
    }
    return true;
}

std::size_t CsvInterface::getVectorSize(std::size_t idx) const
{
    if (_data[idx].type() == typeid(std::vector<std::string>))
        return boost::any_cast<std::vector<std::string>>(_data[idx]).size();
    else if (_data[idx].type() == typeid(std::vector<double>))
        return boost::any_cast<std::vector<double>>(_data[idx]).size();
    else if (_data[idx].type() == typeid(std::vector<int>))
        return boost::any_cast<std::vector<int>>(_data[idx]).size();
    return 0;
}

void CsvInterface::writeValue(std::size_t vec_idx, std::size_t in_vec_idx)
{
    if (_data[vec_idx].type() == typeid(std::vector<std::string>))
        _out << boost::any_cast<std::vector<std::string>>(_data[vec_idx])[in_vec_idx];
    else if (_data[vec_idx].type() == typeid(std::vector<double>))
        _out << boost::any_cast<std::vector<double>>(_data[vec_idx])[in_vec_idx];
    else if (_data[vec_idx].type() == typeid(std::vector<int>))
        _out << boost::any_cast<std::vector<int>>(_data[vec_idx])[in_vec_idx];
}

} // end namespace IO
} // end namespace BaseLib
