/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ROWCOLUMNINDICES_H_
#define ROWCOLUMNINDICES_H_

#include <iostream>
#include <vector>

namespace MathLib
{

template <typename IDX_TYPE>
struct RowColumnIndices
{
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_)
		: rows(rows_), columns(columns_)
	{ }

	LineIndex const& rows;
	LineIndex const& columns;

	void write(std::ostream &os) const
	{
		os << "rows: "; for (auto i : rows) os << i << " "; os << ", ";
		os << "cols: "; for (auto i : columns) os << i << " "; os << "\n";
	}
};

template <typename IDX_TYPE>
std::ostream& operator<< (std::ostream &os, const RowColumnIndices<IDX_TYPE> &p)
{
	p.write (os);
	return os;
}

} // MathLib

#endif  // ROWCOLUMNINDICES_H_
