/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TINInterface.h"

#include <fstream>
#include <limits>

#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/Surface.h"


namespace FileIO
{

GeoLib::Surface* TINInterface::readTIN(std::string const& fname, std::vector<GeoLib::Point*> &pnt_vec, std::vector<std::string>* errors)
{
	// open file
	std::ifstream in(fname.c_str());
	if (!in) {
		WARN("readTIN(): could not open stream from %s.", fname.c_str());
		if (errors) errors->push_back ("readTINFile error opening stream from " + fname);
		return nullptr;
	}

	GeoLib::Surface* sfc = new GeoLib::Surface(pnt_vec);
	std::size_t id;
	double p0[3], p1[3], p2[0];
	while (in.good())
	{
		// read id
		if (!(in >> id))
			continue;
		// read first point
		if (!(in >> p0[0] >> p0[1] >> p0[2])) {
			ERR("Could not read coords of 1st point of triangle %d.", id);
			if (errors)
				errors->push_back (std::string("readTIN error: ") +
					std::string("Could not read coords of 1st point in triangle ") +
					std::to_string(id));
			in.close();
			delete sfc;
			return nullptr;
		}
		// read second point
		if (!(in >> p1[0] >> p1[1] >> p1[2])) {
			ERR("Could not read coords of 2nd point of triangle %d.", id);
			if (errors)
				errors->push_back (std::string("readTIN error: ") +
					std::string("Could not read coords of 2nd point in triangle ") +
					std::to_string(id));
			in.close();
			delete sfc;
			return nullptr;
		}
		// read third point
		if (!(in >> p2[0] >> p2[1] >> p2[2])) {
			ERR("Could not read coords of 3rd point of triangle %d.", id);
			if (errors)
				errors->push_back (std::string("readTIN error: ") +
					std::string("Could not read coords of 3rd point in triangle ") +
					std::to_string(id));
			in.close();
			delete sfc;
			return nullptr;
		}
		// determine size pnt_vec to insert the correct ids
		std::size_t pnt_pos(pnt_vec.size());
		pnt_vec.push_back(new GeoLib::Point(p0));
		pnt_vec.push_back(new GeoLib::Point(p1));
		pnt_vec.push_back(new GeoLib::Point(p2));
		// create new Triangle
		sfc->addTriangle(pnt_pos, pnt_pos + 1, pnt_pos + 2);
	}

	if (sfc->getNTriangles() == 0) {
		WARN("readTIN(): No triangle found.", fname.c_str());
		if (errors)
			errors->push_back ("readTIN error because of no triangle found");
		delete sfc;
		return nullptr;
	}

	return sfc;
}

void TINInterface::writeSurfaceAsTIN(GeoLib::Surface const& surface, std::string const& file_name)
{
	std::ofstream os (file_name.c_str());
	os.precision(std::numeric_limits<double>::digits10);
	const std::size_t n_tris (surface.getNTriangles());
	for (std::size_t l(0); l < n_tris; l++) {
		GeoLib::Triangle const& tri (*(surface[l]));
		os << l << " " << *(tri.getPoint(0)) << " " << *(tri.getPoint(1)) << " " << *(tri.getPoint(2)) << "\n";
	}
	os.close();
}

} // end namespace GeoLib
