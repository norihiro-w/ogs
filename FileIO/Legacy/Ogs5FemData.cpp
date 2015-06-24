/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FemData.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "Ogs5FemData.h"


#include "BaseLib/CodingTools.h"

namespace BaseLib
{
template <typename T>
void releaseObjectsInStdVector(std::vector<T*> &vec)
{
	for (auto p :vec) delete p;
}
}

namespace ogs5
{

Ogs5FemData::~Ogs5FemData()
{
    BaseLib::releaseObjectsInStdVector(pcs_vector);
    BaseLib::releaseObjectsInStdVector(mfp_vector);
    BaseLib::releaseObjectsInStdVector(msp_vector);
    BaseLib::releaseObjectsInStdVector(mmp_vector);
    BaseLib::releaseObjectsInStdVector(cp_vector);
    BaseLib::releaseObjectsInStdVector(bc_vector);
    BaseLib::releaseObjectsInStdVector(st_vector);
    BaseLib::releaseObjectsInStdVector(ic_vector);
    BaseLib::releaseObjectsInStdVector(out_vector);
    BaseLib::releaseObjectsInStdVector(time_vector);
    BaseLib::releaseObjectsInStdVector(num_vector);
    //geo, msh objects are passed
	BaseLib::releaseObjectsInStdVector(KinReact_vector);
	BaseLib::releaseObjectsInStdVector(KinReactData_vector);
	BaseLib::releaseObjectsInStdVector(KinBlob_vector);

	BaseLib::releaseObjectsInStdVector(kurven_vector);
}

}
