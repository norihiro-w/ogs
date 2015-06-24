/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FemData.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "rf_pcs.h"
#include "rf_mfp_new.h"
#include "rfmat_cp.h"
#include "rf_bc_new.h"
#include "rf_st_new.h"
#include "rf_ic_new.h"
#include "rf_out_new.h"
#include "rf_tim_new.h"
#include "rf_msp_new.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"
#include "Output.h"
#include "rf_kinreact.h"
#include "rf_cur.h"

namespace ogs5
{

class Ogs5FemData
{

public:
    Ogs5FemData() : geo_obj(nullptr) {};
    ~Ogs5FemData();

public:
    std::vector<CRFProcess*> pcs_vector;
    std::vector<CFluidProperties*> mfp_vector;
    std::vector<CompProperties*> cp_vector;
    std::vector<CBoundaryCondition*> bc_vector;
    std::vector<CSourceTerm*> st_vector;
    std::vector<CInitialCondition*> ic_vector;
    std::vector<COutput*> out_vector;
    std::vector<CTimeDiscretization*> time_vector;
    std::vector<CSolidProperties*> msp_vector;
    std::vector<CMediumProperties*> mmp_vector;
    std::vector<CNumerics*>num_vector;
    std::vector<Kurven*> kurven_vector;
    GeoLib::GEOObjects* geo_obj;
    std::string geo_unique_name;
    std::vector<MeshLib::Mesh*> list_mesh;
	// this block is for kinetic reactions ----------------------------------------------------
	std::vector<CKinReact*>     KinReact_vector;               // declare instance CKinReact_vector
	std::vector<CKinReactData*> KinReactData_vector;
    std::vector<CKinBlob*>      KinBlob_vector;    // declare extern instance of class Blob
	// ----------------------------------------------------------------------------------------
	
};

}
