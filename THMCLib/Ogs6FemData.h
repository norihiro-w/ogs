/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <string>

#include "BaseLib/OrderedMap.h"
#include "MathLib/LinAlg/LinAlgLibType.h"
#include "MathLib/DataType.h"
#include "GeoLib/GEOObjects.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"

#include "MaterialLib/IMedium.h"
#include "MaterialLib/PorousMediumModel.h"
#include "MaterialLib/FluidModel.h"
#include "MaterialLib/Compound.h"
#include "MaterialLib/SolidModel.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "THMCLib/ProcessLib/Process.h"
#include "THMCLib/Output/OutputController.h"
#include "FeElementData.h"


namespace THMCLib
{

/**
 * \brief Fem data storage
 *
 * This class follows singleton pattern.
 */
class Ogs6FemData
{
public:
    static Ogs6FemData* getInstance();
private:
    static Ogs6FemData* _obj;

private:
    Ogs6FemData(): geo(NULL) {}

public:
    //material data
    std::vector<MaterialLib::IMedium*> list_medium;
    std::vector<MaterialLib::PorousMediumModel*> list_pm;
    std::vector<MaterialLib::SolidModel*> list_solid;
    std::vector<MaterialLib::FluidModel*> list_fluid;
    std::vector<MaterialLib::Compound*> list_compound;
    // geometric data
    std::string geo_unique_name;
    GeoLib::GEOObjects* geo = nullptr;
    // mesh data
    std::vector<MeshLib::Mesh*> list_mesh;
    std::vector<MeshGeoToolsLib::MeshNodeSearcher*> list_nodeSearcher;
    std::vector<MeshGeoToolsLib::BoundaryElementsSearcher*> list_beSearcher;
    // time group data
    std::vector<NumLib::ITimeStepAlgorithm*> list_tim;
    //process
    BaseLib::OrderedMap<std::string, THMCLib::Process*> list_pcs;
    //
    ogs6::OutputController outController;
    //
    std::string project_name;
    std::string project_dir;
    std::string output_dir;
#ifdef USE_LIS
    MathLib::LinAlgLibType linalg_type = MathLib::LinAlgLibType::EigenLis;
#else
    MathLib::LinAlgLibType linalg_type = MathLib::LinAlgLibType::Eigen;
#endif
    //std::vector<std::vector<FeElementData>> fe_ele_data; // msh_id, (ele id, fe data)
    std::vector<FeMeshData> fe_mesh_data;

    void initialize()
    {
        delete geo;
        for (auto p : list_medium) delete p;
        for (auto p : list_solid) delete p;
        for (auto p : list_fluid) delete p;
        for (auto p : list_compound) delete p;
        for (auto p : list_mesh) delete p;
        for (auto p : list_tim) delete p;
        for (auto p : list_pcs) delete p.second;
        for (auto p : list_nodeSearcher) delete p;
        for (auto p : list_beSearcher) delete p;
    }

    ~Ogs6FemData()
    {
        initialize();
    }
};

} //THMCLib
