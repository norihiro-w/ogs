/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Ogs5ToOgs6.h"

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

#include "NumLib/Fem/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Function/TXFunctionConstant.h"
#include "NumLib/Function/TXFunctionType.h"
#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"

#include "MaterialLib/FluidModel.h"
#include "MaterialLib/FluidDensity.h"
#include "MaterialLib/FluidViscosity.h"
#include "MaterialLib/Fracture.h"
#include "MaterialLib/SolidModel.h"

#include "GeoProcessBuilder.h"

using namespace ogs5;
using namespace boost::property_tree;

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertTimeStepping(const std::vector<CTimeDiscretization*> &td_vector, std::vector<NumLib::ITimeStepAlgorithm*> &tf_vector)
{
//    bool done_shared = false;
    for (size_t i=0; i<td_vector.size(); i++)
    {
        const CTimeDiscretization* td = td_vector[i];
//        if (td->time_independence || !done_shared) {
            NumLib::ITimeStepAlgorithm* tf = new NumLib::FixedTimeStepping(td->time_start, td->time_end, td->time_step_vector);
            tf_vector.push_back(tf);
//            if (!td->time_independence)
//                done_shared = true;
//        }
    }
}

void convertFluidProperty(const CFluidProperties &mfp, MaterialLib::FluidModel &fluid)
{
    MaterialLib::FluidDensity::Type density_type = MaterialLib::FluidDensity::Type::Constant;
    std::vector<double> density_parameters;
    if (mfp.density_model==1) {
        density_type = MaterialLib::FluidDensity::Type::Constant;
        density_parameters.push_back(mfp.rho_0);
        fluid.density = new MaterialLib::FluidDensity(density_type, density_parameters);
    }

    MaterialLib::FluidViscosity::Type viscosity_type = MaterialLib::FluidViscosity::Type::Constant;
    std::vector<double> viscosity_parameters;
    if (mfp.viscosity_model==1) {
        viscosity_type = MaterialLib::FluidViscosity::Type::Constant;
        viscosity_parameters.push_back(mfp.my_0);
        fluid.viscosity = new MaterialLib::FluidViscosity(viscosity_type, viscosity_parameters);
    }

    MaterialLib::FluidSpecificHeat::Type cp_type = MaterialLib::FluidSpecificHeat::Type::Constant;
    std::vector<double> cp_parameters;
    if (mfp.heat_capacity_model==1) {
        cp_type = MaterialLib::FluidSpecificHeat::Type::Constant;
        cp_parameters.push_back(mfp.specific_heat_capacity);
        fluid.specific_heat = new MaterialLib::FluidSpecificHeat(cp_type, cp_parameters);
    }

    MaterialLib::FluidThermalConductivity::Type lambda_type = MaterialLib::FluidThermalConductivity::Type::Constant;
    std::vector<double> lambda_parameters;
    if (mfp.heat_conductivity_model==1) {
        lambda_type = MaterialLib::FluidThermalConductivity::Type::Constant;
        lambda_parameters.push_back(mfp.heat_conductivity);
        fluid.thermal_conductivity = new MaterialLib::FluidThermalConductivity(lambda_type, lambda_parameters);
    }
}

void convertSolidProperty(const CSolidProperties &msp, MaterialLib::SolidModel &solid)
{
    if (msp.Density_mode==1) {
        solid.density = new MaterialLib::SolidDensity(MaterialLib::SolidDensity::Type::Constant, (*msp.data_Density)(0));
    }

    solid.poisson_ratio = new MaterialLib::PoissonRatio(MaterialLib::PoissonRatio::Type::Constant, msp.PoissonRatio);

    if (msp.Youngs_mode==1) {
        solid.Youngs_modulus = new MaterialLib::YoungsModulus(MaterialLib::YoungsModulus::Type::Constant, (*msp.data_Youngs)(0));
    }

    if (msp.Capacity_mode==1) {
        solid.specific_heat = new MaterialLib::SolidSpecificHeat(MaterialLib::SolidSpecificHeat::Type::Constant, (*msp.data_Capacity)(0));
    }

    if (msp.Conductivity_mode==1) {
        solid.thermal_conductivity = new MaterialLib::SolidThermalConductivity(MaterialLib::SolidThermalConductivity::Type::Constant, msp.data_Conductivity->Rows(), (*msp.data_Conductivity)(0));
    }
}

void convertPorousMediumProperty(const CMediumProperties &mmp, MaterialLib::PorousMediumModel &pm)
{
    if (mmp.permeability_model==1) {
        pm.hydraulic_conductivity = new MaterialLib::HydraulicConductivity(MaterialLib::HydraulicConductivity::Type::Constant, mmp.geo_dimension, mmp.permeability_tensor[0]);
        pm.permeability = new MaterialLib::Permeability(MaterialLib::Permeability::Type::Constant, mmp.geo_dimension, mmp.permeability_tensor[0]);
    }

    if (mmp.porosity_model==1) {
        pm.porosity = new MaterialLib::Porosity(MaterialLib::Porosity::Type::Constant, mmp.porosity_model_values[0]);
    }

    if (mmp.storage_model==1) {
        pm.storage = new MaterialLib::SpecificStorage(MaterialLib::SpecificStorage::Type::Constant, mmp.storage_model_values[0]);
    }

    if (mmp.mass_dispersion_model==1) {
        pm.dispersivity_long  = new MaterialLib::MassDispersivityLong(MaterialLib::MassDispersivityLong::Type::Constant, mmp.mass_dispersion_longitudinal);
        pm.dispersivity_trans = new MaterialLib::MassDispersivityTrans(MaterialLib::MassDispersivityTrans::Type::Constant, mmp.mass_dispersion_transverse);
    }

    pm.geo_area = new MaterialLib::GeometricArea(MaterialLib::GeometricArea::Type::Constant, mmp.geo_area);

    if (mmp.heat_capacity_model==1) {
        std::vector<double> para;
        pm.heat_capacity  = new MaterialLib::PorousMediumHeatCapacity(MaterialLib::PorousMediumHeatCapacity::Type::ARITHMETIC, para);
    }

    if (mmp.heat_conductivity_model==1) {
        std::vector<double> para;
        pm.thermal_conductivity = new MaterialLib::PorousMediumThermalConductivity(MaterialLib::PorousMediumThermalConductivity::Type::ARITHMETIC, para);
    }
}


void convertFractureProperty(const CMediumProperties &mmp, MaterialLib::Fracture &pm)
{
    if (mmp.permeability_model==1) {
        pm.hydraulic_conductivity = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
        pm.permeability = new NumLib::TXFunctionConstant(mmp.permeability_tensor[0]);
    }

    if (mmp.porosity_model==1) {
        pm.porosity = new NumLib::TXFunctionConstant(mmp.porosity_model_values[0]);
    }

    if (mmp.storage_model==1) {
        pm.storage = new NumLib::TXFunctionConstant(mmp.storage_model_values[0]);
    }

    pm.geo_area = new NumLib::TXFunctionConstant(mmp.geo_area);
}
void convertCompoundProperty(const CompProperties &mcp, MaterialLib::Compound &new_cp)
{
    new_cp.name = mcp.compname;
    new_cp.is_mobile = (mcp.mobil == 1);

    if (mcp.diffusion_model==1) {
        new_cp.molecular_diffusion = new MaterialLib::MolecularDiffusion(MaterialLib::MolecularDiffusion::Type::Constant, mcp.diffusion_model_values[0]);
    }
}


std::string convertLinearSolverType(int ls_method)
{
    std::string ls_sol_type = "";
    switch (ls_method) {
    case 1:
        ls_sol_type = "GAUSS";
        break;
    case 2:
        ls_sol_type = "BiCGSTAB";
        break;
    case 3:
        ls_sol_type = "BICG";
        break;
    case 5:
        ls_sol_type = "CG";
        break;
    case 805:
        ls_sol_type = "PARDISO";
        break;
    default:
        break;
    } 
    return ls_sol_type;
}

std::string convertLinearSolverPreconType(int ls_precon)
{
    std::string str = "";
    switch (ls_precon) {
    case 0:
        str = "NONE";
        break;
    case 1:
        str = "JACOBI";
        break;
    case 100:
        str = "ILU";
        break;
    default:
        break;
    } 
    return str;
}

NumLib::TXFunctionType::type convertDistributionType(FiniteElement::DistributionType ogs5_type)
{
    switch (ogs5_type) {
        case FiniteElement::CONSTANT:
            return NumLib::TXFunctionType::CONSTANT;
        case FiniteElement::LINEAR:
            return NumLib::TXFunctionType::GEOSPACE;
        case FiniteElement::FUNCTION:
            return NumLib::TXFunctionType::ANALYTICAL;
        default:
            return NumLib::TXFunctionType::INVALID;
    }

}

bool convert(const Ogs5FemData &ogs5fem, THMCLib::Ogs6FemData &ogs6fem, boost::property_tree::ptree & option)
{

    // -------------------------------------------------------------------------
    // Materials
    // -------------------------------------------------------------------------
    // MFP
    const size_t n_fluid = ogs5fem.mfp_vector.size();
    for (size_t i=0; i<n_fluid; i++)
    {
        CFluidProperties* mfp = ogs5fem.mfp_vector[i];

        MaterialLib::FluidModel* fluid = new MaterialLib::FluidModel();
        ogs6fem.list_fluid.push_back(fluid);

        convertFluidProperty(*mfp, *fluid);
    }

    // MSP
    for (size_t i=0; i<ogs5fem.msp_vector.size(); i++)
    {
        CSolidProperties* msp = ogs5fem.msp_vector[i];

        MaterialLib::SolidModel* solid = new MaterialLib::SolidModel();
        ogs6fem.list_solid.push_back(solid);

        convertSolidProperty(*msp, *solid);
    }

    // MMP
    for (size_t i=0; i<ogs5fem.mmp_vector.size(); i++)
    {
        CMediumProperties* mmp = ogs5fem.mmp_vector[i];
        MaterialLib::IMedium* mmp_ogs6 = NULL;
        if (mmp->is_fracture) {
            MaterialLib::Fracture* frac = new MaterialLib::Fracture();
            convertFractureProperty(*mmp, *frac);
            ogs6fem.list_pm.push_back(NULL);
            mmp_ogs6 = frac;
        } else {
            MaterialLib::PorousMediumModel* pm = new MaterialLib::PorousMediumModel();
            convertPorousMediumProperty(*mmp, *pm);
            ogs6fem.list_pm.push_back(pm);
            mmp_ogs6 = pm;
        }
        ogs6fem.list_medium.push_back(mmp_ogs6);
    }

    // MCP
    for (size_t i=0; i<ogs5fem.cp_vector.size(); i++)
    {
        CompProperties* mcp = ogs5fem.cp_vector[i];
        MaterialLib::Compound* new_cp = new MaterialLib::Compound();
        ogs6fem.list_compound.push_back(new_cp);
        convertCompoundProperty(*mcp, *new_cp);
    }

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------
    ogs6fem.geo_unique_name = ogs5fem.geo_unique_name;
    ogs6fem.geo = ogs5fem.geo_obj;

    // -------------------------------------------------------------------------
    // Mesh
    // -------------------------------------------------------------------------
    if (ogs5fem.list_mesh.size()==0) {
        ERR("***Error: no mesh found in ogs5");
        return false;
    }
    for (auto msh : ogs5fem.list_mesh) {
        ogs6fem.list_mesh.push_back(msh);
        ogs6fem.list_nodeSearcher.push_back(new MeshGeoToolsLib::MeshNodeSearcher(*msh));
        ogs6fem.list_beSearcher.push_back(new MeshGeoToolsLib::BoundaryElementsSearcher(*msh, *ogs6fem.list_nodeSearcher.back()));
    }

    // -------------------------------------------------------------------------
    // Time group
    // -------------------------------------------------------------------------
    convertTimeStepping(ogs5fem.time_vector, ogs6fem.list_tim);
    for (std::size_t i=ogs6fem.list_tim.size(); i<ogs5fem.pcs_vector.size(); i++) {
        ogs6fem.list_tim.push_back(ogs6fem.list_tim.front()->clone());
    }

    // -------------------------------------------------------------------------
    // User-defined curves
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Coupling
    // -------------------------------------------------------------------------
    auto &optCoupling = option.put_child("coupling", ptree());
    auto &optP = optCoupling.put_child("P", ptree());
    optP.add("algorithm", "Serial");
    optP.add("convergence", "FemFunctionConvergenceCheck");
    optP.add("max_itr", 1);
    optP.add("epsilon", 1e-4);
    for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
    {
        CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
        std::vector<std::string>& var_name = rfpcs->primary_variable_name;
        for (auto& s : var_name)
            optP.add("out", s);
        if (rfpcs->getProcessType()==FiniteElement::ProcessType::LIQUID_FLOW)
            optP.add("out", "VELOCITY1");
    }
    auto &optProblems = optP.put_child("problems", ptree());
    for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
    {
        CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
        std::string pcs_name = FiniteElement::convertProcessTypeToString(rfpcs->getProcessType());
        std::vector<std::string>& var_name = rfpcs->primary_variable_name;
        auto &optM = optProblems.add_child("M", ptree());
        optM.add("type", pcs_name);
        if (rfpcs->getProcessType()==FiniteElement::ProcessType::HEAT_TRANSPORT)
            optM.add("in", "VELOCITY1");
        for (auto& s : var_name)
            optM.add("out", s);
        if (rfpcs->getProcessType()==FiniteElement::ProcessType::LIQUID_FLOW)
            optM.add("out", "VELOCITY1");
    }


    // -------------------------------------------------------------------------
    // Individual process and IVBV
    // -------------------------------------------------------------------------
    // PCS
    auto optionalPcsData = option.get_child_optional("processList");
    ptree* optPcsData = nullptr;
    if (optionalPcsData)
    	optPcsData = &*optionalPcsData;
    else
        optPcsData = &option.put_child("processList", ptree());
    size_t masstransport_counter = 0;
    if (ogs5fem.pcs_vector.size()==0) {
        ERR("***Error: no PCS found in ogs5");
        return false;
    }

   for (size_t i=0; i<ogs5fem.pcs_vector.size(); i++)
    {
        CRFProcess* rfpcs = ogs5fem.pcs_vector[i];
        std::string pcs_name = FiniteElement::convertProcessTypeToString(rfpcs->getProcessType());
        std::vector<std::string>& var_name = rfpcs->primary_variable_name;

        auto &optPcs = optPcsData->add_child("process", ptree());
        optPcs.add("type", pcs_name);
        optPcs.add("name", pcs_name);

        //Mesh
        optPcs.add("MeshID", rfpcs->mesh_id);

        //Time
        optPcs.add("TimeGroupID", i); //TODO rfpcs->timegroup_id

        // IC
        auto &optIcList = optPcs.put_child("ICList", ptree());
        for (size_t i=0; i<ogs5fem.ic_vector.size(); i++)
        {
            CInitialCondition* rfic = ogs5fem.ic_vector[i];
            std::string ic_pcs_name = FiniteElement::convertProcessTypeToString(rfic->getProcessType());
            if ( ic_pcs_name.compare(pcs_name)==0 ) {
                for (size_t j=0; j<var_name.size(); j++) {
                    if ( rfic->primaryvariable_name.find(var_name[j])!=std::string::npos) {
                        auto &optIc = optIcList.add_child("IC", ptree());
                        optIc.add("Variable", rfic->primaryvariable_name);
                        optIc.add("GeometryType", rfic->geo_type_name);
                        optIc.add("GeometryName", rfic->geo_name);
                        optIc.add("DistributionType", FiniteElement::convertDisTypeToString(rfic->getProcessDistributionType()));
                        optIc.add("DistributionValue", rfic->geo_node_value);
                    }  // end of if rfic
                }
            }
        }  // end of for i

        // BC
        auto& optBcList = optPcs.put_child("BCList", ptree());
        for (size_t i=0; i<ogs5fem.bc_vector.size(); i++)
        {
            CBoundaryCondition* rfbc = ogs5fem.bc_vector[i];
            std::string bc_pcs_name = FiniteElement::convertProcessTypeToString(rfbc->getProcessType());
            if ( bc_pcs_name.compare(pcs_name)==0 )
            {
                for ( size_t j=0; j<var_name.size(); j++ )
                {
                    if ( rfbc->primaryvariable_name.find(var_name[j])!=std::string::npos)
                    {
                        auto &optBc = optBcList.add_child("BC", ptree());
                        optBc.add("Variable", rfbc->primaryvariable_name);
                        optBc.add("GeometryType", rfbc->geo_type_name);
                        optBc.add("GeometryName", rfbc->geo_name);
                        NumLib::TXFunctionType::type ogs6dis_type = convertDistributionType(rfbc->getProcessDistributionType());
                        optBc.add("DistributionType", NumLib::convertTXFunctionTypeToString(ogs6dis_type));
                        switch (rfbc->getProcessDistributionType())
                        {
                            case FiniteElement::CONSTANT:
                                optBc.add("DistributionValue", rfbc->geo_node_value);
                                break;
                            case FiniteElement::LINEAR:
                                {
                                    const size_t n_pt = rfbc->_PointsHaveDistribedBC.size();
                                    optBc.add("PointSize", n_pt);
                                    for (size_t k=0; k<n_pt; k++) {
                                        auto& optPtList = optBc.add_child("PointValueList", ptree());
                                        optPtList.add("PointID", rfbc->_PointsHaveDistribedBC[k]);
                                        optPtList.add("Value", rfbc->_DistribedBC[k]);
                                    }
                                }
                                break;
                            case FiniteElement::FUNCTION:
                                {
                                    optBc.add("DistributionFunction", rfbc->function_exp);
                                }
                                break;
                            default:
                                //error
                                break;
                        }
                    }  // end of if rfbc
                }
            }
        }  // end of for i

        //ST
        auto& optStList = optPcs.put_child("STList", ptree());
        for (size_t i=0; i<ogs5fem.st_vector.size(); i++)
        {
            CSourceTerm* rfst = ogs5fem.st_vector[i];
            std::string st_pcs_name = FiniteElement::convertProcessTypeToString(rfst->getProcessType());
            if (st_pcs_name.compare(pcs_name)==0 )
                for (size_t j=0; j<var_name.size(); j++ ) {
                    if ( rfst->primaryvariable_name.find(var_name[j])!=std::string::npos) {
                        auto& optSt = optStList.add_child("ST", ptree());
                        optSt.add("Variable", rfst->primaryvariable_name);
                        optSt.add("GeometryType", rfst->geo_type_name);
                        optSt.add("GeometryName", rfst->geo_name);
                        std::string ogs5distype = FiniteElement::convertDisTypeToString(rfst->getProcessDistributionType());
                        std::string ogs6sp_distype = ogs5distype;
                        bool isNeumannBC = (ogs5distype.find("_NEUMANN")!=std::string::npos);
                        if (isNeumannBC) {
                            optSt.add("STType", "NEUMANN");
                            ogs6sp_distype = BaseLib::replaceString("_NEUMANN", "", ogs6sp_distype);
                        } else {
                            optSt.add("STType", "SOURCESINK");
                        }
                        bool isTransientST = (rfst->CurveIndex>0);
                        std::string ogs6distype;
                        if (isTransientST) {
                            // spatial
                            auto &opStSpace = optSt.add_child("Spatial", ptree());
                            opStSpace.add("DistributionType", ogs6sp_distype);
                            opStSpace.add("DistributionValue", rfst->geo_node_value);
                            // temporal
                            auto &opStTime = optSt.add_child("Temporal", ptree());
                            std::string ogs6tim_distype = NumLib::convertTXFunctionTypeToString(NumLib::TXFunctionType::T_LINEAR);
                            opStTime.add("DistributionType", ogs6tim_distype);
                            auto &opStTimeValue = opStTime.add_child("TimeValue", ptree());
                            ogs5::Kurven* curv = ogs5fem.kurven_vector[rfst->CurveIndex-1];
                            for (size_t i_pt=0; i_pt<curv->stuetzstellen.size(); i_pt++) {
                                auto& opStTimeValuePoint = opStTimeValue.add_child("Point", ptree());
                                opStTimeValuePoint.add("Time", curv->stuetzstellen[i_pt]->punkt);
                                opStTimeValuePoint.add("Value", curv->stuetzstellen[i_pt]->wert);
                            }

                            ogs6distype = NumLib::convertTXFunctionTypeToString(NumLib::TXFunctionType::TX);
                        } else {
                            ogs6distype = ogs6sp_distype;
                            optSt.add("DistributionValue", rfst->geo_node_value);
                        }
                        optSt.add("DistributionType", ogs6distype);
                    }  // end of if rfst
                }
        }  // end of for i

        // NUM
        auto& optNum = optPcs.put_child("Numerics", ptree());
        for (size_t i=0; i<ogs5fem.num_vector.size(); i++)
        {
            CNumerics* rfnum = ogs5fem.num_vector[i];
            if (rfnum->pcs_type_name.compare(pcs_name)==0) {
                // linear solver
                auto& optLS = optNum.put_child("LinearSolver", ptree());
                optLS.add("solver_type", convertLinearSolverType(rfnum->ls_method));
                optLS.add("precon_type", convertLinearSolverPreconType(rfnum->ls_precond));
                //optLS.addOption("error_type", "NONE");
                optLS.add("error_tolerance", rfnum->ls_error_tolerance);
                optLS.add("max_iteration_step", rfnum->ls_max_iterations);

                // nonlinear solver
                auto &optNS = optNum.put_child("NonlinearSolver", ptree());
                optNS.add("solver_type", rfnum->nls_method_name);
                optNS.add("error_type", rfnum->nls_error_method);
                optNS.add("error_tolerance", rfnum->nls_error_tolerance[0]);
                optNS.add("max_iteration_step", rfnum->nls_max_iterations);

                // other stuff
                optNum.add("TimeTheta", rfnum->ls_theta);
                optNum.add("GaussPoint", rfnum->ele_gauss_points);
                optNum.add("FEM_FCT", rfnum->fct_method);
            }
        }

        // some process specific
        if (pcs_name.find("MASS_TRANSPORT")!=std::string::npos) {
            optPcs.add("CompoundID", masstransport_counter);
            masstransport_counter++;
        }
		
    }  // end of for ogs5fem.pcs_vector.size()


    // -------------------------------------------------------------------------
    // Coupling
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------

    // OUT
    auto& optOut = option.put_child("outputList", ptree());
    for (size_t i=0; i<ogs5fem.out_vector.size(); i++)
    {
        COutput* rfout = ogs5fem.out_vector[i];

        // convert *_X1 to *_X (e.g. DISPLACEMENT_X1 to DISPLACEMENT_X)
        for (size_t j=0; j<rfout->_nod_value_vector.size(); j++) {
            const size_t len_postfix = 3;
            std::string str = rfout->_nod_value_vector[j];
            if (str.length()<len_postfix+1) continue;
            std::string last3 = str.substr(str.length()-len_postfix, len_postfix);
            if (last3.compare("_X1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_X");
            } else if (last3.compare("_Y1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_Y");
            } else if (last3.compare("_Z1")==0) {
                rfout->_nod_value_vector[j].replace(str.length()-len_postfix, len_postfix, "_Z");
            }
        }

        auto& opt = optOut.add_child("output", ptree());
        opt.add("dataType", rfout->dat_type_name);
        opt.add("meshID", "0"); //TODO
        opt.add("geoType", rfout->geo_type);
        opt.add("geoName", rfout->geo_name);
        opt.add("timeType", rfout->tim_type_name);
        if (rfout->tim_type_name == "STEPS") {
            opt.add("timeSteps", rfout->nSteps);
        } else {
            auto& optTimeList = opt.add_child("timeList", ptree());
            for (auto t : rfout->time_vector)
                optTimeList.add("time", t);
        }
        for (size_t j=0; j<rfout->_nod_value_vector.size(); j++) {
            auto& optVal = opt.add_child("nodeValue", ptree());
            optVal.add("name", rfout->_nod_value_vector[j]);
        }
        std::vector<std::string> vecEleVelocity;
        for (auto vname : rfout->_ele_value_vector) {
            const std::string keyword = "VELOCITY";
            if (vname.find(keyword)!=std::string::npos) { //TODO velocity2
                std::string vname2;
                if (vname.size()>keyword.size())
                    vname2 = vname.erase(vname.find(keyword)+8+1, vname.size());
                else
                    vname2 = keyword;
                vecEleVelocity.push_back(vname2);
                continue;
            }
            auto& optVal = opt.add_child("elementValue", ptree());
            optVal.add("name", vname);
        }
        std::sort(vecEleVelocity.begin(), vecEleVelocity.end());
        vecEleVelocity.erase(std::unique(vecEleVelocity.begin(), vecEleVelocity.end()), vecEleVelocity.end());
        for (auto vname : vecEleVelocity) {
            auto& optVal = opt.add_child("elementValue", ptree());
            optVal.add("name", vname);
        }
//        opt.addOptionAsArray("MMPValues", rfout->mmp_value_vector);
//        opt.addOptionAsArray("MFPValues", rfout->mfp_value_vector);
    }

    ogs6fem.outController.initialize(option, ogs6fem.output_dir, ogs6fem.project_name,
            ogs6fem.list_mesh, ogs6fem.list_nodeSearcher, *ogs6fem.geo, ogs6fem.geo_unique_name);

    return true;
}

};
} //end
