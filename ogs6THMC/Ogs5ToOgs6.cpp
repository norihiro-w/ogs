/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5ToOgs6.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "Ogs5ToOgs6.h"

#include "logog.hpp"

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "MaterialLib/Fluid.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "GeoProcessBuilder.h"
#include "ChemLib/chemReactionKin.h"
#include "ChemLib/chemReductionKin.h"

using namespace ogs5;

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertTimeStepping(const std::vector<CTimeDiscretization*> &td_vector, std::vector<NumLib::ITimeStepFunction*> &tf_vector)
{
    bool done_shared = false;
    for (size_t i=0; i<td_vector.size(); i++)
    {
        const CTimeDiscretization* td = td_vector[i];
        if (td->time_independence || !done_shared) {
            NumLib::ITimeStepFunction* tf = new NumLib::TimeStepFunctionVector(td->time_start, td->time_end, td->time_step_vector);
            tf_vector.push_back(tf);
            if (!td->time_independence)
                done_shared = true;
        }
    }
}

void convertFluidProperty(const CFluidProperties &mfp, MaterialLib::Fluid &fluid)
{
    if (mfp.density_model==1) {
        fluid.density = new NumLib::TXFunctionConstant(mfp.rho_0);
    }

    if (mfp.viscosity_model==1) {
        fluid.dynamic_viscosity = new NumLib::TXFunctionConstant(mfp.my_0);
    }

    if (mfp.heat_capacity_model==1) {
        fluid.specific_heat = new NumLib::TXFunctionConstant(mfp.specific_heat_capacity);
    }

    if (mfp.heat_conductivity_model==1) {
        fluid.thermal_conductivity = new NumLib::TXFunctionConstant(mfp.heat_conductivity);
    }
}

void convertSolidProperty(const CSolidProperties &msp, MaterialLib::Solid &solid)
{
    if (msp.Density_mode==1) {
        solid.density = new NumLib::TXFunctionConstant((*msp.data_Density)(0));
    }

    solid.poisson_ratio = new NumLib::TXFunctionConstant(msp.PoissonRatio);

    if (msp.Youngs_mode==1) {
        solid.Youngs_modulus = new NumLib::TXFunctionConstant((*msp.data_Youngs)(0));
    }

    if (msp.Capacity_mode==1) {
        solid.specific_heat = new NumLib::TXFunctionConstant((*msp.data_Capacity)(0));
    }

    if (msp.Conductivity_mode==1) {
        solid.thermal_conductivity = new NumLib::TXFunctionConstant((*msp.data_Conductivity)(0));
    }
}

void convertPorousMediumProperty(const CMediumProperties &mmp, MaterialLib::PorousMedia &pm)
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

    if (mmp.mass_dispersion_model==1) {
        pm.dispersivity_long  = new NumLib::TXFunctionConstant(mmp.mass_dispersion_longitudinal); 
        pm.dispersivity_trans = new NumLib::TXFunctionConstant(mmp.mass_dispersion_transverse);
    }

    pm.geo_area = new NumLib::TXFunctionConstant(mmp.geo_area);
}

void convertCompoundProperty(const CompProperties &mcp, MaterialLib::Compound &new_cp)
{
    new_cp.name = mcp.compname;
    new_cp.is_mobile = (mcp.mobil == 1);

    if (mcp.diffusion_model==1) {
        new_cp.molecular_diffusion = new NumLib::TXFunctionConstant(mcp.diffusion_model_values[0]);
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
        ls_sol_type = "BICGSTAB";
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

bool convert(const Ogs5FemData &ogs5fem, Ogs6FemData &ogs6fem, BaseLib::Options &option)
{

    // -------------------------------------------------------------------------
    // Materials
    // -------------------------------------------------------------------------
    // MFP
    const size_t n_fluid = ogs5fem.mfp_vector.size();
    for (size_t i=0; i<n_fluid; i++)
    {
        CFluidProperties* mfp = ogs5fem.mfp_vector[i];

        MaterialLib::Fluid* fluid = new MaterialLib::Fluid();
        ogs6fem.list_fluid.push_back(fluid);

        convertFluidProperty(*mfp, *fluid);
    }

    // MSP
    for (size_t i=0; i<ogs5fem.msp_vector.size(); i++)
    {
        CSolidProperties* msp = ogs5fem.msp_vector[i];

        MaterialLib::Solid* solid = new MaterialLib::Solid();
        ogs6fem.list_solid.push_back(solid);

        convertSolidProperty(*msp, *solid);
    }

    // MMP
    for (size_t i=0; i<ogs5fem.mmp_vector.size(); i++)
    {
        CMediumProperties* mmp = ogs5fem.mmp_vector[i];
        MaterialLib::PorousMedia* pm = new MaterialLib::PorousMedia();
        ogs6fem.list_pm.push_back(pm);
        convertPorousMediumProperty(*mmp, *pm);
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
    for (size_t i=0; i<ogs5fem.list_mesh.size(); i++) {
        ogs6fem.list_mesh.push_back(ogs5fem.list_mesh[i]);
//        ogs6fem.list_dis_sys.push_back(new DiscreteLib::DiscreteSystem(*ogs5fem.list_mesh[i]));
    }

    // -------------------------------------------------------------------------
    // Time group
    // -------------------------------------------------------------------------
    convertTimeStepping(ogs5fem.time_vector, ogs6fem.list_tim);

    // -------------------------------------------------------------------------
    // User-defined curves
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Individual process and IVBV
    // -------------------------------------------------------------------------
    // PCS
    BaseLib::Options* optPcsData = option.getSubGroup("processList");
    if (optPcsData==NULL)
        optPcsData = option.addSubGroup("processList");
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

        BaseLib::Options* optPcs = optPcsData->addSubGroup("process");
        optPcs->addOption("type", pcs_name);
        optPcs->addOption("name", pcs_name);

        //Mesh
        optPcs->addOptionAsNum("MeshID", rfpcs->mesh_id);

        //Time
        optPcs->addOptionAsNum("TimeGroupID", rfpcs->timegroup_id);

        // IC
        BaseLib::Options* optIcList = optPcs->addSubGroup("ICList");
        for (size_t i=0; i<ogs5fem.ic_vector.size(); i++)
        {
            CInitialCondition* rfic = ogs5fem.ic_vector[i];
            std::string ic_pcs_name = FiniteElement::convertProcessTypeToString(rfic->getProcessType());
			if ( ic_pcs_name.compare(pcs_name)==0 ) {
				for (size_t j=0; j<var_name.size(); j++) {
					if ( rfic->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optIc = optIcList->addSubGroup("IC");
						optIc->addOption("Variable", rfic->primaryvariable_name);
						optIc->addOption("GeometryType", rfic->geo_type_name);
						optIc->addOption("GeometryName", rfic->geo_name);
						optIc->addOption("DistributionType", FiniteElement::convertDisTypeToString(rfic->getProcessDistributionType()));
						optIc->addOptionAsNum("DistributionValue", rfic->geo_node_value);
					}  // end of if rfic
				}
			}
        }  // end of for i

        // BC
        BaseLib::Options* optBcList = optPcs->addSubGroup("BCList");
        for (size_t i=0; i<ogs5fem.bc_vector.size(); i++)
        {
            CBoundaryCondition* rfbc = ogs5fem.bc_vector[i];
            std::string bc_pcs_name = FiniteElement::convertProcessTypeToString(rfbc->getProcessType());
            if ( bc_pcs_name.compare(pcs_name)==0 )
				for ( size_t j=0; j<var_name.size(); j++ )
					if ( rfbc->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optBc = optBcList->addSubGroup("BC");
						optBc->addOption("Variable", rfbc->primaryvariable_name);
						optBc->addOption("GeometryType", rfbc->geo_type_name);
						optBc->addOption("GeometryName", rfbc->geo_name);
						NumLib::TXFunctionType::type ogs6dis_type = convertDistributionType(rfbc->getProcessDistributionType());
                        optBc->addOption("DistributionType", NumLib::convertTXFunctionTypeToString(ogs6dis_type));
						switch (rfbc->getProcessDistributionType())
						{
                            case FiniteElement::CONSTANT:
                                optBc->addOptionAsNum("DistributionValue", rfbc->geo_node_value);
                                break;
                            case FiniteElement::LINEAR:
                                {
                                    const size_t n_pt = rfbc->_PointsHaveDistribedBC.size();
                                    optBc->addOptionAsNum("PointSize", n_pt);
                                    for (size_t k=0; k<n_pt; k++) {
                                        BaseLib::Options* optPtList = optBc->addSubGroup("PointValueList");
                                        optPtList->addOptionAsNum("PointID", rfbc->_PointsHaveDistribedBC[k]);
                                        optPtList->addOptionAsNum("Value", rfbc->_DistribedBC[k]);
                                    }
                                }
                                break;
                            case FiniteElement::FUNCTION:
                                {
                                    optBc->addOption("DistributionFunction", rfbc->function_exp);
                                }
                                break;
                            default:
                                //error
                                break;
						}
					}  // end of if rfbc
        }  // end of for i

        //ST
        BaseLib::Options* optStList = optPcs->addSubGroup("STList");
        for (size_t i=0; i<ogs5fem.st_vector.size(); i++)
        {
            CSourceTerm* rfst = ogs5fem.st_vector[i];
            std::string st_pcs_name = FiniteElement::convertProcessTypeToString(rfst->getProcessType());
            if (st_pcs_name.compare(pcs_name)==0 )
				for (size_t j=0; j<var_name.size(); j++ )
					if ( rfst->primaryvariable_name.find(var_name[j])!=std::string::npos) {
						BaseLib::Options* optSt = optStList->addSubGroup("ST");
						optSt->addOption("Variable", rfst->primaryvariable_name);
						optSt->addOption("GeometryType", rfst->geo_type_name);
						optSt->addOption("GeometryName", rfst->geo_name);
						std::string ogs5distype = FiniteElement::convertDisTypeToString(rfst->getProcessDistributionType());
						std::string ogs6distype = ogs5distype;
						if (ogs5distype.find("_NEUMANN")!=std::string::npos) {
							optSt->addOption("STType", "NEUMANN");
							ogs6distype = BaseLib::replaceString("_NEUMANN", "", ogs6distype);
						} else {
							optSt->addOption("STType", "SOURCESINK");
						}  // end of if-else
						optSt->addOption("DistributionType", ogs6distype);
						optSt->addOptionAsNum("DistributionValue", rfst->geo_node_value);
					}  // end of if rfst
		}  // end of for i

        // NUM
        BaseLib::Options* optNum = optPcs->addSubGroup("Numerics");
        for (size_t i=0; i<ogs5fem.num_vector.size(); i++)
        {
            CNumerics* rfnum = ogs5fem.num_vector[i];
            if (rfnum->pcs_type_name.compare(pcs_name)==0) {
                // linear solver
                BaseLib::Options* optLS = optNum->addSubGroup("LinearSolver");
                optLS->addOption("solver_type", convertLinearSolverType(rfnum->ls_method));
                optLS->addOption("precon_type", convertLinearSolverPreconType(rfnum->ls_precond));
                //optLS->addOption("error_type", "NONE");
                optLS->addOptionAsNum("error_tolerance", rfnum->ls_error_tolerance);
                optLS->addOptionAsNum("max_iteration_step", rfnum->ls_max_iterations);

                // nonlinear solver
                BaseLib::Options* optNS = optNum->addSubGroup("NonlinearSolver");
                optNS->addOption("solver_type", rfnum->nls_method_name);
                optNS->addOptionAsNum("error_type", rfnum->nls_error_method);
                optNS->addOptionAsNum("error_tolerance", rfnum->nls_error_tolerance[0]);
                optNS->addOptionAsNum("max_iteration_step", rfnum->nls_max_iterations);

                // other stuff
                optNum->addOptionAsNum("TimeTheta", rfnum->ls_theta);
                optNum->addOptionAsNum("GaussPoint", rfnum->ele_gauss_points);
                optNum->addOptionAsNum("FEM_FCT", rfnum->fct_method);
            }
        }

        // some process specific
        if (pcs_name.find("MASS_TRANSPORT")!=std::string::npos) {
            optPcs->addOptionAsNum("CompoundID", masstransport_counter);
            masstransport_counter++;
        }
		
    }  // end of for ogs5fem.pcs_vector.size()

	// HS: add kinetic reactions
	size_t n_KinReactions, n_Compound; 
	n_KinReactions = ogs5fem.KinReact_vector.size(); 
	if ( n_KinReactions > 0 ) {  // only create the "optKRC" option when there are kinetic reactions defined. 
		// HS: convert all compounds to ChemComp data structure
		n_Compound = ogs6fem.list_compound.size(); 
		for ( size_t i=0; i < n_Compound ; i++ )
		{
			ogsChem::ChemComp* mChemComp = NULL; 
			// directly use the type of constructor to initialize the ChemComp data structure
			mChemComp = new ogsChem::ChemComp(); 
			// set index
			mChemComp->setIndex(i); 
			// set name
			mChemComp->set_name(ogs6fem.list_compound[i]->name); 
			// set mobility
			if ( ogs6fem.list_compound[i]->is_mobile )
				mChemComp->set_mobility( ogsChem::MOBILE  ); 
			else 
				mChemComp->set_mobility( ogsChem::MINERAL );

			ogs6fem.map_ChemComp.insert(ogs6fem.list_compound[i]->name, mChemComp); 
			mChemComp = NULL; 
		}

		// adding the kinetic reactions one after another
		for ( size_t i=0; i < n_KinReactions ; i++ )
		{
			// add reactions one after another
			ogs5::CKinReact* rfKinReact = ogs5fem.KinReact_vector[i];
			// creating reaction instance 
			ogsChem::chemReactionKin* mKinReaction; 
			mKinReaction = new ogsChem::chemReactionKin(); 
			// convert the ogs5 KRC data structure into ogs6 Kinetic reactions
			mKinReaction->readReactionKRC( ogs6fem.map_ChemComp, rfKinReact ); 
			// adding the instance of one single kinetic reaction
			ogs6fem.list_kin_reactions.push_back(mKinReaction); 
		}  // end of for

		// at last, initialize the reduction scheme 
		ogsChem::chemReductionKin* mReductionScheme; 
		mReductionScheme = new ogsChem::chemReductionKin(ogs6fem.map_ChemComp, ogs6fem.list_kin_reactions); 
		ogs6fem.m_KinReductScheme = mReductionScheme; 
		mReductionScheme = NULL;
	}  // end of if ( n_KinReactions > 0 )

    // -------------------------------------------------------------------------
    // Coupling
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------

    // OUT
    BaseLib::Options* optOut = option.addSubGroup("outputList");
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

        BaseLib::Options* opt = optOut->addSubGroup("output");
        opt->addOption("dataType", rfout->dat_type_name);
        opt->addOption("meshID", "0"); //TODO
        opt->addOption("geoType", rfout->geo_type);
        opt->addOption("geoName", rfout->geo_name);
        opt->addOption("timeType", rfout->tim_type_name);
        opt->addOptionAsNum("timeSteps", rfout->nSteps);
        for (size_t j=0; j<rfout->_nod_value_vector.size(); j++) {
            BaseLib::Options* optVal = opt->addSubGroup("nodeValue");
            optVal->addOption("name", rfout->_nod_value_vector[j]);
        }
        for (size_t j=0; j<rfout->_ele_value_vector.size(); j++) {
            BaseLib::Options* optVal = opt->addSubGroup("elementValue");
            optVal->addOption("name", rfout->_ele_value_vector[j]);
        }
//        opt->addOptionAsArray("MMPValues", rfout->mmp_value_vector);
//        opt->addOptionAsArray("MFPValues", rfout->mfp_value_vector);
    }

    ogs6fem.outController.initialize(option, ogs6fem.output_dir, ogs6fem.project_name, ogs6fem.list_mesh, *ogs6fem.geo, ogs6fem.geo_unique_name);

    return true;
}

};
} //end
