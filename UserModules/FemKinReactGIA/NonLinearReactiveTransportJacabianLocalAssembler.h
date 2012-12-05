/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearReactiveTransportJacobianLocalAssembler.h
 *
 * Created on 2012-09-24 by Haibing Shao
 */

/**
  * This file is same as the MassTransportTimeODELocalAssembler.h
  * The difference is, the compound molecular diffusion coefficient is disabled, 
  */

#ifndef NON_LINEAR_REACTIVE_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H
#define NON_LINEAR_REACTIVE_TRANSPORT_JACOBIAN_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "ChemLib/chemReductionKin.h"
#include "Concentrations.h"
#include "Ogs6FemData.h"

template <class T_NODAL_FUNCTION_SCALAR, class T_FUNCTION_DATA>
class NonLinearReactiveTransportJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
	NonLinearReactiveTransportJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, ogsChem::chemReductionKin* ReductionScheme, T_FUNCTION_DATA* concentrations)
        : _feObjects(*feObjects), _vel(NULL), _reductionKin(ReductionScheme), _concentrations(concentrations), _xi_mob_rates(NULL), _xi_immob_rates(NULL)
    {
    };

    virtual ~NonLinearReactiveTransportJacobianLocalAssembler()
    {
        _vel              = NULL;
        _reductionKin     = NULL; 
        _concentrations   = NULL;
        _xi_mob_rates     = NULL;
        _xi_mob_rates_old = NULL;
        _xi_immob_rates   = NULL;
        _drates_dxi       = NULL;
    };

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

	void set_xi_mob_rates(   std::vector<T_NODAL_FUNCTION_SCALAR*> * xi_mob_rates )
	{
	    _xi_mob_rates = xi_mob_rates; 
	}
    
	void set_xi_immob_rates( std::vector<T_NODAL_FUNCTION_SCALAR*> * xi_immob_rates )
	{
	    _xi_immob_rates = xi_immob_rates; 
	}

	void set_xi_mob_drate_dxi( std::vector<T_NODAL_FUNCTION_SCALAR*> * drates_dxi )
	{
	    _drates_dxi = drates_dxi; 
	}

	T_FUNCTION_DATA* get_function_data(void)
	{
	    return _concentrations; 
	}

	void assembly(const NumLib::TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &, const MathLib::LocalVector & u1, const MathLib::LocalVector & u0, MathLib::LocalMatrix & localJ)
    {
		size_t i, j, k, m, node_idx; 
  		const size_t n_nodes = e.getNumberOfNodes(); 
		const size_t n_xi_mob = _xi_mob_rates->size(); 
        const size_t mat_id  = e.getGroupID(); 
		const size_t n_dim = e.getDimension();
 
        // clear the local Jacobian matrix
        localJ = MathLib::LocalMatrix::Zero(n_nodes*n_xi_mob, n_nodes*n_xi_mob);

		FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
        double cmp_mol_diffusion = .0;
        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);
        double dt = time.getTimeStepSize();  // time step size
        double theta = 1.0;

        MathLib::LocalMatrix matM(n_nodes, n_nodes);
        MathLib::LocalMatrix matDiff(n_nodes, n_nodes);
        MathLib::LocalMatrix matAdv(n_nodes, n_nodes);
        MathLib::LocalMatrix mat_dR(n_nodes, n_nodes); 
		MathLib::LocalVector vec_drates_dxi_value = MathLib::LocalVector::Zero(n_nodes);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();

		size_t n_sp    = q->getNumberOfSamplingPoints();  // number of sampling points
		double gp_x[3], real_x[3];
        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro = MathLib::LocalMatrix::Zero(3,3);   // HS, change size to 1 by 3
        double disp_l = 0.0; 
        double disp_t = 0.0; 
        MathLib::LocalMatrix d_rate = MathLib::LocalMatrix::Zero(1,1); 
        MathLib::LocalMatrix localJ_tmp = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
        NumLib::ITXFunction::DataType v;
		NumLib::ITXFunction::DataType v2;
		
        for (i=0; i<n_xi_mob; i++) {

			matM    = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
			matDiff = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);
			matAdv  = MathLib::LocalMatrix::Zero(n_nodes, n_nodes);

			for (j=0; j<n_sp; j++) {
				q->getSamplingPoint(j, gp_x);
				fe->computeBasisFunctions(gp_x);
				fe->getRealCoordinates(real_x);
                NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);

				pm->porosity->eval(real_x, poro);
                pm->dispersivity_long->eval(gp_pos, disp_l);
                pm->dispersivity_trans->eval(gp_pos, disp_t);

                d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
                d_poro(1,1) = cmp_mol_diffusion * poro(0,0);
                d_poro(2,2) = cmp_mol_diffusion * poro(0,0);

				_vel->eval(real_x, v);
				v2 = v.topRows(n_dim).transpose();
                
                NumLib::ITXFunction::DataType dispersion_diffusion = NumLib::ITXFunction::DataType::Identity(n_dim, n_dim); 
                dispersion_diffusion *= disp_l * v.norm(); 
                dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm(); 
                dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);
				
                fe->integrateWxN(j, poro, matM);
				fe->integrateDWxDN(j, dispersion_diffusion, matDiff);
				fe->integrateWxDN(j, v2, matAdv); 
            
                MathLib::LocalMatrix &Np = *fe->getBasisFunction(); // HS
				// now dealing with the rate change terms
				// each location has n_xi_mob * n_xi_mob dR/dxi entries
				for (m=0; m<n_xi_mob; m++)  // this is to which derivative term
				{
					mat_dR = MathLib::LocalMatrix::Zero(n_nodes, n_nodes); 
					// loop over all the adjacent nodes to get the right drate/dxi value
           			for (k=0; k<n_nodes; k++)
           			{
                        node_idx = e.getNodeID( k ); 
                        vec_drates_dxi_value( k ) = _drates_dxi->at(i*n_xi_mob+m)->getValue( node_idx ); 
                    }  // end of for k<n_nodes
                    
                    // get mean value
                    // HS, disable now and try using basis function -----------
                    // d_rate(0,0) = vec_drates_dxi_value.mean(); 
                    // --------------------------------------------------------
                    d_rate = Np * vec_drates_dxi_value; 
					fe->integrateWxN(j, d_rate, mat_dR);
					
                    // plugging mat_dR to the corresponding Jacobian matrix position
					localJ.block(n_nodes*i, n_nodes*m, n_nodes, n_nodes) -= mat_dR;
                    
				}  // end of for m
			}  // end of for j (sp)

			matM *= 1.0 / dt;
			matDiff *= theta;
			matAdv *= theta;
			localJ_tmp = matM;
			localJ_tmp += matDiff;
			localJ_tmp += matAdv;
            
            localJ.block(n_nodes*i,n_nodes*i,n_nodes,n_nodes) += localJ_tmp;

        }  // end of for i

    }  // end of function assembly

private:
    /**
      * FEM object
      */ 
    FemLib::LagrangianFeObjectContainer _feObjects;

    /**
      * velocity function
      */ 
    NumLib::ITXFunction* _vel;

    /**
      * pointer to the reduction scheme
      */ 
	ogsChem::chemReductionKin* _reductionKin; 

    /**
      * concentration data
      */ 
	T_FUNCTION_DATA* _concentrations; 

    /**
      * nodal xi_mob rate values
      */ 
	std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_mob_rates;

    /**
      * nodal xi_mob rate old values
      */ 
    std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_mob_rates_old;
    
    /**
      * nodal xi_immob rate values
      */ 
    std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_immob_rates;
	
    /**
      * nodal drates over dxi values
      */ 
    std::vector<T_NODAL_FUNCTION_SCALAR*> * _drates_dxi;
};

#endif  // end of ifndef
