/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonLinearReactiveTransportTimeODELocalAssembler.h
 *
 * Created on 2012-09-24 by Haibing Shao
 */

#ifndef NON_LINEAR_REACTIVE_TRANSPORT_TIME_ODE_LOCAL_ASSEMBLER_H
#define NON_LINEAR_REACTIVE_TRANSPORT_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"
#include "ChemLib/chemReductionKin.h"
#include "Concentrations.h"

#include "Ogs6FemData.h"

/**
 * \brief Local assembly of time ODE components for linear transport in porous media
 * This class is same as the MassTransportTimeODELocalAssembler class 
 * onlye difference is that no compound information is provided and no molecular 
 * diffusion is included in the assembly. 
 */
template <class T, class T_NODAL_FUNCTION_SCALAR>
class NonLinearReactiveTransportTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    NonLinearReactiveTransportTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects, ogsChem::chemReductionKin* ReductionScheme)
        : _feObjects(*feObjects), _vel(NULL), _reductionKin(ReductionScheme), _xi_mob_rates(NULL), _xi_immob_rates(NULL)
    {
    };

    virtual ~NonLinearReactiveTransportTimeODELocalAssembler() {};

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

protected:
    virtual void assembleODE(const NumLib::TimeStep & time, const MeshLib::IElement &e, const LocalVectorType & u1, const LocalVectorType & u0, LocalMatrixType & localM, LocalMatrixType & localK, LocalVectorType & localF)
    {
		size_t i, j, k, node_idx, n_xi_mob, n_nodes, n_sp, n_rows_localK, n_cols_localK;
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		n_rows_localK = localK.rows(); 
		n_cols_localK = localK.cols(); 
		LocalMatrixType localDispersion = LocalMatrixType::Zero( n_rows_localK, n_cols_localK ); 
        LocalMatrixType localAdvection  = LocalMatrixType::Zero( n_rows_localK, n_cols_localK ); 
		LocalMatrixType rate_xi_mob_gp  = LocalMatrixType::Zero( 1            , 1             );  // HS

		double cmp_mol_diffusion = .0;
        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);  // HS

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro = MathLib::LocalMatrix::Zero(3,3);   // HS, change size to 1 by 3
        double disp_l = 0.0; 
        double disp_t = 0.0; 
        NumLib::ITXFunction::DataType v;

		// getting node xi_mob rate values
		n_xi_mob = _xi_mob_rates->size(); 
		n_nodes  = e.getNumberOfNodes(); 
        n_sp     = q->getNumberOfSamplingPoints();  // number of sampling points
		LocalMatrixType node_xi_mob_rate_values = LocalMatrixType::Zero(n_nodes, n_xi_mob ); // TODO
		for (i=0; i<n_nodes; i++)
		{
		    node_idx = e.getNodeID( i ); 
			for (k=0; k<n_xi_mob; k++)
			    node_xi_mob_rate_values( i, k ) = _xi_mob_rates->at(k)->getValue( node_idx ); 
		} // end of for i

        LocalMatrixType localDispersion_tmp = LocalMatrixType::Zero(n_nodes, n_nodes); 
        LocalMatrixType localAdvection_tmp  = LocalMatrixType::Zero(n_nodes, n_nodes); 
        LocalMatrixType localM_tmp          = LocalMatrixType::Zero(n_nodes, n_nodes); 
        LocalMatrixType localK_tmp          = LocalMatrixType::Zero(n_nodes, n_nodes); 

        for (j=0; j<n_sp; j++) 
        {
			q->getSamplingPoint(j, gp_x);
			fe->computeBasisFunctions(gp_x);
			fe->getRealCoordinates(real_x);
			NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e.getID(), j, real_x);

			pm->porosity->eval(gp_pos, poro);
            pm->dispersivity_long->eval(gp_pos, disp_l);
            pm->dispersivity_trans->eval(gp_pos, disp_t);
            
            d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
            d_poro(1,1) = cmp_mol_diffusion * poro(0,0);
            d_poro(2,2) = cmp_mol_diffusion * poro(0,0);

			_vel->eval(gp_pos, v);
			NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

            NumLib::ITXFunction::DataType dispersion_diffusion = NumLib::ITXFunction::DataType::Identity(n_dim, n_dim); 
            dispersion_diffusion *= disp_l * v.norm(); 
            dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm(); 
            dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);

    		// mass matrix
			fe->integrateWxN(j, poro, localM_tmp);
			// dispersion
			fe->integrateDWxDN(j, dispersion_diffusion, localDispersion_tmp);
			// advection
			fe->integrateWxDN(j, v2, localAdvection_tmp);

            LocalMatrixType &Np = *fe->getBasisFunction(); // HS
         	for (k=0; k<n_xi_mob; k++)
            {
                // localF 
                rate_xi_mob_gp = Np * node_xi_mob_rate_values.col(k); 
                // right hand side xi_mob rates
                localF.segment(n_nodes*k,n_nodes).noalias() += Np.transpose() * rate_xi_mob_gp * fe->getDetJ() * q->getWeight(j);
            }
	    }  // end of for j
        localK_tmp = localDispersion_tmp + localAdvection_tmp; 

		for (k=0; k<n_xi_mob; k++)
        {
            // localM
            localM.block(n_nodes*k,n_nodes*k,n_nodes,n_nodes) = localM_tmp; 
            // localK
            localK.block(n_nodes*k,n_nodes*k,n_nodes,n_nodes) = localK_tmp; 
        }  // end of for k

    }

private:
    FemLib::LagrangianFeObjectContainer _feObjects;
    NumLib::ITXFunction* _vel;
	ogsChem::chemReductionKin* _reductionKin; 

	std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_mob_rates;
    std::vector<T_NODAL_FUNCTION_SCALAR*> * _xi_immob_rates;
};



#endif  // end of ifndef