/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearTransportTimeODELocalAssembler.h
 *
 * Created on 2012-09-19 by Haibing Shao
 */

#ifndef LINEAR_TRANSPORT_TIME_ODE_LOCAL_ASSEMBLER_H
#define LINEAR_TRANSPORT_TIME_ODE_LOCAL_ASSEMBLER_H

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Integration/Integration.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TransientAssembler/IElementWiseTimeODELocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Compound.h"

#include "Ogs6FemData.h"

/**
 * \brief Local assembly of time ODE components for linear transport in porous media
 * This class is same as the MassTransportTimeODELocalAssembler class 
 * onlye difference is that no compound information is provided and no molecular 
 * diffusion is included in the assembly. 
 */
template <class T>
class LinearTransportTimeODELocalAssembler: public T
{
public:
    typedef MathLib::LocalVector LocalVectorType;
    typedef MathLib::LocalMatrix LocalMatrixType;

    /**
      * constructor
      */ 
    LinearTransportTimeODELocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects)
        : _feObjects(*feObjects), _vel(NULL)
    {
    };

    /**
      * destructor
      */ 
    virtual ~LinearTransportTimeODELocalAssembler() 
    {
    };

    /**
      * set the velocity
      */ 
    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

protected:

    /**
      * assemble local Mass, Advection and Dispersion matrix
      */ 
    virtual void assembleODE(const NumLib::TimeStep &/*time*/, const MeshLib::IElement &e, const LocalVectorType &/*u1*/, const LocalVectorType &/*u0*/, LocalMatrixType &localM, LocalMatrixType &localK, LocalVectorType &/*localF*/)
    {
        FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
        const size_t n_dim = e.getDimension();
        size_t mat_id = e.getGroupID();
        MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        LocalMatrixType localDispersion(localK);
        LocalMatrixType localAdvection(localK);
        double cmp_mol_diffusion = .0;
        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        MathLib::LocalMatrix poro(1,1);
        MathLib::LocalMatrix d_poro = MathLib::LocalMatrix::Zero(3,3);   // HS, change size to 1 by 3
        double disp_l = 0.0; 
        double disp_t = 0.0; 
        NumLib::ITXFunction::DataType v;

        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
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

            // calculating dispersion tensor according to Benchmark book p219, Eq. 10.15
            // D_{ij} = \alpha_T |v| \delta_{ij} + (\alpha_L - \alpha_T) \frac{v_i v_j}{|v|} + D^{d}_{ii} 
            NumLib::ITXFunction::DataType dispersion_diffusion = NumLib::ITXFunction::DataType::Identity(n_dim, n_dim); 
            dispersion_diffusion *= disp_l * v.norm(); 
            dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm(); 
            dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);
            // --------debugging--------------
            // std::cout << "dispersion_diffusion Matrix" << std::endl;
            // std::cout << dispersion_diffusion << std::endl;
            // --------end of debugging-------

            fe->integrateWxN(j, poro, localM);
            fe->integrateDWxDN(j, dispersion_diffusion, localDispersion); 
            fe->integrateWxDN(j, v2, localAdvection);
        }

        localK = localDispersion + localAdvection;
    }

private:
    /**
      * FEM object
      */ 
    FemLib::LagrangianFeObjectContainer _feObjects;
    
    /**
      * velocity function
      */ 
    NumLib::ITXFunction* _vel;
};



#endif  // end of ifndef