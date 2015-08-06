/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <valarray>
#include <cmath>

#include "MathLib/TemplateVectorX.h"
#include "MathLib/DataType.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

#include "MeshLib/Mesh.h"

#include "NumLib/Function/ITXDiscreteFunction.h"



namespace NumLib
{
/**
 * \brief Template class for FEM integration point-based functions
 */
template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public ITXDiscreteFunction //: public NumLib::ITXDiscreteFunction<MathLib::TemplateVectorX<Tvalue> >
{
public:
    typedef MathLib::TemplateVectorX<Tvalue> IntegrationPointVectorType;
    typedef std::valarray<IntegrationPointVectorType> VectorType;
    typedef TemplateFEMIntegrationPointFunction<Tvalue> MyClassType;

    TemplateFEMIntegrationPointFunction()
    : _msh(nullptr), _values(nullptr)
    {};

    virtual ~TemplateFEMIntegrationPointFunction()
    {
        delete _values;
    }

//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len)
//    {
//        initialize(dis, len);
//    };
//
//    ///
//    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len, Tvalue &v0)
//    {
//        initialize(dis, len);
//        for (size_t i=0; i<len; i++) {
//            (*_values)[i].resize(1);
//            (*_values)[i][0] = v0;
//        }
//    };

    ///
    explicit TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
    : _msh(nullptr), _values(nullptr)
    {
        initialize(src._msh, src._values);
        //TODO (*this->_values) = (*src._values);
    };

    /// vector operation: set data
    virtual TemplateFEMIntegrationPointFunction& operator= (const TemplateFEMIntegrationPointFunction &src)
    {
        initialize(src._msh, src._values);
        return *this;
    }

    ///
    TemplateFEMIntegrationPointFunction<Tvalue>* clone() const
    {
        return new TemplateFEMIntegrationPointFunction<Tvalue>(*this);
    };

//    ///
//    const MeshLib::Mesh* getMesh() const
//    {
//        return this->_discrete_system->getMesh();
//    }

    ///
    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::IntegrationPoint:
            {
                size_t ele_id = x.getId(0);
                size_t gp_id = x.getId(1);
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()<gp_id+1) return;
#ifdef OGS_USE_EIGEN
                v = gp_values[gp_id]; 
#endif
            }
            break;
        case NumLib::TXPosition::Element:
            {
                // calculate mean value
                size_t ele_id = x.getId();
                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
                if (gp_values.size()==0) return;
                Tvalue val = gp_values[0]; 
                val *= .0;
                for (size_t i=0; i<gp_values.size(); i++)
                    val += gp_values[i];
                val /= gp_values.size();
#ifdef OGS_USE_EIGEN
                v = val;
#else
#endif
            }
            break;
        default:
            if (_values->size() == 0) return;
            if ((*_values)[0].size() == 0) return;
#ifdef OGS_USE_EIGEN
            v = (*_values)[0][0];
#else
#endif
            break;
        }
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q )
    {
        assert(ip<(*_values)[i_e].size());
        (*_values)[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n)
    {
        (*_values)[i_e].resize(n);
    }

    const IntegrationPointVectorType& getIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e];
    }

    bool hasIntegrationPointValues(size_t i_e) const
    {
        return (*_values)[i_e].size()>0;
    }

    const VectorType* getDiscreteData() const
    {
        return _values;
    }

    VectorType* getDiscreteData()
    {
        return _values;
    }

    void printout() const
    {
        std::cout << "integration_pt_values = ";
        for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
            const IntegrationPointVectorType &val1 = (*_values)[i];
            std::cout << "(";
            for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    virtual double diff_norm(const ITXDiscreteFunction &v, MathLib::VecNormType normType) const
    {
        const MyClassType &vec(static_cast<const MyClassType&>(v));
        if (!_values) return 0;
        double norm = 0.0;
        if (normType==MathLib::VecNormType::INFINITY_N)
        {
            for (size_t i=0; i<_values->size(); i++) {
                const IntegrationPointVectorType &val1 = (*_values)[i];
                for (size_t j=0; j<val1.size(); j++) {
                    for (int k=0; k<val1[j].size(); k++) {
                        double diff = val1[j][k] - (*vec._values)[i][j][k];
                        norm = std::max(norm, diff);
                    }
                }
            }
        } else if (normType==MathLib::VecNormType::NORM2) {
            for (size_t i=0; i<_values->size(); i++) {
                const IntegrationPointVectorType &val1 = (*_values)[i];
                for (size_t j=0; j<val1.size(); j++) {
                    for (int k=0; k<val1[j].size(); k++) {
                        double diff = val1[j][k] - (*vec._values)[i][j][k];
                        norm += diff*diff;
                    }
                }
            }
            norm = std::sqrt(norm);
        }
        return norm;
    }


public:
    void initialize(const MeshLib::Mesh* msh)
    {
        _msh = msh;
        _values = new VectorType(msh->getNElements());
    }

    void initialize(const MeshLib::Mesh* msh, Tvalue v0)
    {
        initialize(msh);
        for (size_t i=0; i<_values->size(); i++) {
            (*_values)[i].resize(1);
            (*_values)[i][0] = v0;
        }
    }

    void initialize(const MeshLib::Mesh* msh, VectorType* val)
    {
        initialize(msh);
        (*_values) = (*val);
    }

private:
    MeshLib::Mesh const* _msh;
    std::valarray<IntegrationPointVectorType>* _values;
};

//template<>
//class TemplateFEMIntegrationPointFunction<double> : public NumLib::ITXFunction
//{
//public:
//    typedef MathLib::TemplateVectorX<double> IntegrationPointVectorType;
//
//    inline void eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
//    {
//        size_t ele_id = x.getId();
//        IntegrationPointVectorType &gp_values = (*_values)[ele_id];
//        if (gp_values.size()==0) return;
//        double val = .0;
//        for (size_t i=0; i<gp_values.size(); i++)
//            val += gp_values[i];
//        val /= gp_values.size();
//
//        v.resize(1,1);
//        v(0,0) = val;
//        NumLib::ITXFunction::DataType mat(1,1);
//    };
//
//
//    TemplateFEMIntegrationPointFunction()
//    : _msh(nullptr), _values(nullptr)
//    {};
//
//    virtual ~TemplateFEMIntegrationPointFunction()
//    {
//        if (_values!=0) {
//            delete _values;
//            _values = 0;
//        }
//    }
//
////    ///
////    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len)
////    {
////        initialize(dis, len);
////    };
////
////    ///
////    TemplateFEMIntegrationPointFunction(MyDiscreteSystem* dis, size_t len, Tvalue &v0)
////    {
////        initialize(dis, len);
////        for (size_t i=0; i<len; i++) {
////            (*_values)[i].resize(1);
////            (*_values)[i][0] = v0;
////        }
////    };
//
//    ///
//    explicit TemplateFEMIntegrationPointFunction(const TemplateFEMIntegrationPointFunction &src)
//    {
//        initialize(src._msh);
//        (*this->_values) = (*src._values);
//        //_values = 0;
//    };
//
//    ///
//    TemplateFEMIntegrationPointFunction<double>* clone() const
//    {
//        return new TemplateFEMIntegrationPointFunction<double>(*this);
//    };
//
////    ///
////    const MeshLib::IMesh* getMesh() const
////    {
////        return this->_discrete_system->getMesh();
////    }
//
////    ///
////    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
////    {
////        switch (x.getIdObjectType()) {
////        case NumLib::TXPosition::IntegrationPoint:
////            {
////                size_t ele_id = x.getId(0);
////                size_t gp_id = x.getId(1);
////                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
////                if (gp_values.size()==0) return;
////                v = gp_values[gp_id];
////            }
////            break;
////        case NumLib::TXPosition::Element:
////            {
////                // calculate mean value
////                size_t ele_id = x.getId();
////                IntegrationPointVectorType &gp_values = (*_values)[ele_id];
////                if (gp_values.size()==0) return;
////                Tvalue val = gp_values[0];
////                val *= .0;
////                for (size_t i=0; i<gp_values.size(); i++)
////                    val += gp_values[i];
////                val /= gp_values.size();
////                v = val;
////            }
////            break;
////        default:
////            break;
////        }
////    };
//
//    void setIntegrationPointValue( size_t i_e, size_t ip, double &q )
//    {
//        assert(ip<(*_values)[i_e].size());
//        (*_values)[i_e][ip] = q;
//    }
//
//    void setNumberOfIntegationPoints(size_t i_e, size_t n)
//    {
//        (*_values)[i_e].resize(n);
//    }
//
//    const IntegrationPointVectorType& getIntegrationPointValues(size_t i_e) const
//    {
//        return (*_values)[i_e];
//    }
//
//    bool hasIntegrationPointValues(size_t i_e) const
//    {
//        return (*_values)[i_e].size()>0;
//    }
//
//    const MathLib::IVector* getElementValues() const
//    {
//        return _values;
//    }
//
//    void printout() const
//    {
//        std::cout << "integration_pt_values = ";
//        for (size_t i=_values->getRangeBegin(); i<_values->getRangeEnd(); ++i) {
//            const IntegrationPointVectorType &val1 = (*_values)[i];
//            std::cout << "(";
//            for (size_t j=0; j<val1.size(); ++j) std::cout << val1[j] << " ";
//            std::cout << ") ";
//        }
//        std::cout << std::endl;
//    }
//
//public:
//    void initialize(MeshLib::Mesh* const msh)
//    {
//        _msh = msh;
//        _values = new MathLib::DenseVector<IntegrationPointVectorType>(_msh->getNElements());
//    }
//
//    void initialize(MeshLib::Mesh* const msh, double v0)
//    {
//        initialize(msh);
//        for (size_t i=0; i<_values->size(); i++) {
//            (*_values)[i].resize(1);
//            (*_values)[i][0] = v0;
//        }
//    }
//
//private:
//    const MeshLib::Mesh* _msh;
//    MathLib::DenseVector<IntegrationPointVectorType>* _values;
//};
//
//typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar;


typedef TemplateFEMIntegrationPointFunction<MathLib::LocalVector> FEMIntegrationPointFunctionVector;

//typedef TemplateFEMIntegrationPointFunction<LocalVector> FEMIntegrationPointFunctionVector;

} //end
