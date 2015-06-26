/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <algorithm>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/LinAlgBuilder.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "MathLib/DataType.h"

#include "MeshLib/Mesh.h"

#include "NumLib/Function/TXPosition.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "NumLib/Fem/Tools/FemElementObjectContainer.h"
#include "NumLib/Fem/PolynomialOrder.h"


namespace NumLib
{

/**
 * \brief Template class for FEM node-based functions (assuming Lagrangian elements)
 *
 * This class represents the following
 * u^h(x) = N(x)*u_i
 *
 * @tparam Tvalue Nodal value type, e.g. double, vector
 */
template<typename Tvalue>
class TemplateFEMNodalFunction : public NumLib::ITXDiscreteFunction
{
public:
    typedef TemplateFEMNodalFunction<Tvalue> MyClassType;
    typedef MathLib::IVector MyVector;

    ///
    explicit TemplateFEMNodalFunction(MathLib::LinAlgLibType libType)
    : _linAlgLibType(libType), _msh(nullptr), _nodal_values(nullptr), _order(PolynomialOrder::Linear), _feObjects(nullptr)
    {};

    /// @param dis         Discrete system
    /// @param order     Polynomial order
    /// @param v0        initial value
    TemplateFEMNodalFunction(MathLib::LinAlgLibType libType, MeshLib::Mesh *msh, PolynomialOrder order, Tvalue v0)
    : _linAlgLibType(libType)
    {
        initialize(msh, order);
        resetNodalValues(v0);
    }

    /// @param dis         Discrete system
    /// @param order Polynomial order
    TemplateFEMNodalFunction(MathLib::LinAlgLibType libType, MeshLib::Mesh *msh, PolynomialOrder order)
    : _linAlgLibType(libType)
    {
        initialize(msh, order);
    }

    /// @param org source object for copying
    explicit TemplateFEMNodalFunction(const MyClassType &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
    	delete _nodal_values;
        delete _feObjects;
    }

    ///
    TemplateFEMNodalFunction &operator=(const TemplateFEMNodalFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    MyClassType* clone() const
    {
        return new MyClassType(*this);
    };

    ///
    size_t getNumberOfNodes() const {return _nodal_values->size();};


    /// evaluate this function at the given point
    virtual void eval(const NumLib::TXPosition x, NumLib::ITXFunction::DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
            {
                //TODO v = (*_nodal_values)[x.getId()];
            }
            break;
        default:
            break;
        }
    };

    /// get nodal value
    Tvalue& getValue(size_t node_id)
    {
        return (*_nodal_values)[node_id];
    }

    /// get nodal value
    const Tvalue& getValue(size_t node_id) const
    {
        return (*_nodal_values)[node_id];
    }

    ///
    void setValue(size_t node_id, Tvalue &v)
    {
        (*_nodal_values)[node_id] = v;
    }

    double diff_norm(const ITXDiscreteFunction &v, MathLib::VecNormType normType) const
    {
        MathLib::IVector* v_diff = MathLib::LinAlgBuilder::duplicateVector(*_nodal_values);
        *v_diff -= *static_cast<const MyClassType&>(v)._nodal_values;
        delete v_diff;
        switch (normType) {
        case MathLib::VecNormType::NORM1:
            return v_diff->norm1();
        case MathLib::VecNormType::NORM2:
            return v_diff->norm2();
        case MathLib::VecNormType::INFINITY_N:
            return v_diff->norm_max();
        default:
            return 0.0;
        }
    }

    /// get nodal values
    MyVector* getNodalValues() { return _nodal_values; }

    /// get nodal values
    const MyVector* getNodalValues() const { return _nodal_values; }

    /// set nodal values
    void setNodalValues( Tvalue* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            (*_nodal_values)[i+i_start] = x[i];
    }

    /// set nodal values
    void setNodalValues( const MyVector &x )
    {
        *_nodal_values = x;
    }

    /// reset nodal values with the given value
    void resetNodalValues (Tvalue &v)
    {
        *_nodal_values = v;
    }

    /// get Finite element object container
    IFeObjectContainer* getFeObjectContainer() const
    {
        return _feObjects;
    }

    void setFeObjectContainer(IFeObjectContainer *fe)
    {
        _feObjects = fe;
    }

    /// printout internal data for debugging
    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }

    /// initialize
    void initialize(const MeshLib::Mesh* msh, PolynomialOrder order)
    {
        _msh = msh;
        _order = order;
        _nodal_values = MathLib::LinAlgBuilder::generateVector(_linAlgLibType, _msh->getNNodes());
        _feObjects = nullptr;
    }

    void initialize(const MeshLib::Mesh* msh, PolynomialOrder order, Tvalue v0)
    {
        initialize(msh, order);
        resetNodalValues(v0);
    }

    PolynomialOrder getOrder() const {return _order;};

private:
    /// Assign this object from the given object
    void assign(const MyClassType &org)
    {
        initialize(*org._mshs, org._order);
        for (size_t i=org._nodal_values->getRangeBegin(); i<org._nodal_values->getRangeEnd(); ++i)
            (*_nodal_values)[i] = (*org._nodal_values)[i];
    }

private:
    MathLib::LinAlgLibType _linAlgLibType;
    MeshLib::Mesh*_msh;
    MyVector* _nodal_values;
    PolynomialOrder _order;
    IFeObjectContainer* _feObjects;
};

template <>
class TemplateFEMNodalFunction<double> : public NumLib::ITXDiscreteFunction
{
public:
    inline void eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
    {
        NumLib::ITXFunction::DataType val(1,1);
        val(0,0) = (*_nodal_values)[x.getId()];
        v = val;
    };

public:
    typedef TemplateFEMNodalFunction<double> MyClassType;
    typedef MathLib::IVector MyVector;

    ///
    explicit TemplateFEMNodalFunction(MathLib::LinAlgLibType libType)
    : _linAlgLibType(libType), _msh(nullptr), _nodal_values(nullptr), _order(PolynomialOrder::Linear), _feObjects(nullptr)
    {}

    /// @param org source object for copying
    explicit TemplateFEMNodalFunction(const MyClassType &org)
    {
        assign(org);
    }

    ///
    virtual ~TemplateFEMNodalFunction()
    {
        delete _nodal_values;
    }

    ///
    TemplateFEMNodalFunction &operator=(const TemplateFEMNodalFunction &org)
    {
        assign(org);
        return *this;
    }

    /// make a clone of this object
    /// @return MathLib::IFunction*
    MyClassType* clone() const
    {
        return new MyClassType(*this);
    };

    const MeshLib::Mesh& getMesh() const {return *_msh; }

    ///
    size_t getNumberOfNodes() const {return _nodal_values->size();};

    /// get nodal value
    double getValue(size_t node_id) const
    {
        return (*_nodal_values)[node_id];
    }

    ///
    void setValue(size_t node_id, double v)
    {
       _nodal_values->set(node_id, v);
    }

    double diff_norm(const ITXDiscreteFunction &v, MathLib::VecNormType normType) const
    {
        std::unique_ptr<MathLib::IVector> v_diff(MathLib::LinAlgBuilder::duplicateVector(*_nodal_values));
        *v_diff -= *static_cast<const MyClassType&>(v)._nodal_values;
        switch (normType) {
        case MathLib::VecNormType::NORM1:
            return v_diff->norm1();
        case MathLib::VecNormType::NORM2:
            return v_diff->norm2();
        case MathLib::VecNormType::INFINITY_N:
            return v_diff->norm_max();
        default:
            return 0.0;
        }
    }

    /// get nodal values
    MyVector* getNodalValues() { return _nodal_values; }

    /// get nodal values
    const MyVector* getNodalValues() const { return _nodal_values; }

    /// set nodal values
    void setNodalValues( double* x, size_t i_start, size_t n )
    {
        for (size_t i=0; i<n; ++i)
            _nodal_values->set(i+i_start, x[i]);
    }

    /// set nodal values
    void setNodalValues( const MyVector &x )
    {
        *_nodal_values = x;
    }

	/// reset nodal values with the given value
    void resetNodalValues (double v)
    {
        *_nodal_values = v;
    }

    /// get Finite element object container
    const IFeObjectContainer* getFeObjectContainer() const
    {
        return _feObjects;
    }

    void setFeObjectContainer(const IFeObjectContainer *fe)
    {
        _feObjects = fe;
    }

    /// printout internal data for debugging
    void printout() const
    {
        std::cout << "nodal_values = ";
        for (size_t i=_nodal_values->getRangeBegin(); i<_nodal_values->getRangeEnd(); ++i)
            std::cout << (*_nodal_values)[i] << " ";
        std::cout << std::endl;
    }

    /// initialize
    void initialize(MeshLib::Mesh* msh, PolynomialOrder order)
    {
    	_msh = msh;
        _order = order;
        size_t n = _msh->getNNodes();
        _nodal_values = MathLib::LinAlgBuilder::generateVector(_linAlgLibType,n);
        _feObjects = nullptr;
    }

    template <class T_SYS>
    void initialize(T_SYS &dis, PolynomialOrder order, double v0)
    {
        initialize(dis, order, _linAlgLibType);
        resetNodalValues(v0);
    }

    PolynomialOrder getOrder() const {return _order;};

private:
    /// Assign this object from the given object
    void assign(const MyClassType &org)
    {
        _nodal_values = org._nodal_values->duplicate();
        _feObjects = org._feObjects;
    }

private:
    MathLib::LinAlgLibType _linAlgLibType;
    const MeshLib::Mesh*_msh;
    MyVector* _nodal_values;
    PolynomialOrder _order;
    const IFeObjectContainer* _feObjects;
};

///// evaluate this function at the given point
//template <class T_DIS_SYS>
//inline void TemplateFEMNodalFunction<T_DIS_SYS,double>::eval(const NumLib::TXPosition x,  NumLib::ITXFunction::DataType &v) const
//{
//    NumLib::ITXFunction::DataType val(1,1);
//    val(0,0) = (*_nodal_values)[x.getId()];
//    v = val;
//};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
typedef TemplateFEMNodalFunction<MathLib::LocalVector> FemNodalFunctionVector;

//typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
//typedef TemplateFEMNodalFunction<MathLib::LocalVector> FemNodalFunctionVector;

} //end
