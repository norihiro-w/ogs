/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Vector.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
//#include <vector>
#include <valarray>
#include <Eigen/Eigen>

namespace MathLib
{

template<typename T>
class TemplateVectorX
{
public:
    TemplateVectorX()
    {
    };
    explicit TemplateVectorX(size_t n)
    {
        _data.resize(n);
    };
    TemplateVectorX(size_t n, T v)
    {
        _data.resize(n, v);
    }
    TemplateVectorX(const TemplateVectorX &v)
    {
        _data = v._data;
    }

    void resize(size_t n) {_data.resize(n);};
    size_t size() const {return _data.size();};

    const T* getRaw() const {return _data.size()>0 ? &_data[0] : 0;};
    T* getRawRef() {return _data.size()>0 ? &_data[0] : 0;};

    T& operator[] (size_t idx) {
        return _data[idx];
    }
    const T& operator[] (size_t idx) const {
        return _data[idx];
    }

    TemplateVectorX<T>& operator= (T a)
    {
        _data = a;

        return *this;
    }
    TemplateVectorX<T>& operator+= (T a)
    {
        _data += a;

        return *this;
    }
    TemplateVectorX<T>& operator+= (const TemplateVectorX<T> &v)
    {
        if (_data.size()==0)
            _data.resize(v._data.size());
        _data += v._data;

        return *this;
    }
    TemplateVectorX<T>& operator-= (const TemplateVectorX<T> &v)
    {
        if (_data.size()==0)
            _data.resize(v._data.size());
        _data -= v._data;

        return *this;
    }
    TemplateVectorX<T>& operator/= (T a)
    {
        _data /= a;

        return *this;
    }
    TemplateVectorX<T> operator* (T a)
    {
        TemplateVectorX<T> v(_data.size());
        for (size_t i = 0; i < _data.size(); i++)
            v._data[i] = _data[i]*a;

        return v;
    }
    TemplateVectorX<T> operator/ (T a)
    {
        TemplateVectorX<T> v(_data.size());
        for (size_t i = 0; i < _data.size(); i++)
            v._data[i] = _data[i]/a;
        return v;
    }

    TemplateVectorX<T> operator- (const TemplateVectorX<T> &ref) const
    {
        TemplateVectorX<T> v(_data.size());
        for (size_t i = 0; i < _data.size(); i++)
            v._data[i] = _data[i] - ref._data[i];

        return v;
    }

    bool operator<(const TemplateVectorX<T> &ref) const
    {
        return this->magnitude() < ref.magnitude();
    }

    T magnitude() const
    {
        T m = 0;
        for (size_t i = 0; i < _data.size(); i++)
            m += _data[i]*_data[i];
        return std::sqrt(m);
    }

    T max() const {return _data.max();};

private:
    std::valarray<T> _data;
};


template<typename T, unsigned N>
class TemplateVector
{
public:
    TemplateVector() {};
    TemplateVector(T v)
    {
        for (size_t i = 0; i < N; i++)
            _data[i] = v;
    }
    TemplateVector(T v[N])
    {
        for (size_t i=0; i<N; i++)
            _data[i] = v[i];
    }
    //TemplateVector(T &v1, T &v2);
    //TemplateVector(T &v1, T &v2, T&v3);

    TemplateVector<T,N>& operator= (T a);
    void operator-= (const TemplateVector<T,N> &v);

    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
private:
    T _data[N];
};

template<typename T, unsigned N> TemplateVector<T,N>& TemplateVector<T,N>::operator= (T a)
{
    for (size_t i = 0; i < N; i++)
        _data[i] = a;

    return *this;
}

//template<>
//TemplateVector<double,2>::TemplateVector(double &v1, double &v2)
//{
//    getRawRef()[0] = v1;
//    getRawRef()[1] = v2;
//}
//
//template<>
//TemplateVector<double,3>::TemplateVector(double &v1, double &v2, double &v3)
//{
//    getRawRef()[0] = v1;
//    getRawRef()[1] = v2;
//    getRawRef()[2] = v3;
//}

template<typename T>
class TemplateVector<T,2>
{
public:
    TemplateVector() 
    {
        for (size_t i=0; i<2; i++)
            _data[i] = .0;
    };
    TemplateVector(T v)
    {
        for (size_t i=0; i<2; i++)
            _data[i] = v;
    }
    TemplateVector(T v1, T v2)
    {
        getRawRef()[0] = v1;
        getRawRef()[1] = v2;
    }
    TemplateVector(const TemplateVector&v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] = v._data[i];
    }

    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};

    T& operator[] (size_t idx) {
        assert (idx <= 2);
        return _data[idx];
    }
    const T& operator[] (size_t idx) const {
        assert (idx <= 2);
        return _data[idx];
    }

    TemplateVector<T,2>& operator= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] = a;

        return *this;
    }
    TemplateVector<T,2>& operator+= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] += a;

        return *this;
    }
    TemplateVector<T,2>& operator+= (const TemplateVector<T,2> &v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] += v._data[i];

        return *this;
    }
    TemplateVector<T,2>& operator-= (const TemplateVector<T,2> &v)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] -= v._data[i];

        return *this;
    }
    TemplateVector<T,2>& operator/= (T a)
    {
        for (size_t i = 0; i < 2; i++)
            _data[i] /= a;

        return *this;
    }
    TemplateVector<T,2> operator* (T a)
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v._data[i] = _data[i]*a;

        return v;
    }
    TemplateVector<T,2> operator/ (T a)
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v._data[i] = _data[i]/a;

        return v;
    }

    TemplateVector<T,2> operator- (const TemplateVector<T,2> &ref) const
    {
        TemplateVector<T,2> v;
        for (size_t i = 0; i < 2; i++)
            v[i] = _data[i] - ref._data[i];

        return v;
    }

    bool operator<(const TemplateVector<T,2> &ref) const
    {
        return this->magnitude() < ref.magnitude();
    }

    T magnitude() const
    {
        T m = 0;
        for (size_t i = 0; i < 2; i++)
            m += _data[i]*_data[i];
        return std::sqrt(m);
    }

private:
    T _data[2];
};

template<typename T>
class TemplateVector<T,3>
{
public:
    TemplateVector(T v1, T v2, T v3)
    {
        getRawRef()[0] = v1;
        getRawRef()[1] = v2;
        getRawRef()[2] = v3;
    }
    const T* getRaw() const {return _data;};
    T* getRawRef() {return _data;};
private:
    T _data[3];
};

typedef TemplateVector<double, 2> Vector2D;
typedef TemplateVector<double, 3> Vector3D;
} //end mathLib

inline std::ostream& operator<<(std::ostream& output, const MathLib::TemplateVector<double,2>& p)
{
    output << "[";
    for (size_t i=0; i<2; i++)
        output << p[i] << " ";
    output << "]";
    return output;  // for multiple << operators.
}

inline std::ostream& operator<<(std::ostream& output, const MathLib::TemplateVectorX<double>& p)
{
    output << "[";
    for (size_t i=0; i<p.size(); i++)
        output << p[i] << " ";
    output << "]";
    return output;  // for multiple << operators.
}

namespace std
{

//inline MathLib::TemplateVectorX<double> abs(const MathLib::TemplateVectorX<double> &v)
//{
//    MathLib::TemplateVectorX<double> r(v.size());
//    for (size_t i=0; i<v.size(); i++)
//        r[i] = std::abs(v[i]);
//    return r;
//}

template <typename T>
inline MathLib::TemplateVectorX<T> abs(const MathLib::TemplateVectorX<T> &v)
{
    MathLib::TemplateVectorX<T> r(v.size());
    for (size_t i=0; i<v.size(); i++)
        r[i] = std::abs(v[i]);
    return r;
}


inline MathLib::TemplateVector<double,2> abs(const MathLib::TemplateVector<double,2> &v)
{
    MathLib::TemplateVector<double,2> r;
    for (size_t i=0; i<2; i++)
        r[i] = std::abs(v[i]);
    return r;
}

inline double max(double arg0, const MathLib::TemplateVector<double,2> &arg1)
{
    double r = arg0;
    for (size_t i=0; i<2; i++)
        r = std::max(r, arg1[i]);
    return r;
}

//template <typename T>
//inline double max(double arg0, const MathLib::TemplateVectorX<T> &arg1)
//{
//    double r = arg0;
//    for (size_t i=0; i<arg1.size(); i++)
//        r = std::max(r, arg1[i]);
//    return r;
//}

inline double max(double arg0, const MathLib::TemplateVectorX<Eigen::VectorXd> &arg1)
{
    double r = arg0;
    for (size_t i=0; i<arg1.size(); i++) {
        r = std::max(r, arg1[i].maxCoeff());
    }

    return r;
}

} //end std
