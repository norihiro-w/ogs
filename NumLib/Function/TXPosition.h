/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXPosition.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

namespace NumLib
{

/**
 * \brief Position in space-time domain
 */
class TXPosition
{
public:
    enum IdObjectType
    {
        INVALID,
        Node,
        Element,
        IntegrationPoint
    };

    explicit TXPosition(const double* x_) : _t(.0), _x(x_), _obj_type(INVALID) {};
    TXPosition(double t_, const double* x_) : _t(t_), _x(x_), _obj_type(INVALID) {};
    TXPosition(IdObjectType objType, size_t n) : _t(.0), _x(0), _array_id(1, n), _obj_type(objType) {};
    TXPosition(IdObjectType objType, size_t n, const double* p) : _t(.0), _x(p), _array_id(1, n), _obj_type(objType) {};
    TXPosition(IdObjectType objType, size_t id1, size_t id2, const double* p) : _t(.0), _x(p), _obj_type(objType) 
    {
        _array_id.push_back(id1);
        _array_id.push_back(id2);
    };
    

    double getTime() const {return _t;};
    const double* getSpace() const {return _x;};
    size_t getId(size_t i=0) const {return _array_id.size()<i+1 ? BaseLib::index_npos : _array_id[i];};
    IdObjectType getIdObjectType() const {return _obj_type;};

private:
    double _t;
    const double* _x;
    std::vector<size_t> _array_id;
    IdObjectType _obj_type;

};


} //end
