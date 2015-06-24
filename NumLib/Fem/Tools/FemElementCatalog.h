/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElementFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "NumLib/Fem/FiniteElement/IFemElement.h"

#include "FemElementFactoryBase.h"
#include "FemElementFactoryImpl.h"

namespace NumLib
{

/**
 * \brief Catalog of registered finite element types
 *
 * Finite element types are distinguished by unique id which is integer.
 * Unique id can be specified by outside or can be automatically assigned.
 */
class FemElementCatalog
{
public:
    FemElementCatalog() : _max_id(0) {};

    ~FemElementCatalog()
    {
        for (auto p : _map_fe_factories)
            delete p.second;
    }

    /// return if a FE type with the given id exists or not
    bool exist(const int fe_id) const
    {
        return (_map_fe_factories.count(fe_id)>0);
    };

    /**
     * register new finite element type with the given id
     * @param fe_id     Unique ID for this finite element type
     * @return the given ID if succeed, otherwise -1
     */
    template<class T_FE>
    int registerFeType(const int fe_id);

    /**
     * register new finite element type and automatically assign id
     * @return the automatically assigned ID if succeed, otherwise -1
     */
    template<class T_FE>
    int registerFeType();

    /**
     * create an object of the given finite element type
     * @param fe_id
     * @param msh
     * @return a pointer to IFiniteElement
     */
    IFiniteElement* createFeObject(const int fe_id) const;

private:
    std::map<int, FemElementFactoryBase*> _map_fe_factories;
    int _max_id;
};


/**
 * register new finite element type with the given id
 * @param fe_id     Unique ID for this finite element type
 * @return the given ID if succeed, otherwise -1
 */
template<class T_FE>
int FemElementCatalog::registerFeType(const int fe_id)
{
    if (exist(fe_id))
        return -1;

    _map_fe_factories[fe_id] = new FemElementFactoryImpl<T_FE>;
    _max_id = std::max(_max_id, fe_id);

    return fe_id;
}

/**
 * register new finite element type and automatically assign id
 * @return the automatically assigned ID if succeed, otherwise -1
 */
template<class T_FE>
int FemElementCatalog::registerFeType()
{
    const int fe_id = _max_id + 1;
    assert(_map_fe_factories.count(fe_id)==0);
    _map_fe_factories[fe_id] = new FemElementFactoryImpl<T_FE>;
    _max_id = fe_id;

    return fe_id;
}

}
