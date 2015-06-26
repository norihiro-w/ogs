/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractMonolithicProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cassert>
#include <vector>
#include <queue>
#include "BaseLib/CodingTools.h"
#include "NumLib/IOSystem/NamedIOSystem.h"

namespace NumLib
{

/**
 * \brief MonolithicSolution
 */
template<class T_BASE>
class AbstractMonolithicSystem : public NamedIOSystem<T_BASE>
{
public:
    ///
    AbstractMonolithicSystem() {};

    ///
    virtual ~AbstractMonolithicSystem()
    {
    }

    /// check consistency
    virtual bool check() const {return true;};

private:

    DISALLOW_COPY_AND_ASSIGN(AbstractMonolithicSystem);
};

}
