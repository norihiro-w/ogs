/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionAnalytical.h
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Analytically defined time-space function
 */
class TXFunctionAnalytical : public ITXFunction
{
public:
    explicit TXFunctionAnalytical(const std::string &str_expression)
    : _exp(str_expression), _a0(.0), _b0(.0), _c0(.0), _d0(.0)
    {
        interpret(_exp);
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionAnalytical() {};

    virtual void eval(const TXPosition x, double &val) const
    {
        val = getValue(x.getSpace()[0], x.getSpace()[1], x.getSpace()[2]);
    }

    virtual TXFunctionAnalytical* clone() const
    {
        return new TXFunctionAnalytical(_exp);
    }

private:
    double getValue(double x, double y, double z) const
    {
        return _a0 + _b0 * x + _c0 * y + _d0 * z;
    }

    void interpret(const std::string &str_buff)
    {
        const char seps[] = "+\n";
        const char seps1[] = "*";
        char* pch;
        double f_buff;

        _a0 = _b0 = _c0 = _d0 = 0.0;

        std::vector<std::string> tokens;
        std::stringstream buff;

        pch = strtok(const_cast<char*> (str_buff.c_str()), seps);
        buff << pch;
        buff >> _a0;
        buff.clear();
        while (pch != NULL)
        {
            pch = strtok(NULL, seps);
            if (pch == NULL)
                break;
            std::string token = pch;
            tokens.push_back(token);
        }
        for (size_t k = 0; k < tokens.size(); k++)
        {
            pch = strtok(const_cast<char*> (tokens[k].c_str()), seps1);
            buff << pch;
            buff >> f_buff;
            buff.clear();
            pch = strtok(NULL, seps1);
            switch (pch[0])
            {
            case 'x':
                _b0 = f_buff;
                break;
            case 'y':
                _c0 = f_buff;
                break;
            case 'z':
                _d0 = f_buff;
                break;
            }
        }
    }

private:
    std::string _exp;
    double _a0, _b0, _c0, _d0;
};

} //end


