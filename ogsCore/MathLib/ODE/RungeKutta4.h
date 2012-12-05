/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file RungeKutta4.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

namespace MathLib
{

/**
 * \brief 4th order Runge Kutta ODE solver
 *
 * \tparam F_DXDT   Function class evaluating F(t,X)
 */
template <class F_DXDT, class T_VALUE>
class RungeKutta4
{
public:
    /**
     * solve dX/dt = F(t,X)
     */
    template <class T_TIM_VECTOR, class T_Y_VECTOR>
    void solve(F_DXDT &f, const T_TIM_VECTOR &t, int n, T_Y_VECTOR &y) const
    {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            solve(f, t[i], dt, y[i], y[i+1]);
        }
    };

    /**
     * solve dX/dt = F(t,X)
     */
    void solve(F_DXDT &f, double t0, double dt, T_VALUE &y0, T_VALUE &y) const
    {
        T_VALUE k1 = dt*f(t0,y0);
        T_VALUE k2 = dt*f(t0+dt/2.0, y0+k1/2.0);
        T_VALUE k3 = dt*f(t0+dt/2.0, y0+k2/2.0);
        T_VALUE k4 = dt*f(t0+dt, y0+k3);
        y = y0 + (k1+2.0*k2+2.0*k3+k4)/6.0;
    };
};

} //end
