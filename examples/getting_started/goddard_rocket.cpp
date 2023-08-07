/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/**
 *    \file   examples/getting_started/simple_ocp.cpp
 *    \author Boris Houska, Hans Joachim Ferreau
 *    \date   2009
 */

#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>

int main()
{

    USING_NAMESPACE_ACADO

    DifferentialState h, v, m;      // the differential states
    Control T;                      // the control input u (thrust)
    Parameter t;                    // the time horizon t
    DifferentialEquation f(0.0, t); // the differential equation

    //  -------------------------------------
    OCP ocp(0.0, t);          // time horizon of the OCP: [0,y]
    ocp.maximizeMayerTerm(t); // the time T should be optimized

    double h_0 = 1; // Initial height
    double v_0 = 0; // Initial velocity
    double m_0 = 1; // Initial mass
    double g_0 = 1; // Gravity at the surface
    double v_c = 620;
    double T_c = 3.5;
    double T_max = T_c * g_0 * m_0;
    double h_c = 500;

    double D_c = 0.5 * v_c * m_0 / g_0;
    // double D_h_v = D_c * exp(-h_c*((h-h_0)/h_0));
    // double g_h = g_0*pow(h_0/h,2);
    double c = 0.5 * sqrt(g_0 * h_0);  //Thrust-to-fuel mass

    f << dot(h) == v;                                                                           // an implementation
    f << dot(v) == ((T - (D_c * exp(-h_c * ((h - h_0) / h_0)))) / m) - (g_0 * pow(h_0 / h, 2)); // of the model equations
    f << dot(m) == -T/c;                                                               // for the rocket.

    ocp.subjectTo(f);                  // maximize t s.t. the model,
    ocp.subjectTo(AT_START, h == h_0); // the initial values for s,
    ocp.subjectTo(AT_START, v == v_0); // bv,
    ocp.subjectTo(AT_START, m == m_0); // and m,
    ocp.subjectTo(AT_START, T == T_max);

    ocp.subjectTo(AT_END, T == 0); // the terminal constraints for s

    ocp.subjectTo(0 <= T <= T_max); // the control input u,

    //  -------------------------------------

    GnuplotWindow window;
    window.addSubplot(h, "THE ALTITUDE s");
    window.addSubplot(v, "THE VELOCITY v");
    window.addSubplot(m, "THE MASS m");
    window.addSubplot(T, "THE THRUST T");

    OptimizationAlgorithm algorithm(ocp); // the optimization algorithm
    algorithm << window;
    algorithm.solve(); // solves the problem.

    return 0;
}
