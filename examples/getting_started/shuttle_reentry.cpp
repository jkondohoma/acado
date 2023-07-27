/**
 *    \file   examples/getting_started/shuttle_reentry.cpp
 *    \author Jaelle Kondohoma
 *    \date   2023
 */

#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>
#include <cmath>
#include <math.h>

int main()
{

    USING_NAMESPACE_ACADO

    DifferentialState h, phi, theta, v, y, psi; // the differential states
    Control alpha, beta, q;                     // the control inputs
    Parameter T;                                // the time horizon T
    DifferentialEquation f(0.0, T);             // the differential equation

    // Initial conditions
    double h_s = 260000;
    double phi_s = 0;
    double theta_s = 0;
    double alpha_s = 0;
    double beta_s = 0;
    double v_s = 25600;
    double y_s = -1;
    double psi_s = 90;

    // Final conditions, Terminal Area Energy Management (TAEM)
    double h_t = 80000;
    double v_t = 2500;
    double y_t = -5;

    // constants & expressions
    double a_hat = (180 * alpha_s) / M_PI;
    double R_e = 20902900;
    double S = 2690;
    double a_0 = -0.20704;
    double a_1 = 0.029244;
    double miu = 0.14076539 * pow(10, 17);
    double b_0 = 0.07854;
    double b_1 = -0.61592 * pow(10, -2);
    double b_2 = -0.621408 * pow(10, 3);
    double ro_0 = 0.002378;
    double h_r = 23800;
    double c_0 = 1.0672181;
    double c_1 = -0.19213774 * pow(10, -1);
    double c_2 = 0.21286289 * pow(10, -3);
    double c_3 = -0.10117249 * pow(10, -5);

    double r = R_e + h_s;
    double ro = ro_0 * exp((-h_s) / h_r);
    double expression = pow(0.0001 * v_s, 3.07);
    double q_r = 17700 * sqrt(ro * expression);
    double q_a = c_0 + (c_1 * a_hat) + (c_2 * pow(a_hat, 2)) + (c_3 * pow(a_hat, 3));
    double c_D = b_0 + (b_1 * a_hat) + (b_2 * sqrt(a_hat));
    double c_L = a_0 + (a_1 * a_hat);
    double D = (1 / 2) * c_D * S * ro * pow(v_s, 2);
    double L = (1 / 2) * c_L * S * ro * pow(v_s, 2);

    double w = 203000.0; // weight (lb)
    double g = 32.174;   // acceleration (ft/sec^2)
    double m = w / g;    // mass (slug)

    // double q = q_a * q_r;
    // double q_U = 70;

    //  -------------------------------------
    OCP ocp(0.0, T);             // time horizon of the OCP: [0,T]
    ocp.maximizeLagrangeTerm(T); // the time T should be optimized

    //  implementation of the model equations for the rocket
    f << dot(h) == v * sin(y);
    f << dot(phi) == ((v / r) * cos(y) * sin(psi)) / cos(theta);
    f << dot(theta) == (v / r) * cos(y) * cos(psi);
    f << dot(v) == (-D / m) - (g * sin(y));
    f << dot(y) == ((L / (m * v)) * cos(beta)) + (cos(y * ((v / r) - (g / v))));
    f << dot(psi) == ((1 / m * v * cos(y)) * L * sin(beta)) + ((v / (r * cos(theta))) * cos(y) * sin(psi) * cos(theta));

    ocp.subjectTo(f); // maximize T s.t. the model,
    ocp.subjectTo(AT_START, h == h_s);
    ocp.subjectTo(AT_START, phi == phi_s);
    ocp.subjectTo(AT_START, theta == theta_s);
    ocp.subjectTo(AT_START, v == v_s);
    ocp.subjectTo(AT_START, y == y_s);
    ocp.subjectTo(AT_START, psi == psi_s);

    ocp.subjectTo(AT_END, h == h_t);
    ocp.subjectTo(AT_END, v == v_t);
    ocp.subjectTo(AT_END, y == y_t);

    // boundary restrictionss
    ocp.subjectTo(0 <= h);
    ocp.subjectTo(1 <= v);
    ocp.subjectTo(-90 <= alpha <= 90);
    ocp.subjectTo(-80 <= theta <= 89);
    ocp.subjectTo(-80 <= y <= 89);
    ocp.subjectTo(-80 <= beta <= 89);
    // ocp.subjectTo(q <= q_U);
    //  -------------------------------------

    GnuplotWindow window;
    window.addSubplot(h, "THE ALTITUDE h");
    window.addSubplot(v, "THE VELOCITY v");
    window.addSubplot(phi, "THE LONGITUDE phi");
    window.addSubplot(y, "THE FLIGHT PATHT y");

    OptimizationAlgorithm algorithm(ocp); // the optimization algorithm
    algorithm << window;
    algorithm.solve(); // solves the problem.

    return 0;
}
