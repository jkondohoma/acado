
/**
 * Author: Jaelle D. Kondohoma
*/

#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>

int main() {
    USING_NAMESPACE_ACADO

    DifferentialState   x1,x2,x3 ;     // the differential states
    Control u          ;     // the control input u
    Parameter  T          ;     // the time horizon T
    DifferentialEquation  f( 0.0, T );     // the differential equation  

    // Create an optimal control problem
     OCP ocp( 0.0, T );   // time horizon of the OCP: [0,T]
     ocp.minimizeMayerTerm(T);

    // Dynamics constraints
    f << dot(x1) == x2;
    f << dot(x2) == -x2 + u;
    f << dot(x3) == (x1*x1) + (x2) + (0.005*(u));
    // f << dot(x3) == (x1*x1) + (x2*x2) + (0.005*(u*u));

    // Initial conditions
    ocp.subjectTo( f ); 
    ocp.subjectTo( AT_START, x1 == 0 );
    ocp.subjectTo( AT_START,x2 == -1 );
    ocp.subjectTo( AT_START, x3 == 0 );

    //Terminal Conditions
    ocp.subjectTo( AT_END, T == 0 );


    //Boundary Conditions
    // ocp.subjectTo(2 <= T);
   //ocp.subjectTo((x2 - 8) * ((T-0.5) * (T-0.5)) <= 0);
    ocp.subjectTo(x1 <= 0);


    

    // Solve the optimization problem
    GnuplotWindow window;
    window.addSubplot( x1, "x1" );
    window.addSubplot( x2, "x2" );
    window.addSubplot( x3, "x3" );
    window.addSubplot( u, "THE CONTROL INPUT u" );

    // Display the results
    OptimizationAlgorithm algorithm(ocp);
    algorithm << window;


    //modify settings
   // OptionsName INTEGRATOR_TYPE, INTEGRATOR_TOLERANCE, DiscretizationType, KKTtolerance;
    // algorithm.set(INITIAL_INTEGRATOR_STEPSIZE,1e10 );
    // algorithm.set(MIN_INTEGRATOR_STEPSIZE, 1e50 );
    // algorithm.set(KKT_TOLERANCE, 0.1 );
    algorithm.set(MAX_NUM_ITERATIONS, 10000 );
    algorithm.set(INTEGRATOR_TOLERANCE,1e-4 );
    algorithm.set(HESSIAN_APPROXIMATION, FULL_BFGS_UPDATE);
    algorithm.set(DISCRETIZATION_TYPE, MULTIPLE_SHOOTING);
    algorithm.set(INTEGRATOR_TYPE, INT_BDF);
    

    const clock_t begin_time = clock();
    algorithm.solve();  // solves the problem.
    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

    return 0;
}
