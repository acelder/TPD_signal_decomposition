

#include "unit_tests.h"
#include "sensitivity_analysis.h"
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include "quadrature.h"


using namespace std;


int main( int argc , char **argv ) {

    clock_t t;
    t = clock();



//    skew_sensitivity_analysis("skew_sensitivity_1.csv");

    // set up variables
    gml::variable site ("site", 1);
    gml::variable v ("v", 1);
    gml::variable Ea ("Ea", 1);

    gml::canonical_form site1;
    site1.g = log(1);

    gml::canonical_form v1;
    v1.mu.resize(1);
    v1.mu(0) = 1e13;
    v1.Sigma.resize(1,1);
    v1.Sigma(0,0) = 1;
    v1.scope.insert("Hz");

    gml::skew_canonical_form Ea1;
    Ea1.mu.resize(1);
    Ea1.mu(0) = 5e4;
    Ea1.Sigma.resize(1,1);
    Ea1.Sigma(0,0) = 1;
    Ea1.scope.insert("energy");
    Ea1.skew.resize(1);
    Ea1.skew(0) = 0.0;

    site_system system = site_system(&Ea1, &v1, 0.25, 120);
    site_system_jacobi jacobi = site_system_jacobi(&system);

    ublas::vector<double> C (1, exp(site1.g));
    ublas::vector<double> dC (1,0);
//    system(C, dC, 300.0);

    ublas::vector<double> r (2);
    r(0)= 5e4; r(1) =  1e13;
    ublas::vector<double> x (1, 5e4);


    double rate_contribution = integrate_rate_contributions(&Ea1, &v1, 120); // why does this equal 4.2 at 120 K?
    cout << rate_contribution << endl;



    t = clock() - t;
    cout << "process finished in " << (float)t/CLOCKS_PER_SEC << " seconds." << endl;


    return 0;
}