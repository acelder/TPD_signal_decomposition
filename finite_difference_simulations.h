//
// Created by Alexander Elder on 11/1/17.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_FINITE_DIFFERENCE_SIMULATIONS_H
#define TPD_SIGNAL_DECOMPOSITION_FINITE_DIFFERENCE_SIMULATIONS_H

#include <iostream>
#include <fstream>
#include <utility>// is this needed?
#include <boost/numeric/odeint.hpp>
//#include <boost/math/constants>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace boost::numeric::odeint;
using namespace std;
namespace phoenix = boost::phoenix;

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

/*
 * Define a data structure to hold the system parameters for internal diffusion processes
 */
struct microporous_parameters{
    // parameter values
    double rho;             // site density
    double nu;              // pre-exponential
    double Ea;              // activation energy of desorption
    double steric_factor;
    double adsorbate_mass;
    double site_mass;
    double adsorbate_radius
    double site_radius;
    double D_0;             // Diffusion pre-exponential
    double Ea_D;            // Diffusion activation energy
    double initial_temperature;
    double ramp_rate;
    double characteristic_length;
};

//boost::numeric::ublas::matrix< double > compute_dynamic_concentration(parameters params);

/*
 * boost function wrappers for dx/dt
 */
struct microporous_plate_system{
    microporous_parameters p;

    void operator()(const vector_type &x, vector_type &dxdt, double t);

    microporous_plate_system(const microporous_parameters &p);
};

struct microporous_sphere_system{
    microporous_parameters p;
    void operator()(const vector_type &x, vector_type &dxdt, double t);
};

struct mesoporous_film_system{
    mesoporous_parameters p;
    void operator()(const vector_type &x, vector_type &dxdt, double t);
};

struct surface_system{
    external_parameters p;
    void operator()(const vector_type &x, vector_type &dxdt, double t);
};


#endif //TPD_SIGNAL_DECOMPOSITION_FINITE_DIFFERENCE_SIMULATIONS_H
