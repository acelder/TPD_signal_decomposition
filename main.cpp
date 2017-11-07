/*
 * rosenbrock4.cpp
 *
 * Copyright 2009-2012 Karsten Ahnert
 * Copyright 2009-2012 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "finite_difference_simulations.h"


using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;


int main( int argc , char **argv ) {

    microporous_parameters p;
    p.initial_temp = 120; // K
    p.final_temp   = 300; // K
    p.ramp_rate    = 0.25;// K/s
    p.rho          = 1e15;
    p.nu           = 1e13;
    p.Ea           = 5e4;
    p.adsorbate_mass = 2.99e-26;
    p.site_mass    = 1.0552e-25;
    p.site_radius  = 0.071e-9;
    p.adsorbate_mass = 0.14e-9;
    p.steric_factor = 1.0;
    p.characteristic_length = 5e-8;
    p.D_0 = 1e-11;
    p.Ea_D = 3e4;

    boost_vector v (20);
    boost_vector dvdx(20);

    int i,j;
    for(i=0; i<20; ++i){
        v[i] = 1;
    }
    vector< boost_vector > states;
    vector<double> times;


    size_t num_steps = integrate_const(make_dense_output< runge_kutta_dopri5<boost_vector>>(1.0e-6, 1.0e-6),
                                       microporous_plate_system(p), v, 0.0, 900.0, 1.0,
                                       push_back_state_and_time(states ,times));


    cout << num_steps;
//    ofstream out_file;
//    out_file.open("test.csv");
//
//    for(i=0; i<)

    return 0;
}