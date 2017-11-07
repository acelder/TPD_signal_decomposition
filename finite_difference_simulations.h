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

typedef boost::numeric::ublas::vector< double > boost_vector;
typedef boost::numeric::ublas::matrix< double > matrix_type;

/*
 * Data structure to hold the system parameters for internal diffusion processes in microporous particles
 */
struct microporous_parameters{
    // parameter values
    double rho;             // site density
    double nu;              // pre-exponential
    double Ea;              // activation energy of desorption
    double steric_factor;
    double adsorbate_mass;
    double site_mass;
    double adsorbate_radius;
    double site_radius;
    double D_0;             // Diffusion pre-exponential
    double Ea_D;            // Diffusion activation energy
    double ramp_rate;
    double initial_temp;
    double final_temp;
    double characteristic_length;
};


/*
 * Data structure to hold system parameters for internal diffusion in mesoporous flims
 */

struct mesoporous_parameters{
    double rho;           // site density
    double nu;            // pre-exponential
    double Ea;            // activation energy of desorption
    double steric_factor;
    double adsorbate_mass;// kg
    double site_mass;     // kg
    double adsorbate_radius;// m
    double site_radius;   // m
    double pore_diameter;
    double initial_temperature;
    double ramp_rate;
    double initial_temp;
    double final_temp;
    double characteristic_length;
};

/*
 * Data structure to hold system parameters for external surfaces
 */

struct external_parameters{
    double nu;          // pre-exponential
    double Ea;          // activation energy of desorption
    double initial_temperature;
    double ramp_rate;   //
    double initial_temp;
    double final_temp;
};

/*
 * Data structure to hold simulation results
 */

struct sim_results{
    vector<boost_vector> concentration_profile;
    vector<double> time;
    vector<double> temperature;
};

//boost::numeric::ublas::matrix< double > compute_dynamic_concentration(parameters params);

/*
 * boost function wrappers for dx/dt
 */
struct microporous_plate_system{
    microporous_parameters p;

    void operator()(const boost_vector &x, boost_vector &dxdt, double t);

    microporous_plate_system(const microporous_parameters &p);
};

struct microporous_sphere_system{
    microporous_parameters p;

    void operator()(const boost_vector &x, boost_vector &dxdt, double t);

    microporous_sphere_system(const microporous_parameters &p);
};

struct mesoporous_film_system{
    mesoporous_parameters p;
    void operator()(const boost_vector &x, boost_vector &dxdt, double t);

    mesoporous_film_system(const mesoporous_parameters &p);
};

struct surface_system{
    external_parameters p;
    void operator()(const boost_vector &x, boost_vector &dxdt, double t);

    surface_system(const external_parameters &p);
};

/*
 * observer function
 */

struct push_back_state_and_time{
    vector< boost_vector >& m_states;
    vector<double>& m_times;

    push_back_state_and_time(vector<boost_vector> &m_states, vector<double> &m_times);

    void operator()(const boost_vector &x, double t);
};

sim_results simulate_process(const microporous_parameters& p_internal,
                             bool particle_is_sphere, const boost_vector& initial_coverage);



#endif //TPD_SIGNAL_DECOMPOSITION_FINITE_DIFFERENCE_SIMULATIONS_H
