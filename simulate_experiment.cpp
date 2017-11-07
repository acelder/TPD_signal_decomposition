//
// Created by Alexander Elder on 11/6/17.
//

#include "finite_difference_simulations.h"

push_back_state_and_time::push_back_state_and_time(vector<boost_vector> &m_states, vector<double> &m_times) : m_states(
        m_states), m_times(m_times) {}

void push_back_state_and_time::operator()(const boost_vector &x, double t){
    m_states.push_back(x);
    m_times.push_back(t);
}

sim_results simulate_process(const microporous_parameters& p_internal,
                                bool particle_is_sphere, const boost_vector& initial_coverage){
    double t_f = (p_internal.initial_temp - p_internal.final_temp)/p_internal.ramp_rate;
    vector<boost_vector> states;
    vector<double> times;
    size_t num_of_steps;
    if(particle_is_sphere){
         num_of_steps = integrate_const(make_dense_output< runge_kutta_dopri5< boost_vector > >( 1.0e-6 , 1.0e-6 ),
                                        microporous_sphere_system(p_internal), initial_coverage,
                                        0.0, t_f, 1, push_back_state_and_time(states, times));

    }
    else{
         num_of_steps = integrate_const(make_dense_output< runge_kutta_dopri5< boost_vector > >( 1.0e-6 , 1.0e-6 ),
                                        microporous_plate_system(p_internal), initial_coverage,
                                        0.0, t_f, 1, push_back_state_and_time(states, times));
    }

    sim_results results;
    results.concentration_profile = states;
    results.time = times;
    results.temperature.resize(num_of_steps);

    for(int i=0; i<num_of_steps; ++i){
        results.temperature[i] = p_internal.initial_temp + p_internal.ramp_rate * times[i];
    }

    return results;
}
