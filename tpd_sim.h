//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_TPD_SIM_H
#define TPD_SIGNAL_DECOMPOSITION_TPD_SIM_H

#endif //TPD_SIGNAL_DECOMPOSITION_TPD_SIM_H


#include "lin_alg.h"
#include "factor_operations.h"
#include "probabilistic_graph.h"
#include "storage_adaptors.h"
#include "interpolation.h"
#include "optimization.h"
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace lin_alg;
using namespace boost::numeric::odeint;
namespace ublas = boost::numeric::ublas;


/**
 * @brief: Structure that represents a surface chemical system.
 */
struct site_system{

    double mu_E;  // mean activation energy
    double sigma_squared_E;// variance in activation energy

    double mu_v; // mean pre-exponential, this does not vary with mu_E or time so variance is not needed

    double ramp;
    double T_i;

    void operator()( double C, double dC, double t );

    site_system() = default;
    site_system(double mu_E, double sigma_squared_E, double mu_v, double ramp, double T_i) :
            mu_E(mu_E), sigma_squared_E(sigma_squared_E),
            mu_v(mu_v), ramp(ramp), T_i(T_i) {}
};

/**
 * @brief overload of () operator for site system, computes the reaction rate dC as a function of concentration C and time t
 * @param C : concentration at time t
 * @param dC : desorption rate at time t
 * @param t : time
 */

void site_system::operator()(const double C, double dC, double t) {
    // compute temperature
    double T = T_i + ramp * t;
    const double R = 8.3144598; // J/molK

    // compute the integral over the Arrhenius term multiplied by the activation energy distribution
    double term1 = sigma_squared_E/(2.0 * pow(R*T,2.0));
    double term2 = mu_E/(R*T);
    double Ea_int = exp(term1 - term2);// <- this is time dependent, needs to be computed at each step


    // compute the rate of surface concentration change
    dC = C*Ea_int*mu_v;
}

/**
 * @brief: observer object for storing state and time of system during ode integration
 */

struct push_back_state_and_time{
    vector<double>& m_states;
    vector<double>& m_times;

    push_back_state_and_time(vector<double> &m_states, vector<double> &m_times) : m_states(m_states), m_times(m_times) {}

    void operator()(const double &x, double t){
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

/**
 * @brief: Data structure for storing results of simulation
 */

struct composite_signal{
    double T_i;
    double ramp;
    vector<ublas::vector<double>> states;
    vector<ublas::vector<double>> rates;
    vector<ublas::vector<double>> times;

    double flow_rate = 1;
    double T = 300;
    const double R = 8.3144598;
    vector<ublas::vector<double>> partial_pressures;
    ublas::vector<double> total_pressure;

    ublas::vector<double> compute_temperature(int i) const;
    void compute_pressures();
    void interpolate_total_pressure();
};

/**
 * @brief: computes the temperature vector corresponding to simulation i
 * @param i : index of desired site simulation
 * @return : vector of temperatures corresponding to times[i]
 */
ublas::vector<double> composite_signal::compute_temperature(int i) const{
    long N = rates[i].size();
    ublas::vector<double> T (N);
    for(int j=0; j<N; ++j){
        T[j] = T_i + ramp*times[i][j];
    }
}

/**
 * @brief: Computes partial pressures from simulated desorption rates
 */
void composite_signal::compute_pressures() {
    long N = rates.size();
    partial_pressures.resize(N);
    for(int i=0; i<N; ++i){
        for(int j=0; j<rates[i].size(); ++j){
            partial_pressures[i][j] = R*T*rates[j]/flow_rate;
        }
    }
}

/**
 * @brief: Interpolates total pressure from partial pressures computed from simulation
 */

void composite_signal::interpolate_total_pressure() {
    long N = partial_pressures.size();
    total_pressure = partial_pressures[0];

    for(int i=1; i<N; ++i){
        total_pressure += analysis::lin_interp(times[i], partial_pressures[i], times[0]);
    }
}

/**
 * @brief: Computes the TPD signal from a kinetically dominated surface system
 * @param surface_system : probabilistic model representing a surface system
 * @param T_i : (double) initial temperature [K]
 * @param T_f : (double) final temperature [K]
 * @param ramp : (double) ramp rate [K/s]
 * @return
 */

composite_signal kinetic_signal(gml::probabilistic_graph surface_system, double T_i, double T_f, double ramp){
    // compute the signal parameters for each site
    gml::variable *site = surface_system.get_variable("site");
    gml::factor *phi_Ea = surface_system.get_factor("Ea");// copy the activation energy factor
    gml::factor *phi_v = surface_system.get_factor("v"); // copy the pre-exponential factor

    // initialize vector of signal vectors, each signal vector associated with a site
    composite_signal s;
    s.states.resize(site->cardinality);
    s.rates.resize(site->cardinality);
    s.times.resize(site->cardinality);

    // compute initial and final time
    double t_i = 0;
    double t_f = (T_f - T_i)/ramp;

    // iterate over each site, this should only be tabulated over the number of sites-update for multiple gases
    for(int i=0; i<site->cardinality; ++i){

        double mu_E = phi_Ea->canonical_table[i].mu(0,0);// should be one-dimensional
        double s_E = phi_Ea->canonical_table[i].Sigma(0,0);
        double mu_v = phi_Ea->canonical_table[i].mu(0,0);

        double C_i = exp(site->canonical_column[i].g); // assumes site is represented in surface concentration space, not probability space

        // temporarily store signal and time in these vectors
        vector<double> state;
        vector<double> time;

        // run surface site simulation
        size_t num_of_steps = integrate_const(make_dense_output<runge_kutta_dopri5<ublas::vector<double>>> (1.0e-3, 1.0e-3),
                                              site_system(mu_E, s_E, mu_v, ramp, T_i), C_i,
                                              0, t_f, 1.0, push_back_state_and_time(state, time));

        // compute the desorption rate from the integrated state and time values
        site_system s_sys (mu_E, s_E, mu_v, ramp, T_i);
        vector<double> rate (state.size());
        for(int i=0; i<state.size(); ++i){
            s_sys()(state, rate, time);// double check this, might want to test site_system struct in general
        }
        s.states[i] = ublas::make_vector_from_pointer(state.size(), &state);
        s.times[i]  = ublas::make_vector_from_pointer(time.size(), &time);
        s.rates[i]  = ublas::make_vector_from_pointer(rate.size(), &rate);
        s.ramp = ramp;
        s.T_i = T_i;
    }

    return s;
}

struct experimental_signal{
    double T_i, T_f, ramp;
    ublas::vector<double> time;
    ublas::vector<double> temperature;
    ublas::vector<double> signal;
    double calibration;

    probabilistic_graph fit_model;

    /**
     * @brief: computes the residual vector between the stored experimental data and modeled data
     * @param s
     * @return
     */
    ublas::vector<double> compute_residuals( const composite_signal s ){
        long N = s.total_pressure.size();
        ublas::vector<double> T_sim = s.compute_temperature(0);
        ublas::vector<double> p_hat = analysis::lin_interp(temperature, signal, T_sim);
        ublas::vector<double> residuals (N);
        for(int i=0; i<N; ++i){
            residuals[i] = pow(p_hat[i] - s.total_pressure[i], 2.0);
        }
        return residuals;
    }

    /**
     * @brief: computes the residuals given a vector of model parameters
     * @param kinetic_parameters : reaction kinetic parameters
     * @return : redidual vector
     */

    ublas::vector<double> score_kinetic_model( ublas::vector<double> kinetic_parameters ){
        // assumes first N kinetic parameters are mean activation energies, second N are activation energy variances,
        // and third N are mean vibrational frequencies, and fourth N are total site coverages
        long N = kinetic_parameters.size()/4;

        gml::variable site ("site", N);
        gml::variable Ea ("Ea", N);
        gml::variable v ("v", N);

        for(int i=0; i<N; ++i){
            Ea.canonical_column[i].mu.resize(1,1);
            Ea.canonical_column[i].mu(0,0) = kinetic_parameters[i];
            Ea.canonical_column[i].Sigma.resize(1,1);
            Ea.canonical_column[i].Sigma(0,0) = kinetic_parameters[N+i];

            v.canonical_column[i].mu.resize(1,1);
            v.canonical_column[i].mu(0,0) = kinetic_parameters[2*N+i];

            site.canonical_column[i].g = log(kinetic_parameters[3*N+i]);
        }

        gml::factor phi_1 ("site", &site);
        gml::factor phi_2 ("Ea", &Ea);
        gml::factor phi_3 ("v", &v);

        vector<gml::factor*> factors = {
                &phi_1, &phi_2, &phi_3
        };

        gml::probabilistic_graph surface_system (factors);
        composite_signal s = kinetic_signal(surface_system, T_i, T_f, ramp);
        ublas::vector<double> residuals = compute_residuals(s);
        return residuals;
    }

    opt_lib::optimum fit_kinetic_model(const ublas::vector<double> initial_guess, opt_lib::learning_parameters params){
        // fit model with gradient descent
        opt_lib::optimum lsq = opt_lib::gradient_descent( &score_kinetic_model, initial_guess, params );


        // store parameters in probabilistic graph
        ublas::vector<double> kinetic_parameters = lsq.arg_max;
        long N = kinetic_parameters.size()/4;

        gml::variable site ("site", N);
        gml::variable Ea ("Ea", N);
        gml::variable v ("v", N);

        for(int i=0; i<N; ++i){
            Ea.canonical_column[i].mu.resize(1,1);
            Ea.canonical_column[i].mu(0,0) = kinetic_parameters[i];
            Ea.canonical_column[i].Sigma.resize(1,1);
            Ea.canonical_column[i].Sigma(0,0) = kinetic_parameters[N+i];

            v.canonical_column[i].mu.resize(1,1);
            v.canonical_column[i].mu(0,0) = kinetic_parameters[2*N+i];

            site.canonical_column[i].g = log(kinetic_parameters[3*N+i]);
        }

        gml::factor phi_1 ("site", &site);
        gml::factor phi_2 ("Ea", &Ea);
        gml::factor phi_3 ("v", &v);

        vector<gml::factor*> factors = {
                &phi_1, &phi_2, &phi_3
        };

        gml::probabilistic_graph surface_system (factors);

        // return residuals and optimal parameters
        return lsq;
    }
};



