//
// Functions and classes for simulating temperature programmed desorption signals.
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
#include "quadrature.h"
#include "gaussian_operations.h"
#include <iostream>
#include <numeric>
#include <fstream>
#include <vector>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/odeint.hpp>
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

    double fractional_coverage;
    double Ea;
    double nu;

    double ramp;
    double T_i;
    double P_reaction;

    void operator()( const ublas::vector<double> &C, ublas::vector<double> &dC, double t );

    site_system() = default;
    site_system( double coverage, double Ea, double nu, double ramp, double T_i) :
            fractional_coverage(coverage), Ea(Ea), nu(nu), ramp(ramp), T_i(T_i) {}
};

/**
 * @brief overload of () operator for site system, computes the reaction rate dC as a function of concentration C and time t
 * @param C : concentration at time t
 * @param dC : desorption rate at time t
 * @param t : time
 */

void site_system::operator()(const ublas::vector<double> &C, ublas::vector<double> &dC, double t) {
    // compute temperature
    double T = T_i + ramp * t;
    const double R = 8.3144598; // J/molK

    // compute the integral over the joint distribution of vibrational frequency, activation energy, and reaction probability
    double P_activation = exp(-Ea/(R*T)); // activation energy stored in mu(0)
    double P_configuration = nu; // assume zero variance in the vibrational frequency
    P_reaction = P_activation * P_configuration;

    // compute the rate of surface concentration change
    dC[0] = -C[0]*P_reaction;
}

struct site_system_jacobi{
    double Ea;
    double nu;

    double ramp;
    double T_i;
    double P_reaction;

    void operator()( const ublas::vector<double> &C, ublas::matrix<double> &J, double t);//, ublas::vector<double> dfdt );

    site_system_jacobi() = default;
    site_system_jacobi(const site_system &deriv) : Ea(deriv.Ea), nu(deriv.nu), ramp(deriv.ramp), T_i(deriv.T_i) {}
};

void site_system_jacobi::operator()(const ublas::vector<double> &C, ublas::matrix<double> &J, const double t){
    const double R = 8.3144598;
    double T = T_i + ramp*t;
    double P_activation = exp(-Ea/(R*T)); // activation energy stored in mu(0)
    double P_configuration = nu; // assume zero variance in the vibrational frequency
    J(0,0) = P_activation * P_configuration;
}

/**
 * @brief: observer object for storing state and time of system during ode integration
 */

struct push_back_state_and_time{
    vector<ublas::vector<double>>& m_states;
    vector<double>& m_times;

    push_back_state_and_time(vector<ublas::vector<double>> &m_states, vector<double> &m_times)
            : m_states(m_states), m_times(m_times) {}

    void operator()(const ublas::vector<double> &x, double t){
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
    double R = 8.3144598;

    composite_signal() = default;
    vector<ublas::vector<double>> partial_pressures;
    ublas::vector<double> total_pressure;

    ublas::vector<double> compute_temperature(int i) const;
    void interpolate_rates();
    void compute_pressures();
    void interpolate_total_pressure();
    void normalize();
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
    return T;
}

void composite_signal::interpolate_rates(){
    long N = rates.size();
    for(int j=0; j<N; ++j){
        rates[j] = analysis::lin_interp(times[j], rates[j], times[0]);
    }
}

/**
 * @brief: Computes partial pressures from simulated desorption rates
 */
void composite_signal::compute_pressures() {
    long N = rates.size();
    partial_pressures.resize(N);
    for(int i=0; i<N; ++i){
        partial_pressures[i].resize(rates[i].size());
        for(int j=0; j<rates[i].size(); ++j){
            partial_pressures[i][j] = R*T*rates[i][j]/flow_rate;
        }
        // interpolate
        partial_pressures[i] = analysis::lin_interp(times[i], rates[i], times[0]);
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
 * @brief: Normalizes the pressure for fitting uncalibrated signals
 */

void composite_signal::normalize(){
    auto it = max_element(total_pressure.begin(), total_pressure.end());
    double z = (*it);
    total_pressure = total_pressure*(1/z);// compiler might not like this

    for(int i=0; i<partial_pressures.size(); ++i){
        partial_pressures[i] = partial_pressures[i]*(1/z);
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

composite_signal kinetic_signal(vector<double> site_distribution, vector<double> activation_energies,
                                vector<double> prefactors, double T_i, double T_f, double ramp){
    // compute initial and final time
    double t = 0.0;
    double t_f = (T_f - T_i)/ramp;

    // initialize vector of signal vectors, each signal vector associated with a site
    size_t num_sites = site_distribution.size();
    composite_signal sig;
    sig.states.resize(num_sites);//
    sig.rates.resize(num_sites);
    sig.times.resize(num_sites);
    sig.ramp = ramp;
    sig.T_i = T_i;

    // iterate over each site, this should only be tabulated over the number of sites-update for multiple gases
    for(int i=0; i<num_sites; ++i){

        double Ea = activation_energies[i];
        double v  = prefactors[i];
        double site_coverage = site_distribution[i];

        ublas::vector<double> C (1, site_coverage);

        // temporarily store signal and time in these vectors


        list<pair<double, ublas::vector<double>>> state_log;

        // simulation loop
        typedef boost::numeric::odeint::implicit_euler< double > booststepper;
        auto stepper = booststepper();
        site_system system = site_system(site_coverage, Ea, v, ramp, T_i);
        site_system_jacobi jacobi = site_system_jacobi(system);
        double dt = 1;
        auto ode = make_pair(system, jacobi);
        const double min_threshold = 0.001*C[0];// might want to make this more flexible

        while(t < t_f){
            pair<double, ublas::vector<double>> state;
            state.first = t;
            state.second = C;
            state_log.push_back(state);

            stepper.do_step(ode, C, t, dt);
            t+= dt;

            // if min threshold is reached, set value at t_f and break
            if(C[0] < min_threshold){
                state.first = t_f;
                C(0) = 0;
                state.second = C;
                state_log.emplace_back(make_pair(t_f,C));

                break;
            }
        }

        // compute the desorption rate from the integrated state and time values
        vector<double> rate (state_log.size());



        // load everything to signal data structure
        sig.states[i].resize(state_log.size());
        sig.times[i].resize(state_log.size());
        sig.rates[i].resize(state_log.size());

        int j = 0;
        for(auto L: state_log){
            sig.times[i][j] = L.first;
            sig.states[i][j] = L.second(0);

            // compute and store the rate
            ublas::vector<double> temp_rate (1);
            system(L.second, temp_rate, L.first);
            sig.rates[i][j] = -temp_rate[0];
            ++j;
        }
    }

    return sig;
}

class experimental_signal{
public:
    double T_i, T_f, ramp;
    ublas::vector<double> time;
    ublas::vector<double> temperature;
    ublas::vector<double> signal;
    double calibration;

    gml::probabilistic_graph model;

    /**
     * @brief: Load cleaned experimental data from csv. Example format can be found in project files.
     * @param file
     */

    void load_csv(string file, long n_points){

        fstream infile;
        infile.open(file);
        string str;
        vector<string> spl;

        time.resize(n_points);
        temperature.resize(n_points);
        signal.resize(n_points);

        int line = 0;
        while(getline(infile, str)){
            boost::split(spl, str, boost::is_any_of(","));
            switch(line){
                case 0: T_i = stod(spl[1]);
                    break;
                case 1: T_f = stod(spl[1]);
                    break;
                case 2: ramp = stod(spl[1]);// beware of spaces
                    break;
                case 3: // do nothing on the third line
                    break;
                default:
                    int i = line - 4;
                    time[i] = stod(spl[0]);
                    temperature[i] = stod(spl[1]);
                    signal[i] = stod(spl[2]);
            }
            ++line;
        }
        infile.close();
    }

    /**
     * @brief: normalizes the signal
     */
    void normalize_signal(){
        auto it = max_element(signal.begin(), signal.end());
        double z = (*it);
        signal = signal*(1/z);
    }

    /**
     * @brief: zeros the signal
     */
    void zero_signal(){
        auto it = min_element(signal.begin(), signal.end());
        for(int i=0; i<signal.size(); ++i){

            signal[i] = signal[i]-(*it);
        }
    }

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


};





