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
#include "quadrature.h"
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
 * @brief compute the joint density over reaction probability, vibrational frequency, and activation energy
 * @param fv : distribution function for vibrational frequency
 * @param fEa : distribution function for activaiton energy
 * @param v : vibrational frequency
 * @param Ea : activation energy
 * @param T : temperature
 * @return
 */

double joint_density( canonical_form* fEa, canonical_form* fv, ublas::vector<double> r, double T){
    double R = 8.3144598;
    ublas::vector<double> Ea_vec (1, r(0));
    ublas::vector<double> v_vec (1, r(1));
    double density = fv->density(v_vec)*fEa->density(Ea_vec)*exp(-Ea_vec(0)/(R*T));
    return density;
}

double rate_integrand(canonical_form* fEa, canonical_form* fv, ublas::vector<double> r, double T){
    double R = 8.3144598;
    ublas::vector<double> Ea_vec (1, r(0));
    ublas::vector<double> v_vec (1, r(1));
    double rate_contribution = r(1)*fv->density(v_vec)*fEa->density(Ea_vec)*exp(-Ea_vec(0)/(R*T));
    return rate_contribution;
}

/**
 * @brief A silly implementation of skew bivariate distribution density. Keep for now given deadline constraints.
 * @param fEa
 * @param fv
 * @param r a vector in Energy-Frequency space
 * @return
 */

double skew_bivariate_density(canonical_form* fEa, canonical_form* fv,
                                             ublas::vector<double> r){
    ublas::vector<double> Ea_vec (1, r(0));
    ublas::vector<double> v_vec (1, r(1));
    double density = fv->density(v_vec)*fEa->density(Ea_vec);
    return density;
}

/**
 * Integrates over bivariate joint density by rotating and centering the distribution overlaying bivariate skew normal distribution
 * @param fv
 * @param fEa
 * @param T
 * @param error : maximum allowable error
 * @return
 *
 * still needs to be tested
 */

double integrate_rate_contributions( canonical_form* fEa, canonical_form* fv, double T, double error = 0.01){
    // This subroutine will make extensive use of function bindings so let's go ahead and use the placeholders namespace.
    using namespace std::placeholders;

    // Begin by finding the peak of the unmodified bivariate distribution. This is necessary to establish the rectangular
    // boundary of integration about the distribution.

    // Bind the skew_bivariate_density function to the input Gaussians.
    // density_func(Ea_vec, v_vec) returns the density of the input Gaussians at Ea_vec, v_vec
    function<double(ublas::vector<double>)> density_func = bind(skew_bivariate_density, fEa, fv, _1);

    // Find the peak using gradient ascent.
    ublas::vector<double> initial_guess (2); // set initial guess
    initial_guess(0) = fEa->mu(0);
    initial_guess(1) = fv->mu(0);
    opt_lib::learning_parameters params(0.1, 1000, 0.1);// set learning parameters
    opt_lib::optimum skew_peak = opt_lib::gradient_ascent(density_func, initial_guess, params);




    // In order to determine the correct bounds of integration, first locate the direction of slowest descent about the
    // center. This will be accomplished with the following rotary search algorithm.
    int n_angles = 360; // Let's just brute force this for now and search 360 degrees about the center
    double max_density = skew_peak.opt;
    vector<double> theta (n_angles,0);
    vector<double> secants (n_angles, 0);
    for(int i=0; i<n_angles; ++i){
        theta[i] = 2*pi*i/n_angles;
        // measure the slope of the secant line between the optimum and arg_opt + 0.1 at angle theta
        ublas::vector<double> r = skew_peak.arg_opt;
        r(0) += cos(theta[i])*0.1;
        r(0) += sin(theta[i])*0.1;
        double shifted_density = skew_bivariate_density(fEa, fv, r);
        secants[i] = max_density - shifted_density; // could compute as euclidean distance but not necessary here
    }
    long ind = distance(secants.begin(), min_element(secants.begin(), secants.end()));// get minimum descent direction
    double search_angle = theta[ind];


    // Now perform a one-dimensional line search with unit step size. Might want to optimize this later.
    ublas::vector<double> corner_one (2,0);
    ublas::vector<double> line_search_position = skew_peak.arg_opt;
    for(int i=0; i<1000; ++i){
        ublas::vector<double> new_line_search_position (2);
        new_line_search_position(0) = line_search_position(0) + cos(search_angle); // stick with a unit step length for now,
        new_line_search_position(1) = line_search_position(1) + sin(search_angle); // maybe implement adaptive step lengh later

        double new_density = skew_bivariate_density(fEa, fv, new_line_search_position);

        if(new_density < error){
            corner_one = new_line_search_position;
        }
        else{
            line_search_position = new_line_search_position;
        }
    }

    // With the first corner and the search angle in hand, we can compute the other three corners of the integration
    // square. Start by reflecting the first corner across the vertical axis.
    ublas::vector<double> corner_two = corner_one;
    corner_two(0) *= -1;

    // Find the next corner by reflecting the first through the origin.
    ublas::vector<double> corner_three = corner_one;
    corner_three(0) *= -1;
    corner_three(1) *= -1;

    // Find the final corner by reflecting the first across the horizontal axis
    ublas::vector<double> corner_four = corner_one;
    corner_four(1) *= -1;

    double Ea_coords [4] = {corner_one(0), corner_two(0), corner_three(0), corner_four(0)};
    double v_coords [4] = {corner_one(1), corner_two(1), corner_three(1), corner_four(1)};

    double* Ea_lower_bound = min_element(Ea_coords, Ea_coords+4);
    double* Ea_upper_bound = max_element(Ea_coords, Ea_coords+4);
    double* v_lower_bound = min_element(v_coords, v_coords+4);
    double* v_upper_bound = max_element(v_coords, v_coords+4);

    // Finally, bind the full joint distribution with the specific parameters and integrate over the selected region.
    // joint_density_func(r) returns joint probability density over frequency, activation energy, and reaction probability

    ublas::matrix<double> integration_range (2,2);
    integration_range(0,0) = *Ea_lower_bound;
    integration_range(0,1) = *Ea_upper_bound;
    integration_range(1,0) = *v_lower_bound;
    integration_range(1,1) = *v_upper_bound;

    function<double(ublas::vector<double>)> contribution_func = bind(rate_integrand, fEa, fv, _1, T);//

    ublas::vector<long> n_domains (2, 100);

    double integral = lin_alg::analysis::dbl_newton_cotes(contribution_func, n_domains, integration_range, "simpson" );
    return integral;
}

/**
 * @brief: Structure that represents a surface chemical system.
 */
struct site_system{

    canonical_form* fEa;
    canonical_form* fv;

    double ramp;
    double T_i;
    double P_reaction;

    void operator()( const ublas::vector<double> &C, ublas::vector<double> &dC, double t );

    site_system() = default;
    site_system(canonical_form* fEa, canonical_form* fv, double ramp, double T_i) :
            fEa(fEa), fv(fv), ramp(ramp), T_i(T_i) {}
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
    double P_reaction = integrate_rate_contributions(fEa, fv, T);

    // compute the rate of surface concentration change
    dC[0] = -C[0]*P_reaction;
}

struct site_system_jacobi{
    site_system* deriv;


    void operator()( const ublas::vector<double> &C, ublas::matrix<double> &J, double t);//, ublas::vector<double> dfdt );

    site_system_jacobi() = default;
    site_system_jacobi(site_system* deriv) : deriv(deriv) {}
};

void site_system_jacobi::operator()(const ublas::vector<double> &C, ublas::matrix<double> &J, const double t){
    J(0,0) = deriv->P_reaction;

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

composite_signal kinetic_signal(gml::probabilistic_graph *surface_system, double T_i, double T_f, double ramp){
    // compute the signal parameters for each site
    gml::variable *site = surface_system->variable_map["site"];
    gml::factor *phi_Ea = surface_system->factor_map["Ea"];// extract the activation energy factor
    gml::factor *phi_v = surface_system->factor_map["v"]; // extract the pre-exponential factor
    //


    // compute initial and final time
    double t = 0.0;
    double t_f = (T_f - T_i)/ramp;

    // initialize vector of signal vectors, each signal vector associated with a site
    composite_signal sig;
    sig.states.resize((size_t)site->cardinality);//
    sig.rates.resize((size_t)site->cardinality);
    sig.times.resize((size_t)site->cardinality);
    sig.ramp = ramp;
    sig.T_i = T_i;

    // iterate over each site, this should only be tabulated over the number of sites-update for multiple gases
    for(int i=0; i<site->cardinality; ++i){

        double mu_E = phi_Ea->canonical_table[i]->mu(0);
        double s_E = phi_Ea->canonical_table[i]->Sigma(0,0);
        double mu_v = phi_Ea->canonical_table[i]->mu(0);

        ublas::vector<double> C (1, exp(site->canonical_column[i]->g));

        // temporarily store signal and time in these vectors


        list<pair<double, ublas::vector<double>>> state_log;

        // simulation loop
        typedef boost::numeric::odeint::implicit_euler< double > booststepper;
        auto stepper = booststepper();
        site_system system = site_system(phi_Ea->canonical_table[i], phi_v->canonical_table[i], ramp, T_i);
        site_system_jacobi jacobi = site_system_jacobi(&system);
        double dt = 0.1;
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
        site_system s_sys (phi_Ea->canonical_table[i], phi_v->canonical_table[i], ramp, T_i);
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
            s_sys(L.second, temp_rate, L.first);
            sig.rates[i][j] = temp_rate[0];
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

    /**
     * @brief: computes the residuals given a vector of model parameters
     * @param kinetic_parameters : reaction kinetic parameters
     * @return : redidual vector
     */

    double score_kinetic_model( ublas::vector<double> kinetic_parameters ){
        // assumes first N kinetic parameters are mean activation energies, second N are activation energy variances,
        // and third N are mean vibrational frequencies, and fourth N are total site coverages
        long N = kinetic_parameters.size()/4;

        gml::variable site ("site", N);
        gml::variable Ea ("Ea", N);
        gml::variable v ("v", N);

        for(int i=0; i<N; ++i){
            Ea.canonical_column[i]->mu.resize(1,1);
            Ea.canonical_column[i]->mu(0) = kinetic_parameters[i];
            Ea.canonical_column[i]->Sigma.resize(1,1);
            Ea.canonical_column[i]->Sigma(0,0) = kinetic_parameters[N+i];

            v.canonical_column[i]->mu.resize(1,1);
            v.canonical_column[i]->mu(0) = kinetic_parameters[2*N+i];

            site.canonical_column[i]->g = log(kinetic_parameters[3*N+i]);
        }

        gml::factor phi_1 ("site", &site);
        gml::factor phi_2 ("Ea", &Ea);
        gml::factor phi_3 ("v", &v);

        vector<gml::factor*> factors = {
                &phi_1, &phi_2, &phi_3
        };

        gml::probabilistic_graph surface_system (factors);
        composite_signal s = kinetic_signal(&surface_system, T_i, T_f, ramp);
        ublas::vector<double> residuals = compute_residuals(s);
        double ssr = 0.0;
        for(auto resid: residuals){
            ssr += pow(resid, 2);
        }
        return ssr;// return sum of squared residuals
    }

    /**
     * @brief: Fits kinetic model to the loaded data given an initial guess
     * @param initial_guess
     * @param params
     */

    void fit_kinetic_model(const ublas::vector<double> initial_guess, opt_lib::learning_parameters params){


        // fit model with gradient descent

        long N = initial_guess.size();
        ublas::vector<double> x = initial_guess;
        ublas::vector<double> x_prev;
        ublas::vector<double> e (N, 0);

        for(int i=0; i<params.num_iter; ++i){

            for(int i=0; i<params.num_iter; ++i){
                ublas::vector<double> central_difference (N, 0);
                // compute the central difference approximation of the gradient
                for(int j=0; j<N; ++j){
                    ublas::vector<double> e (N, 0);
                    e(j) = opt_lib::eps;
                    ublas::vector<double> forward_difference = x + e;
                    ublas::vector<double> backward_difference = x - e;
                    double forward_eval = score_kinetic_model(forward_difference);
                    double backward_eval= score_kinetic_model(backward_difference);
                    central_difference(j) = (forward_eval - backward_eval)/(2*opt_lib::eps);

                }

                x_prev = x;
                x = x - params.rate*central_difference;// subtract the gradient to descend the valley
                if(ublas::norm_1(x - x_prev) <= params.convergence_criterion){
                    break;
                }
            }

        }

        opt_lib::optimum lsq;
        lsq.arg_opt = x;
        lsq.opt = score_kinetic_model(x);


        // store parameters in probabilistic graph
        ublas::vector<double> kinetic_parameters = lsq.arg_opt;
        N = kinetic_parameters.size()/4;// NOTE: N redefined here

        // N_1 = Mu_Ea, N_2 = Sigma_Ea^2, N_3 = Mu_v, N_4 = site concentration

        // store the variables on the heap
        gml::variable *site = new gml::variable("site", N);
        gml::variable *Ea = new gml::variable("Ea", N);
        gml::variable *v = new gml::variable("v", N);

        for(int i=0; i<N; ++i){
            Ea->canonical_column[i]->mu.resize(1,1);
            Ea->canonical_column[i]->mu(0) = kinetic_parameters[i];
            Ea->canonical_column[i]->Sigma.resize(1,1);
            Ea->canonical_column[i]->Sigma(0,0) = kinetic_parameters[N+i];

            v->canonical_column[i]->mu.resize(1,1);
            v->canonical_column[i]->mu(0) = kinetic_parameters[2*N+i];

            site->canonical_column[i]->g = log(kinetic_parameters[3*N+i]);
        }
        // need to store these factors on the heap, not the stack
        gml::factor *phi_1 = new gml::factor("site", site);
        gml::factor *phi_2 = new gml::factor("Ea", Ea);
        gml::factor *phi_3 = new gml::factor("v", v);

        vector<gml::factor*> factors = {
                phi_1, phi_2, phi_3
        };

        gml::probabilistic_graph surface_system (factors);

        // store model
        model = surface_system;
    }

    /**
     * @brief: Save the model fit to a text file.
     * @param name : name of fit file
     */

    void report_model( string name ){

        // compute the signal based on the simulated reaction rates
        composite_signal s_model = kinetic_signal(&model, T_i, T_f, ramp);// this isn't an efficient way to do this but leave it for now
        s_model.compute_pressures();
        s_model.interpolate_total_pressure();

        // interpolate the reaction rates
        vector<ublas::vector<double>> interpolated_signals (s_model.partial_pressures.size());
        // interpolate to align all rows of output
        for(int i=0; i<interpolated_signals.size(); ++i){
            interpolated_signals[i] = analysis::lin_interp(s_model.times[i], s_model.partial_pressures[i], time);
        }

        ublas::vector<double> interpolated_sum = analysis::lin_interp(s_model.times[0], s_model.total_pressure, time);

        // write to file
        ofstream outfile;
        outfile.open(name);

        outfile << "Kinetic model fit parameters and interpolated data.\n";
        outfile << "Site specific parameters:\n";
        for(int i=0; i<model.variable_map["site"]->canonical_column.size(); ++i){
            outfile << "Site " << i << ": \n";
            outfile << "Concentration: ," << exp(model.variable_map["site"]->canonical_column[i]->g) << endl;
            outfile << "Mean activation energy: ," << model.variable_map["Ea"]->canonical_column[i]->mu(0) << endl;
            outfile << "Activation energy variance: ," << model.variable_map["Ea"]->canonical_column[i]->Sigma(0,0) << endl;
            outfile << "Mean vibrational frequency: ," << model.variable_map["v"]->canonical_column[i]->mu(0) << endl;
            outfile << endl;
        }

        outfile << "time,temperature,signal,";
        for(int i=0; i<interpolated_signals.size(); ++i){
            outfile << "site" << i << ",";
        }
        outfile << "model_composite" << endl;

        for(int i=0; i<time.size(); ++i){
            outfile << time[i] << ",";
            outfile << temperature[i] << ",";
            outfile << signal[i] << ",";
            for(int j=0; j<interpolated_signals.size(); ++j){
                outfile << interpolated_signals[j][i] << ",";
            }
            outfile << interpolated_sum[i] << endl;
        }

        outfile.close();
    }
};





