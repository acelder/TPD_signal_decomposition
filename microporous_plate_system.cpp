//
// Created by Alexander Elder on 11/1/17.
//

#include "finite_difference_simulations.h"

// define pi
const double pi = boost::math::constants::pi<double>();

/*
 * constructor for microporous plate system objects
 */

microporous_plate_system::microporous_plate_system(const microporous_parameters &p) : p(p) {}

/*
 * Function wrapper for microporous plate system.
 */
void microporous_plate_system::operator()(const boost_vector &x, boost_vector &dxdt, double t) {
    double rho   = p.rho; // site density
    double nu    = p.nu; // pre-exponential (bond vibrational frequency)
    double Ea    = p.Ea; // activation energy of desorption
    double s     = p.steric_factor;
    double am    = p.adsorbate_mass;
    double sm    = p.site_mass;
    double ar    = p.adsorbate_radius;
    double sr    = p.site_radius;
    double D_0   = p.D_0; // diffusion pre-exponential
    double Ea_D  = p.Ea_D; // diffusion activation energy
    double T_0   = p.initial_temp;// initial temperature
    double alpha = p.ramp_rate;// ramp rate
    double H     = p.characteristic_length;// plate thickness

    // kinetic readsorption coefficients
    double R     = 8.314;
    double kb    = 1.38064852 * pow(10,-23);
    double mu    = am * sm / (am + sm); // reduced mass
    double d     = ar + sr;
    double T     = T_0 + alpha*t;
    double k_rev = s * pi * pow(d,2) * (8 * kb * T)/pow((pi * mu),0.5);

    int n_nodes  = x.size()/2; // extract bin count
    int i;

    // compute temperature dependent coefficients


    // compute change in surface concentration
    for(i=0; i<(n_nodes); ++i){
        dxdt[i] = k_rev * x[i+n_nodes]* (rho - x[i]) // second order readsorption term
                  - x[i] * nu * exp(-Ea / (R * T));  // first order desorption term
    }


    // compute mass transfer coefficients
    double D     = D_0 * exp(-Ea_D / (R * T)); // activated diffucion coefficient
    double dz    = H/n_nodes;                 // finite z-axis differential
    double sigma = D/pow(dz,2);               // differential-scaled diffusion coefficient
    // compute the change in internal pore concentration
    int i_start = n_nodes;// may cause off by one!
    int i_end   = x.size()-1;

    //   start with edge cases
    dxdt[i_start] = -2 * sigma * x[i_start] + sigma * x[i_start+1] // mass transfer component, no flux boundary condition
                    - dxdt[i_start-n_nodes]; // reactive component
    dxdt[i_end]   = -2 * sigma * x[i_end]   + sigma * x[i_end-1] // vacuum boundary condition
                    - dxdt[i_end-n_nodes];

    for(i=i_start+1; i<i_end; i++){
        dxdt[i] = sigma * (x[i-1] + x[i+1]) - 2 * sigma * x[i] - dxdt[i-n_nodes];
    }
}

