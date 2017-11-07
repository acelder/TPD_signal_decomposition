//
// Created by Alexander Elder on 11/3/17.
//


#include "finite_difference_simulations.h"

const double pi = boost::math::constants::pi<double>();


void mesoporous_film_system::operator()(const boost_vector &x, boost_vector &dxdt, double t) {
    //          re-name parameters
    double rho   = p.rho; // site density
    double nu    = p.nu; // pre-exponential (bond vibrational frequency)
    double Ea    = p.Ea; // activation energy of desorption
    double s     = p.steric_factor;
    double am    = p.adsorbate_mass; // expressed as kg
    double sm    = p.site_mass;      // also kg
    double ar    = p.adsorbate_radius;// meteres
    double sr    = p.site_radius;     // meters
    double d_pore= p.pore_diameter;    // meters
    double T_0   = p.initial_temp;// initial temperature
    double alpha = p.ramp_rate;// ramp rate
    double H     = p.characteristic_length;// plate thickness

    // compute kinetic readsorption coefficients
    double R     = 8.314;
    double k_b    = 1.38064852 * pow(10,-23);
    double mu    = am * sm / (am + sm); // reduced mass
    double d     = ar + sr;
    double T     = T_0 + alpha*t;
    double k_rev = s * pi * pow(d,2) * (8 * k_b * T)/pow((pi * mu),0.5);

    int n_nodes  = (x.size() - 12)/2; // extract bin count
    int i;

    //          compute temperature dependent diffusion coefficients


    // set the first 13 bins of dxdt to zero
    for(i=0; i<13; ++i){
        dxdt[i] = 0;
    }

    // compute change in surface concentration
    for(i=0; i<(n_nodes); ++i){
        dxdt[i] = k_rev * x[i+n_nodes]* (rho - x[i]) // second order readsorption term
                  - x[i] * nu * exp(-Ea / (R * T));  // first order desorption term
    }


    // compute mass transfer coefficients, in the mesorporous system the mass transfer coefficient
    // may vary as a function of concentration
    double mts   = sqrt(8*k_b*T/(pi * am));          // m/s
    vector<double> D (n_nodes);                      // diffusion coefficient
    vector<double> l (n_nodes);                      // mean free path (meters)

    for(i=0; i<n_nodes; ++i){
        l[i] = 1/(pi*pow(2*ar,2)*x[i]);              // mean free path = 1/(pi*adsorbate_diameter^2*concentration)
    }

    for(i=0; i<n_nodes; ++i){                        // initialize the diffusion coefficient according to kinetic theory
        if(d_pore < l[i]){                           // or knudsen diffusion
            D[i] = d_pore*mts/3;
        }
        else{
            D[i] = l[i]*mts/3;
        }
    }


    double dz    = H/n_nodes;                        // finite z-axis differential
    vector<double> sigma (n_nodes);                  // differential-scaled diffusion coefficient
    for(i=0; i<n_nodes; i++){
        sigma[i] =  D[i]/pow(dz,2);
    }
    // compute the change in internal pore concentration
    int i_start = n_nodes;// may cause off by one!
    int i_end   = x.size()-1;

    //   start with edge cases
    dxdt[i_start] = -2 * sigma[i_start] * x[i_start] + sigma[i_start+1] * x[i_start+1] // mass transfer component, vacuum boundary condition
                    - dxdt[i_start-n_nodes]; // reactive component
    dxdt[i_end]   = -2 * sigma[i_end] * x[i_end]   + sigma[i_end-1] * x[i_end-1]  // vacuum boundary condition
                    - dxdt[i_end-n_nodes];

    for(i=i_start+1; i<i_end; i++){
        dxdt[i] = sigma[i-1] * x[i-1] + sigma[i+1] * x[i+1] - 2 * sigma[i] * x[i] - dxdt[i-n_nodes];
    }

}