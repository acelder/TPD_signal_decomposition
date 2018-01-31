//
// Created by Alexander Elder on 11/3/17.
//

#include "deprecated/finite_difference_simulations.h"

surface_system::surface_system(const external_parameters &p) : p(p) {}

void surface_system::operator()(const boost_vector &x, boost_vector &dxdt, double t){
    /*
     * Compute the time derivative of the external surface coverage at time t.
     * @param    x: Vector of length 1 containing the surface concentration of adsorbate in units molecules/cm^2
     * @param dxdt: Vector of length 1 to which the derivative will be assigned
     * @param    t: time
     */

    double T_0     = p.initial_temperature;
    double alpha   = p.ramp_rate;
    double Ea      = p.Ea;
    double nu      = p.nu;

    const double R = 8.314;
    double T       = T_0 + alpha*t;
    double k_des   = nu*exp(-Ea/(R*T));

    dxdt[0] = x[0]*k_des;
}