

#include "unit_tests.h"
#include "tpd_sim.h"
#include "probabilistic_graph.h"
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <fstream>
#include "quadrature.h"


using namespace std;

int main( int argc , char **argv ) {

    clock_t t;
    t = clock();

    vector<double> site_dist = {1.0, 1.0};
    vector<double> activation_energies = {5.0e4, 5.3e4};
    vector<double> prefactors = {1.0e13, 1.0e13};
    double T_i = 120;
    double T_f = 300;
    double ramp = 0.25;
    composite_signal s = kinetic_signal(site_dist, activation_energies, prefactors, T_i, T_f, ramp);

    ublas::vector<double> T = s.compute_temperature(0);
    s.interpolate_rates();
    s.compute_pressures();
    s.interpolate_total_pressure();

    // report test
    ofstream outfile;
    outfile.open("2_site_run.csv");
    outfile << "time, temperature, rate_0, rate_1, total_pressure\n";

    for(int i=0; i<s.times[0].size(); ++i){
        outfile << s.times[0][i] << "," << T[i] << "," << s.partial_pressures[0][i] << ","
                << s.partial_pressures[1][i] << "," << s.total_pressure[i] << endl;
    }

    outfile.close();
    t = clock() - t;
    cout << "Process finished in " << (float)t/CLOCKS_PER_SEC << " seconds." << endl;

    return 0;
}