//
// Created by Alexander Elder on 2/26/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_SENSITIVITY_ANALYSIS_H
#define TPD_SIGNAL_DECOMPOSITION_SENSITIVITY_ANALYSIS_H

#include "tpd_sim.h"
#include "interpolation.h"
#include "probabilistic_graph.h"
#include <string>
#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

/**
 * @brief Analyze the sensitivity of the signal with respect to the activation energy skew parameter.
 */

void skew_sensitivity_analysis(string filename){
    cout << "Performing skew sensitivity analysis. File will be saved to: \n";
    cout << filename << endl;

    // set up variabls
    gml::variable site ("site", 1);
    gml::variable v ("v", 1);
    gml::variable Ea ("Ea", 1);

    gml::canonical_form site1;
    site1.g = log(1);

    gml::canonical_form v1;
    v1.mu.resize(1,1);
    v1.mu(0) = 1e13;
    v1.Sigma.resize(1,1);
    v1.Sigma(0,0) = 1e2;
    v1.scope.insert("Hz");
    v1.compute_canonical_parameters();

    gml::skew_canonical_form Ea1;
    Ea1.mu.resize(1,1);
    Ea1.mu(0) = 5e4;
    Ea1.Sigma.resize(1,1);
    Ea1.Sigma(0,0) = 1e4;
    Ea1.scope.insert("energy");
    Ea1.skew.resize(1);
    Ea1.compute_canonical_parameters();

    site.canonical_column[0] = &site1;
    v.canonical_column[0] = &v1;
    Ea.canonical_column[0] = &Ea1;

    gml::factor phi1 ("site", &site);
    gml::factor phi2 ("Ea", &Ea);
    gml::factor phi3 ("v", &v);

    vector<gml::factor*> factors = {
            &phi1, &phi2, &phi3
    };

    // initialize graphical model
    gml::probabilistic_graph model(factors);


    vector<composite_signal> signals (10);
    vector<double> skews = {
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0
    };

    for(int i=0; i<10; ++i){
        cout << "Running test " << i << ".\n";
        Ea1.skew[0] = skews[i];
        signals[i] = kinetic_signal(&model, 120, 300, 0.25);
        signals[i].compute_pressures();
        cout << "Test complete.\n";
    }

    vector<ublas::vector<double>> interpolated_pressures (9);

    // Interpolate computed pressures to allow direct comparisons between simulations.
    for(int i=1; i<10; ++i){
        interpolated_pressures[i] = lin_alg::analysis::lin_interp(signals[i].times[0], signals[i].total_pressure, signals[0].times[0]);
    }

    cout << "Tests complete. Writing to file.\n";
    // Write the interpolated pressures to a .csv file with the input name.
    fstream outfile;
    outfile.open(filename);

    for(int i=0; i<signals[0].times.size(); ++i){
        outfile << signals[0].times[i] << "," << signals[0].total_pressure[i] << ",";

        for(int j=1; j<interpolated_pressures.size(); ++j){
            outfile << interpolated_pressures[j][i] << ",";
        }

        outfile << endl;
    }


}

#endif //TPD_SIGNAL_DECOMPOSITION_SENSITIVITY_ANALYSIS_H
