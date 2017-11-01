//
// Created by Alexander Elder on 9/2/17.
//

#ifndef TPD_DECONVOLUTION_GENETIC_ALGORITHM_H
#define TPD_DECONVOLUTION_GENETIC_ALGORITHM_H

#include <vector>

/*
 * Class for individuals in the population
 */

class individual{
public:
    std::vector<double> parameters;
    int n_parameters;

    double fitness;

    bool feasible;
    double constraint_violation;

    individual(int n_parameters);

    individual(double *params, int num_params, double fitness);

};

/*
 * Class encoding a generation of the population
 */

class generation{
public:
    int size;
    std::vector<individual> members;

    generation(int size);

};

#endif //TPD_DECONVOLUTION_GENETIC_ALGORITHM_H
