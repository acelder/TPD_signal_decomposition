//
// Created by Alexander Elder on 9/2/17.
//

#include "genetic_algorithm.h"

individual::individual(double *params, int num_params, double fitness) : fitness(fitness) {
    n_parameters = num_params;
    int i;
    parameters.resize(n_parameters);
    for(i=0; i<n_parameters; i++){
        parameters[i] = params[i];
    }

}

individual::individual(int n_parameters) : n_parameters(n_parameters) {
    parameters.resize(n_parameters);
}



generation::generation(int size) : size(size) {
    members.resize(size);

}
