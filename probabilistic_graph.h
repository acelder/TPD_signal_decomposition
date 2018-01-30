//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H
#define TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H

#include "gaussian_operations.h"
#include "factor_operations.h"
#include <vector>
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;

namespace gml{
    class probabilistic_graph{
        vector<variable*>variables;

        vector<factor*> factors;
    public:
        unordered_map<string, variable*> variable_map;
        unordered_map<string, factor*> factor_map;

        // constructor
        probabilistic_graph() = default;
        probabilistic_graph( vector<factor*> factor_list );

        // get attribute functions
        int get_num_variables();
        variable *get_variable(int variable_index);
        variable *get_variable(string variable_name);

        int get_num_factors();
        factor *get_factor(int factor_index);
        factor *get_factor(string factor_name);

        // member functions
        void observe_evidence(variable *target_variable, int var_assign);


    };

    probabilistic_graph::probabilistic_graph( vector<factor*> factor_list ){
        for(auto fac : factor_list){
            // add factor to graph factor list
            factors.push_back(fac);
            factor_map[fac->name] = fac;

            // add variables to variable list if not already stored
            for(auto var : fac->scope){
                if(variable_map.find(var->name) == variable_map.end()){
                    variables.push_back(var);
                    variable_map[var->name] = var;
                }
            }
        }
    }

    int probabilistic_graph::get_num_variables(){
        return variables.size();
    }

    variable *probabilistic_graph::get_variable(int variable_index) {
        return variables[variable_index];
    }

    variable *probabilistic_graph::get_variable(string variable_name) {
        return variable_map[variable_name];
    }

    int probabilistic_graph::get_num_factors() {
        return factors.size();
    }

    factor *probabilistic_graph::get_factor(int factor_index){
        return factors[factor_index];
    }

    factor *probabilistic_graph::get_factor(string factor_name){
        return factor_map[factor_name];
    }

    void probabilistic_graph::observe_evidence(variable *target_variable, int var_assign) {
        for(int i=0; i<factors.size(); ++i){
            // if the target variable is in a given factors scope, perform observe_evidence operation
            if(factors[i]->scope.find(target_variable) != factors[i]->scope.end()){
                factors[i]->observe_evidence(target_variable, var_assign);
            }
        }
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_PROBABILISTIC_GRAPH_H
