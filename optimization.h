//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H
#define TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace ublas = boost::numeric::ublas;


namespace opt_lib{
    const double eps = 5e-16;

    struct optimum{
        ublas::vector<double> arg_max;
        ublas::vector<double> max;
    };

    struct learning_parameters{
        double rate;
        long num_iter;
        double convergence_criterion;
    };

    optimum gradient_descent( ublas::vector<double> (*function)(ublas::vector<double>),
                              ublas::vector<double> initial_guess, learning_parameters params ){
        long N = initial_guess.size();
        ublas::vector<double> x = initial_guess;
        ublas::vector<double> e (N);
        for(int i=0; i<N; ++i){
            e[i] = eps;
        }

        for(int i=0; i<params.num_iter; ++i){
            ublas::vector<double> forward_difference = x + eps;
            ublas::vector<double> backward_difference= x - eps;
            ublas::vector<double> forward_eval = function(forward_difference);
            ublas::vector<double> backward_eval = function(backward_difference);
            ublas::vector<double> central_difference = forward_eval - backward_eval;
            for(int j=0; j<central_difference.size(); ++j){
                central_difference[j] = central_difference[j]/(2*eps);
            }
            ublas::vector x_prev = x;
            ublas::vector x = x - params.rate*central_difference;
            if(ublas::norm_1(x) - ublas::norm_1(x_prev) <= params.convergence_criterion){
                break;
            }
        }

        optimum opt;
        opt.arg_max = x;
        opt.max = function(x);
        return opt;
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H
