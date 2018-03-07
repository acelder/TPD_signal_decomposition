//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H
#define TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <functional>
#include <cstdlib>


namespace ublas = boost::numeric::ublas;
using namespace std;

namespace opt_lib{
    const double eps = 5e-16;

    struct optimum{
        ublas::vector<double> arg_opt;
        double opt;
    };

    struct learning_parameters{
        double rate;
        long num_iter;
        double convergence_criterion;

        learning_parameters() = default;
        learning_parameters(double rate, long num_iter, double convergence_criterion) :
                rate(rate), num_iter(num_iter), convergence_criterion(convergence_criterion) {}
    };

    optimum gradient_descent( function<double(ublas::vector<double>)> f,
                              ublas::vector<double> initial_guess, learning_parameters params ){
        long N = initial_guess.size();
        ublas::vector<double> x = initial_guess; // only supports ublas::vector pass and return types at this time
        ublas::vector<double> x_prev;


        for(int i=0; i<params.num_iter; ++i){
            ublas::vector<double> central_difference (N, 0);
            // compute the central difference approximation of the gradient
            for(int j=0; j<N; ++j){
                ublas::vector<double> e (N, 0);
                e(j) = eps;
                ublas::vector<double> forward_difference = x + e;
                ublas::vector<double> backward_difference = x - e;
                double forward_eval = f(forward_difference);
                double backward_eval= f(backward_difference);
                central_difference(j) = (forward_eval - backward_eval)/(2*eps);

            }

            x_prev = x;
            x = x - params.rate*central_difference;// subtract the gradient to descend the valley
            if(ublas::norm_1(x - x_prev) <= params.convergence_criterion){
                break;
            }
        }

        optimum opt;
        opt.arg_opt = x;
        opt.opt = f(x);
        return opt;
    }

    optimum gradient_ascent( function<double(ublas::vector<double>)> f,
                             ublas::vector<double> initial_guess, learning_parameters params ){
        long N = initial_guess.size();
        ublas::vector<double> x = initial_guess; // only supports ublas::vector pass and return types at this time
        ublas::vector<double> x_prev;


        for(int i=0; i<params.num_iter; ++i){// gradient descent needs to be updated
            ublas::vector<double> central_difference (N, 0);
            // compute the central difference approximation of the gradient
            for(int j=0; j<N; ++j){
                ublas::vector<double> e (N, 0);
                e(j) = eps;
                ublas::vector<double> forward_difference = x + e;
                ublas::vector<double> backward_difference = x - e;
                double forward_eval = f(forward_difference);
                double backward_eval= f(backward_difference);
                central_difference(j) = (forward_eval - backward_eval)/(2*eps);

            }

            x_prev = x;
            x = x + params.rate*central_difference;// add the gradient to climb the hill
            if( ublas::norm_1(x - x_prev) <= params.convergence_criterion){
                break;
            }
        }

        optimum opt;
        opt.arg_opt = x;
        opt.opt = f(x);
        return opt;
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_OPTIMIZATION_H
