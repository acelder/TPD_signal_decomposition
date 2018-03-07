//
// Created by Alexander Elder on 2/15/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_SOLVERS_H
#define TPD_SIGNAL_DECOMPOSITION_SOLVERS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <string>
#include "lin_alg.h"

using namespace std;
namespace ublas = boost::numeric::ublas;

namespace lin_alg{
    namespace analysis{


        struct solver_opt{
            bool approximate;
            bool reduced_step;
            double rate = 1;
            double threshold;
            long max_iter;
            ublas::matrix<double> J;
        };

        ublas::vector<double> fsolve(ublas::vector<double> (*f)(ublas::vector<double>),
                                     ublas::vector<double> initial, solver_opt options){
            double threshold = options.threshold;
            long dim = initial.size();
            const double eps = 5e-16;
            double alpha = options.rate;

            // initialize identity matrix and state vector
            ublas::vector<double> x = initial;

            // if approximate method, use Broyden rank-one update rule
            if(options.approximate){
                ublas::vector<double> fx = f(x);
                long f_dim;

                // compute initial estimate of B
                ublas::matrix<double> B (f_dim, dim);
                for(int i=0; i<f_dim; ++i){
                    for(int j=0; j<dim; ++j){
                        ublas::vector<double> e (dim, 0.0);
                        e[j] = eps;
                        ublas::vector<double> inc = x + e;
                        ublas::vector<double> f_diff = f(inc) - f(x);
                        B(i,j) = f_diff(i)/eps;
                    }
                }
                // enter main loop, terminate at max iterations or when threshold is achieved
                for(int i=0; i<options.max_iter; ++i){
                    LUP_decomposition LU (B);
                    ublas::vector<double> pf = ublas::prod(LU.P, -fx);
                    ublas::vector<double> c = substitution(LU.L, pf, false);// probably won't compile
                    ublas::vector<double> dx = substitution(LU.U, c, true);
                    if(ublas::norm_1(dx) < options.threshold){
                        break;
                    }
                    else{
                        fx = f(x + alpha*dx);
                        ublas::matrix<double> O = ublas::outer_prod(fx, x);
                        B = B + O / pow(ublas::norm_1(x), 2.0);// update Broyden matrix
                        x = x + alpha*dx;
                    }
                }

            }
            else{
                ublas::vector<double> fx = f(x);
                ublas::matrix<double> J = options.J;
                for(int i=0; i<options.max_iter; ++i){
                    LUP_decomposition LU (J);
                    ublas::vector<double> pf = ublas::prod(LU.P, -fx);
                    ublas::vector<double> c = substitution(LU.L, pf, false);
                    ublas::vector<double> dx = substitution(LU.U, c, true);
                    if(ublas::norm_1(dx) < options.threshold){
                        break;
                    }
                    else{
                        x = x + alpha*dx;
                        fx = f(x);
                    }
                }
            }
            return x;
        }
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_SOLVERS_H
