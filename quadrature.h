//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H
#define TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;

namespace lin_alg{
    namespace analysis{

        /**
         * \brief: composite newton-cotes numerical integration
         * @param function
         * @param n_supports
         * @param range
         * @return
         */
        double newton_cotes(double (*function)(double), int n_domains, ublas::vector<double> range, string rule){
            double s = 0;
            double a = range(0), b = range(1);
            double step_size = (b-a)/n_domains;

            if(rule == "trap"){
                for(int i=0; i<n_domains-1; ++i){
                    double x_0 = a + i*step_size;
                    double x_1 = a + (i+1)*step_size;
                    double f_0 = function(x_0);
                    double f_1 = function(x_1);

                    s += (x_1 - x_0)*(f_1 + f_0)/2;
                }
            }
            else if(rule == "simpson"){

                for(int i=0; i<n_domains-1; ++i){
                    double x_0 = a + i*step_size;
                    double x_1 = a + (i+1)*step_size/2;
                    double x_2 = a + (i+1)*step_size;
                    double f_0 = function(x_0);
                    double f_1 = function(x_1);
                    double f_2 = function(x_2);

                    s += (x_1 - x_0)*(f_2 + 4*f_1 + f_0)/6;
                }
            }
            else{
                throw runtime_error("Unknown leading order method.");
            }

            return s;
        }

    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H
