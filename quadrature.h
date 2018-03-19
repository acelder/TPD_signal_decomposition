//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H
#define TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H

#include "storage_adaptors.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <functional>
#include <vector>

using namespace std;
namespace ublas = boost::numeric::ublas;

namespace lin_alg{
    namespace analysis {

        /**
         * @brief: composite newton-cotes numerical integration
         * @param function: function to be integrated
         * @param n_supports: number of support domains
         * @param range: range of integration
         * @return
         */
        double newton_cotes(double (*function)(double), int n_domains, ublas::vector<double> range, string rule = "trap") {
            double s = 0;
            double a = range(0), b = range(1);
            double step_size = (b - a) / n_domains;

            if (rule == "trap") {
                for (int i = 0; i < n_domains; ++i) {
                    double x_0 = a + i * step_size;
                    double x_1 = a + (i + 1) * step_size;
                    double f_0 = function(x_0);
                    double f_1 = function(x_1);

                    s += (x_1 - x_0) * (f_1 + f_0) / 2;
                }
            } else if (rule == "simpson") {

                for (int i = 0; i < n_domains; ++i) {
                    double x_0 = a + i * step_size;
                    double x_2 = a + (i + 1) * step_size;
                    double x_1 = (x_2 + x_0)/2;
                    double f_0 = function(x_0);
                    double f_1 = function(x_1);
                    double f_2 = function(x_2);

                    s += (x_2 - x_0) * (f_2 + 4 * f_1 + f_0) / 6;
                }
            } else {
                throw runtime_error("Unknown leading order method.");
            }

            return s;
        }
        /**
         * @brief Two-dimensional newton-cotes numerical integration. Add flexiblity in intermediate points later.
         * @param f : function to be integrated
         * @param n_domains : ublas::vector<double> number of composite domains in each dimension
         * @param range : ublas::matrix<double> Range of integration in each dimension
         * @param rule : integration rule
         * @return
         */



        double dbl_newton_cotes( function<double(ublas::vector<double>)> f, ublas::vector<long> &n_domains,
                                ublas::matrix<double> &range, string rule = "trap") {
            double s = 0;
            double a0 = range(0, 0);
            double b0 = range(0, 1);
            double n0 = n_domains(0);
            double a1 = range(1, 0);
            double b1 = range(1, 1);
            double n1 = n_domains(1);

            double step_size1 = (b0 - a0) / n0;
            double step_size2 = (b1 - a1) / n1;

            if (rule == "trap") {
                for (int i = 0; i < n0; ++i) {
                    for (int j = 0; j < n1; ++j) {
                        // composite trapezoidal rule for 0 intermediate points
                        ublas::vector<double> r00(2);
                        ublas::vector<double> r01(2);
                        ublas::vector<double> r10(2);
                        ublas::vector<double> r11(2);

                        r00(0) = a0 + i * step_size1;
                        r00(1) = a1 + j * step_size2;

                        r01(0) = a0 + i * step_size1;
                        r01(1) = a1 + (j + 1) * step_size2;

                        r10(0) = a0 + (i + 1) * step_size1;
                        r10(1) = a1 + j * step_size2;

                        r11(0) = a0 + (i + 1) * step_size1;
                        r11(1) = a1 + (j + 1) * step_size2;

                        double f00 = f(r00);
                        double f01 = f(r01);
                        double f10 = f(r10);
                        double f11 = f(r11);

                        double h = (r10(0) - r00(0));
                        double k = (r01(1) - r00(1));

                        s += 0.25 * h * k * (f00 + f01 + f10 + f11);
                    }
                }
            }
            else if(rule == "simpson"){
                // composite simpson's rule for 5 intermediate points

                // Simpson's weight matrix
                const double W [5][5] = {
                        {1, 4, 2, 4, 1},
                        {4, 16, 8, 16, 4},
                        {2, 8, 4, 8, 2},
                        {4, 16, 8, 16, 4},
                        {1, 4, 2, 4, 1}
                };

                double H = step_size1/4;
                double K = step_size2/4;

                for(int i=0; i<n0; ++i){
                    for(int j=0; j<n1; ++j){
                        double I = 0;
                        double x = a0 + i * step_size1;
                        double y = a1 + j * step_size2;



                        for(int k=0; k<5; ++k){

                            double x_k = x + (k)*H;


                            for(int l=0; l<5; ++l){
                                double y_l = y + (l)*K;

                                ublas::vector<double> r (2);
                                r(0) = x_k; r(1) = y_l;
                                double f_val = f(r);
                                I += W[k][l]*f_val;
                            }
                        }

                        s += H*K*I/9;

                    }
                }
            }
            return s;
        }

        /**
         * @brief Multiple Guassian quadrature over 1 or 2 dimensions-not yet implementet
         * @param function
         * @param n_domains
         * @param range
         * @return
         */


        double gaussian_composite(function<double(ublas::vector<double>)> f, ublas::vector<long> &n_domains,
                             ublas::matrix<double> &range, long order){


            // weights and positions for order 1 Gaussian quadrature
            std::vector<double> w_1 = {2.0};
            std::vector<double> zeta_1 = {0.0};

            // weights and positions for order 2 Gaussian quadrature
            std::vector<double> w_2 = {1.0, 1.0};
            std::vector<double> zeta_2 = { -1.0/sqrt(3.0), 1.0/sqrt(3.0) };

            // weights and positions for order 3 Gaussian quadrature
            std::vector<double> w_3 = {5.0/9.0, 8.0/9.0, 5.0/9.0};
            std::vector<double> zeta_3 = { -sqrt(3.0)/sqrt(5.0), 0.0, sqrt(3.0)/sqrt(5.0)};

            // weights and positions for order 4 Gaussian quadrature
            std::vector<double> w_4 = {
                    (18.0 - sqrt(30.0))/36.0,
                    (18.0 + sqrt(30.0))/36.0,
                    (18.0 + sqrt(30.0))/36.0,
                    (18.0 - sqrt(30.0))/36.0
            };

            std::vector<double> zeta_4 = {
                    -sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)),
                    -sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
                    sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)),
                    sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))
            };

            std::vector<std::vector<double>> w = {w_1, w_2, w_3, w_4};
            std::vector<std::vector<double>> zeta = {zeta_1, zeta_2, zeta_3, zeta_4};

            if(order > 4){
                cout << "Warning, Gaussian quadrature not supported beyond order 4. Applying order 4 quadrature rule.\n";
                order = 4;
            }
            if(order < 1){
                order = 1;
            }

            long N = n_domains(0);
            long M = n_domains(1);
            double composite_integral = 0, local_integral;

            ublas::vector<double> a (2);
            ublas::vector<double> b (2);

            double H = (range(0,1) - range(0,0))/n_domains(0);
            double K = 0;
            if(M > 0){
                K = (range(1,1) - range(1,0))/n_domains(1);
            }
            ublas::vector<double> r (2,0);


            // iterate over composite domains, mapping each into [-1, 1]^n and approximating the integral to the given order
            for(long i=0; i<N; ++i){
                for(long j=0; j<M; ++j){
                    local_integral = 0;

                    if(M <= 1){
                        a(0) = range(0,0) + i*H;
                        b(0) = a(0) + H;


                        for(long k=0; k<order; ++k){
                            r(0) = a(0) + (b(0) - a(0))/2 * (1 + zeta[order-1][k]);
                            local_integral += w[order-1][k] * f(r);
                        }
                        local_integral *= (b(0) - a(0))/2.0;
                        composite_integral += local_integral;
                    }
                    else{
                        a(0) = range(0,0) + i*H;
                        b(0) = a(0) + H;
                        a(1) = range(1,0) + i*K;
                        b(1) = a(1) + K;

                        for(long k=0; k<order; ++k){
                            r(0) = a(0) + (b(0) - a(0))/2 * (1 + zeta[order-1][k]);
                            for(long l=0; l<order; ++l){
                                r(1) = a(1) + (b(1) - a(1))/2 * (1 + zeta[order-1][l]);
                                local_integral += w[order-1][k]*w[order-1][l]*f(r);
//                                if(i==0 and j==0){
//                                    cout << "Order "<< order << endl;
//                                    cout << "w_" << k << " w_" << l << " = " << w[order-1][k] << "*" << w[order-1][l] << "  ";
//                                    cout << "zeta_" << k << ", zeta_" << l << " = " << zeta[order-1][k] << ", " << zeta[order-1][k] << "  \n";
//                                }
                            }


                        }
                        local_integral *= (b(0) - a(0))*(b(1) - a(1))/4;
                        composite_integral += local_integral;
                    }
                }
            }
            return composite_integral;
        }

    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_QUADRATURE_H
