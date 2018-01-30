//
// Created by Alexander Elder on 1/28/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_INTERPOLATION_H
#define TPD_SIGNAL_DECOMPOSITION_INTERPOLATION_H

#include <map>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;

namespace lin_alg{
    namespace analysis{
        ublas::vector<double> lin_interp(ublas::vector<double> x, ublas::vector<double> f_x, ublas::vector<double> x_q);
        double lagrange_polynomial(ublas::vector<double> x, double x_q, int j);
        ublas::vector<double> lagrange_interp(ublas::vector<double> x, ublas::vector<double> f_x, ublas::vector<double> x_q);


        /**
     * \brief: compute the 1-D linear interpolation of f(x) at points x_q
     * @param x : support points
     * @param f_x : support values
     * @param x_q : query points
     * @return : f_x_q query values
     * Note: query points don't need to be ordered
     */
        ublas::vector<double> lin_interp(ublas::vector<double> x, ublas::vector<double> fx, ublas::vector<double> x_q){
            //initialize data map and interpolated value vector
            map<double, double> data;
            ublas::vector<double> fx_q (x_q.size());

            for(int j=0; j<x.size(); ++j){
                data[x[j]] = fx[j];
            }

            typedef map<double, double>::iterator it;
            it i;

            for(int j=0; j<x_q.size(); ++j){
                i = data.upper_bound(x_q[i]);
                // edge cases
                if(i==data.end()){
                    fx_q[i] = (--i)->second;
                }
                else if(i==data.begin()){
                    fx_q[i] = i->second;
                }
                else{
                    it l = i; --l;
                    const double delta = (x_q[j] - l->first)/(i->first - l->first);
                    fx_q[j] = delta*i->second + (1-delta)*l->second;
                }
            }
            return fx_q;
        }

    /**
     * \brief: compute the lagrange polynomial of x for support point j
     * @param x : support vector
     * @param f_x : support values
     * @param x_q : query value
     * @param j :
     * @return L_j(x_q)
     * NEEDS TO BE TESTED
     */
        double lagrange_polynomial(ublas::vector<double> x, double x_q, int j){
            double L_j = 1; int N = x.size();

            for(int k=0; k<N; ++k){
                if(k != j ){
                    L_j *= (x_q - x[k])/(x[j] - x[k]);
                }
            }

            return L_j;
        }

    /**
     * \brief: compute the lagrange interpolation of the query points x+q for support vector x and values fx
     * @param x : support points
     * @param fx : support values
     * @param x_q : query values
     * @return : Lagrange interpolants of x_q
     * NEEDS TO BE TESTED
     */
    ublas::vector<double> lagrange_interp(ublas::vector<double> x, ublas::vector<double> fx, ublas::vector<double> x_q){
            // initialize vector of polynomial values
            ublas::vector<double> p (x_q.size());
            // compute polynomial values, complexity: O(x_q.size() * fx.size())
            for(int i=0; i<p.size(); ++i){
                for(int l=0; l<fx.size(); ++l){
                    p[i] += fx[l] * lagrange_polynomial(x, x_q[i], l);
                }
            }
            return p;
        }
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_INTERPOLATION_H
