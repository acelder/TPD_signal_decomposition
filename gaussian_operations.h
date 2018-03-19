//
// Created by Alexander Elder on 1/26/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H
#define TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H

#include "lin_alg.h"
#include <string>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace lin_alg;
namespace ublas = boost::numeric::ublas;



/**
 * Gaussian distribution expressed as a canonical form, also stores mean and variance.
 */

namespace gml {

    const double pi =  3.1415926535897;
    /*
     * Declarations
     */
    // variable template
    template<class T> struct variable{
        string name;                   // name of the variable
        int cardinality;               // number of parameters
        T assignment;                  // int for discrete variables, double for continuous

        variable() = default;
        variable(string name, int card) : name(name), cardinality(card) {}
        variable(string name, int card, T assign) : name(name), cardinality(card), assignment(assign){}

    };

    struct canonical_form;

    canonical_form vacuous_form(double dimension);
    canonical_form canonical_product( const canonical_form &form1, const canonical_form &form2);
    canonical_form canonical_quotient(const canonical_form &numerator, const canonical_form &denominator);
    canonical_form canonical_marginal(const canonical_form &joint_form, const vector<variable<double>*> &vars);
    canonical_form canonical_reduction(const canonical_form &joint_form, const vector<variable<double>*> &vars,
                                       const ublas::vector<double> &y);
    canonical_form canonical_collapse(const canonical_form &form_A, const canonical_form &form_B);

    // implementations

    struct canonical_form{
        vector<variable<double>> scope;
        ublas::matrix<double> K; // inverse variance if positive definite
        ublas::vector<double> h; // inverse variance left multiplied on mean vector (Kh = mu)
        double g;

        ublas::matrix<double> Sigma;
        ublas::vector<double> mu;


        canonical_form() = default;
        canonical_form(const canonical_form &f){
            for(auto member: f.scope){
                scope.push_back(member);
            }
            K = f.K;
            h = f.h;
            g = f.g;
            Sigma = f.Sigma;
            mu = f.mu;
        }
        explicit canonical_form(double g) : g(g){
            K.resize(1,1);
            h.resize(1,1);
            K(0,0) = 0;
            h(0) = 0;
        }

        canonical_form(set<string> scope, ublas::matrix<double> Sigma, ublas::vector<double> mu, double g)
                : g(g), Sigma(Sigma), mu(mu), K(inverse(Sigma)), h(ublas::prod(K, mu)) {}

        void compute_variance();
        void compute_mean();
        void compute_canonical_parameters();
        virtual double density( const ublas::vector<double> &x);
        virtual double log_density( const ublas::vector<double> &x);
    };

    /**
     * @brief: Compute the variance as the inverse of K
     */

    void canonical_form::compute_variance(){
        Sigma = inverse(K);
    }

    /**
     * @brief: Compute the mean as the left product of sigma with h
     */

    void canonical_form::compute_mean() {
        mu = ublas::prod(Sigma, h);
    }

    /**
     * @brief: Compute K and h from Sigma and mu
     */

    void canonical_form::compute_canonical_parameters() {
        // compute K and h
        K = inverse(Sigma);
        h = ublas::prod(K, mu);

        // compute g
        long n = mu.size();
        double mu_h = ublas::inner_prod(mu, h);
        double det_Sigma = det(Sigma);
        g = -0.5*mu_h - log(pow(2*pi, n/2)*sqrt(det_Sigma));
    }

    /**
     * @brief: Compute the probability density at x
     * @param x : value of random variable
     * @return : probability density of x
     */

    double canonical_form::density( const ublas::vector<double> &x){
        if(K.size1() != Sigma.size1()){
            compute_canonical_parameters();
        }



        ublas::vector<double> d = x - mu;
        ublas::vector<double> Kd = ublas::prod(K, d);
        double dKd = ublas::inner_prod(d, Kd);
        double N = x.size();
        double denom = pow(2*pi, N/2.0)*sqrt(det(Sigma));
        double density = exp(-0.5*dKd)/denom;

        return density;
    }

    double canonical_form::log_density(const ublas::vector<double> &x) {
        ublas::vector<double> Kx = ublas::prod(K, x);
        double xKx = ublas::inner_prod(x, Kx);
        double hx = ublas::inner_prod(h, x);
        return -0.5*xKx + hx + g;
    }

    struct skew_canonical_form : public canonical_form{
        ublas::vector<double> skew;

        skew_canonical_form() = default;
        skew_canonical_form(ublas::vector<double> skew) : skew(skew) {}

        skew_canonical_form(set<string> &scope, ublas::matrix<double> &Sigma,
                            ublas::vector<double> &mu, double g, ublas::vector<double> skew)
                : canonical_form(scope, Sigma, mu, g), skew(skew) {}


        double density( ublas::vector<double> &x);
    };

    /**
     * @brief Compute the density of a multivariate skew-normal distribution
     * @param x
     * @return
     */

    double skew_canonical_form::density( ublas::vector<double> &x) {
        // maintain a bit of notational consistency with MATLAB prototype
        ublas::vector<double> lambda = skew;
        // let N = dimension of the scope
        int N = scope.size();
        // compute delta
        ublas::matrix<double> Delta = eye(N);
        for(int i=0; i<N; ++i){
            Delta(i,i) = sqrt(1 - lambda[i]/sqrt(1 + lambda[i]*lambda[i]));
        }
        // compute alpha
        // numerator
        ublas::matrix<double> S_inv = inverse(Sigma);
        ublas::matrix<double> D_inv = inverse(Delta);
        ublas::matrix<double> SD_inv = ublas::prod(S_inv, D_inv);
        ublas::vector<double> numerator = ublas::prod(ublas::trans(lambda), SD_inv);

        // denominator
        ublas::vector<double> S_invL = ublas::prod(S_inv, lambda);
        double LS_invL = ublas::inner_prod(lambda, S_invL);
        double denominator = sqrt(1 + LS_invL);
        ublas::vector<double> alpha = numerator/denominator;

        // compute omega
        ublas::matrix<double> Lambda = ublas::outer_prod(lambda, lambda);
        ublas::matrix<double> Omega = Sigma + Lambda;
        Omega = ublas::prod(Delta, Omega);
        Omega = ublas::prod(Omega, Delta);

        // compute density
        ublas::matrix<double> O_inv = inverse(Omega);
        ublas::vector<double> r = x - mu;
        ublas::vector<double> O_invR = ublas::prod(O_inv, r);
        double power = ublas::inner_prod(r, O_invR);
        double det_O = det(Omega);
        denominator = sqrt(pow(2*pi, N) * det_O);
        double f = exp(-0.5 * power)/denominator;
        double erf_arg = ublas::inner_prod(alpha, r);
        f = f * (1 + erf(erf_arg));
        return f;
    }


    /**
     * \brief: Returns a vacuous form over the set of variables.
     * @param dimension
     * @return vacuous form
     *
     * Notes:
     * -passes manual test, automated test not necessary at this time
     */
    canonical_form vacuous_form(double dimension){

        canonical_form vac_form;

        vac_form.K.resize(dimension, dimension);
        vac_form.h.resize(dimension, 1);

        for(int i=0; i<dimension; ++i){
            vac_form.h(i) = 0;
            for(int j=0; j<dimension; ++j){
                vac_form.K(i,j) = 0;
            }
        }

        vac_form.g = 0;
        return vac_form;
    }


    /**
 * \brief: Computes factor product of two gaussian factors in canonical form representation.
 * @param form1 : first gaussian
 * @param form2 : second gaussian
 * @return : product of gaussians
 */

    canonical_form canonical_product( const canonical_form &form1, const canonical_form &form2 ){
        if(form1.scope.empty() and form2.scope.empty()){
            canonical_form prod (form1);
            prod.g += form2.g;
            return prod;
        }

        // initialize the scope of the product as the union of the scopes of the factors
        unordered_map<variable<double>*, int> map1, map2;
        set<variable<double>*> scope_check1, scope_check2;

        vector<variable<double>> scope;
        int j,i = 0;
        for( auto var : form1.scope){
            scope.push_back(var);
            scope_check1.insert(&scope.back());
            map1[&scope.back()] = i; ++i;
        }
        i = 0;
        for( auto var : form2.scope){
            if(scope_check1.find(&var) == scope_check1.end()){// this operation should be done using vector name strings
                scope.push_back(var);
            }
            scope_check2.insert(&scope.back());// keeps track of which scope elements come from which factors
            map2[&scope.back()] = i; ++i;
        }


        // compute product h and K matrix-this implementation is somewhat naive but faster wouldn't give dramatically better speed
        ublas::matrix<double> K1 (scope.size(), scope.size());
        ublas::matrix<double> K2 (scope.size(), scope.size());
        ublas::vector<double> h1 (scope.size(), 1);
        ublas::vector<double> h2 (scope.size(), 1);

        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for( auto v: scope){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(scope_check1.find(&v) == scope_check1.end()){
                h1(i) = 0;
            }
            else{
                h1(i) = form1.h(map1[&v]);
            }
            if(scope_check2.find(&v) == scope_check2.end()){
                h2(i) = 0;
            }
            else{
                h2(i) = form2.h(map2[&v]);
            }

            j = 0;
            // map K matrices
            for(auto w: scope){
                if(scope_check1.find(&v) == scope_check1.end() or scope_check1.find(&w) == scope_check1.end()){ K1(i,j) = 0; }
                else{
                    K1(i,j) = form1.K(map1[&v], map1[&w]);
                }
                if(scope_check2.find(&v) == scope_check2.end() or scope_check2.find(&w) == scope_check2.end()){ K2(i,j) = 0; }
                else{
                    K2(i,j) = form2.K(map2[&v], map2[&w]);
                }
                ++j;
            }
            ++i;
        }

        canonical_form product;
        product.K = K1 + K2;
        product.h = h1 + h2;
        product.g = form1.g + form2.g;

        return product;
    }

    /**
     * \brief : computes the quotient of two Gaussian canonical forms.
     * @param numerator : canonical form
     * @param denominator : canonical form
     * @return : quotient canonical form
     */
    canonical_form canonical_quotient( const canonical_form &numerator, const canonical_form &denominator){
        if(numerator.scope.empty() and denominator.scope.empty()){
            canonical_form quotient (numerator);
            quotient.g -= denominator.g;
            return quotient;
        }

        canonical_form form1 = numerator;
        canonical_form form2 = denominator;

        // initialize the scope of the product as the union of the scopes of the factors
        unordered_map<variable<double>*, int> map1, map2;
        set<variable<double>*> scope_check1, scope_check2;
        vector<variable<double>> scope;
        int j,i = 0;
        for( auto v : form1.scope){
            scope.push_back(v);
            scope_check1.insert(&v);
            map1[&scope.back()] = i; ++i;
        }
        i = 0;
        for( auto v : form2.scope){
            scope.push_back(v);
            scope_check2.insert(&v);
            map2[&scope.back()] = i; ++i;
        }


        // compute product h and K matrix-this implementation is somewhat naive but faster wouldn't give dramatically better speed
        ublas::matrix<double> K1 (scope.size(), scope.size());
        ublas::matrix<double> K2 (scope.size(), scope.size());
        ublas::vector<double> h1 (scope.size(), 1);
        ublas::vector<double> h2 (scope.size(), 1);

        set<string>::iterator it1, it2;
        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for(auto v: scope){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(scope_check1.find(&v) == scope_check1.end()){
                h1(i) = 0;
            }
            else{
                h1(i) = form1.h(map1[&v]);
            }
            if(scope_check2.find(&v) == scope_check2.end()){
                h2(i) = 0;
            }
            else{
                h2(i) = form2.h(map2[&v]);
            }

            j = 0;
            // map K matrices
            for(auto w: scope){
                if(scope_check1.find(&v) == scope_check1.end() or scope_check1.find(&w) == scope_check1.end()){ K1(i,j) = 0; }
                else{
                    K1(i,j) = form1.K(map1[&v], map1[&w]);
                }
                if(scope_check2.find(&v) == scope_check2.end() or scope_check2.find(&w) == scope_check2.end()){ K2(i,j) = 0; }
                else{
                    K2(i,j) = form2.K(map2[&v], map2[&w]);
                }
                ++j;
            }
            ++i;
        }

        canonical_form quotient;
        quotient.K = K1 - K2;
        quotient.h = h1 - h2;
        quotient.g = form1.g - form2.g;

        return quotient;
    }

    /**
     * \brief: Computes the marginal of a canonical form over the given variables.
     * @param joint_form : form prior to marginalization
     * @param vars : pointers to variables in the scope of joint form that will be marginalized out
     * @return : marginalizd form
     */
    canonical_form canonical_marginal( const canonical_form &joint_form, const vector<variable<double>*> &vars){
        if(vars.empty()){
            return joint_form;
        }
        set<variable<double>*> vars_to_elim;
        for(auto v: vars){
            vars_to_elim.insert(v);
        }

        // map variables to their order in the joint matrix
        unordered_map<variable<double>*, int> joint_map;
        int j,i = 0;
        for(auto v: joint_form.scope){
            joint_map[&v] = i;
            ++i;
        }
        // we want to refer to the data in the joint form throughout the construction of the marginal scope
        // therefore, we'll store a vector of pointers to the original scope for now, and use them
        // to construct a new canonical form later
        vector<variable<double>*> marginal_scope;
        for(auto v: joint_form.scope){
            if(vars_to_elim.find(&v) == vars_to_elim.end()){
                marginal_scope.push_back(&v);
            }
        }

        // populate sub-matrices and sub-vectors
        ublas::matrix<double> K_xx (marginal_scope.size(), marginal_scope.size());
        ublas::matrix<double> K_yy (vars.size(), vars.size());
        ublas::matrix<double> K_xy (marginal_scope.size(), vars.size());
        ublas::matrix<double> K_yx (vars.size(), marginal_scope.size());

        ublas::vector<double> h_x (marginal_scope.size());
        ublas::vector<double> h_y (vars.size());
        i = 0;
        for(auto pv: marginal_scope){
            h_x(i) = joint_form.h(joint_map[pv]);
            j = 0;
            for(auto pw: marginal_scope){
                K_xx(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: vars){
            h_y(i) = joint_form.h(joint_map[pv]);
            j = 0;
            for(auto pw: vars){
                K_yy(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: marginal_scope){
            j = 0;
            for(auto pw: marginal_scope){
                K_xy(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: marginal_scope){
            j = 0;
            for(auto pw: marginal_scope){
                K_yx(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }

        // compute marginal K matrix- note prec_prod computes higher precision matrix product than prod
        // term 1
        ublas::matrix<double> K_yy_inv = inverse(K_yy, true);// compute matrix inverse
        ublas::matrix<double> K_xy_K_yy_inv = ublas::prod(K_xy, K_yy_inv);
        ublas::matrix<double> K = K_xx - ublas::prod( K_xy_K_yy_inv, K_yx);

        // compute marginal h vector
        ublas::matrix<double> K_xyK_yyInv = ublas::prod(K_xy, K_yy_inv);
        ublas::vector<double> h_modifier = ublas::prod(K_xyK_yyInv, h_y);


        ublas::vector<double> h = h_x - h_modifier;
        // compute marginal g value
        ublas::vector<double> K_yyInvh_y = ublas::prod(K_yy_inv, h_y);
        double hKh = ublas::inner_prod(h_y, K_yyInvh_y);
        double g = joint_form.g+0.5*(log(2*pi*det(K_yy_inv, true)) + hKh);

        canonical_form marginal;
        for( i=0; i< marginal_scope.size(); ++i){
            marginal.scope.push_back(*marginal_scope[i]);
        }
        marginal.K = K;
        marginal.h = h;
        marginal.g = g;
        return marginal;
    }

    /**
     * \brief : Reduces a canonical form for known data
     * @param joint_form : form to be reduced
     * @param vars : variables to pointers to variables to assign
     * @param y : vector of assignments
     * @return : reduced form
     */

    canonical_form canonical_reduction(const canonical_form &joint_form, const vector<variable<double>*> &vars, const ublas::vector<double> &y){

        set<variable<double>*> vars_to_assign;
        for(auto v: vars){
            vars_to_assign.insert(v);
        }

        // map variables to their order in the joint matrix
        unordered_map<variable<double>*, int> joint_map;
        int j,i = 0;
        for(auto v: joint_form.scope){
            joint_map[&v] = i;
            ++i;
        }
        // we want to refer to the data in the joint form throughout the construction of the marginal scope
        // therefore, we'll store a vector of pointers to the original scope for now, and use them
        // to construct a new canonical form later
        vector<variable<double>*> unassigned_scope;
        for(auto v: joint_form.scope){
            if(vars_to_assign.find(&v) == vars_to_assign.end()){
                unassigned_scope.push_back(&v);
            }
        }

        // populate sub-matrices and sub-vectors
        ublas::matrix<double> K_xx (unassigned_scope.size(), unassigned_scope.size());
        ublas::matrix<double> K_yy (vars.size(), vars.size());
        ublas::matrix<double> K_xy (unassigned_scope.size(), vars.size());
        ublas::matrix<double> K_yx (vars.size(), unassigned_scope.size());

        ublas::vector<double> h_x (unassigned_scope.size());
        ublas::vector<double> h_y (vars.size());
        i = 0;
        for(auto pv: unassigned_scope){
            h_x(i) = joint_form.h(joint_map[pv]);
            j = 0;
            for(auto pw: unassigned_scope){
                K_xx(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: vars){
            h_y(i) = joint_form.h(joint_map[pv]);
            j = 0;
            for(auto pw: vars){
                K_yy(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: unassigned_scope){
            j = 0;
            for(auto pw: unassigned_scope){
                K_xy(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(auto pv: unassigned_scope){
            j = 0;
            for(auto pw: unassigned_scope){
                K_yx(i,j) = joint_form.K(joint_map[pv], joint_map[pw]);
                ++j;
            }
            ++i;
        }
        // assign new h vector
        ublas::vector<double> K_xy_y = ublas::prod(K_xy, y);
        ublas::vector<double> h = h_x - K_xy_y;

        // assign new g value
        ublas::vector<double> K_yyY = ublas::prod(y,K_yy);
        double yKy = ublas::inner_prod( K_yyY, y);

        double hy    = ublas::inner_prod(h_y, y);
        double g = joint_form.g + hy - 0.5*yKy;

        canonical_form reduced_form;
        for(auto pv: unassigned_scope){
            reduced_form.scope.push_back(*pv);
        }
        reduced_form.K = K_xx;
        reduced_form.h = h;
        reduced_form.g = g;

        return reduced_form;
    }

    /**
     * @brief : Collapses a mixture of two Gaussian canonical form into a single Gaussian canonical form
     * @param form_A : first form
     * @param form_B : second form
     * @return : collapsed form
     *
     * @note :
     * Necessarily an approximation of the Gaussian mixture but needed for most inference algorithms. Inaccurate
     * for Gaussians with widely separated means.
     * @note :
     * If continuous scopes are empy, converts gs from log to linear space, adds log(g)s and exponentiates.
     * Equivalent to exp(log(g_A) + log(g_B))
     */
    canonical_form canonical_collapse(const canonical_form &form_A, const canonical_form &form_B){
        if(form_A.scope.empty() and form_B.scope.empty()){
            double p_A = exp(form_A.g);
            double p_B = exp(form_B.g);
            double p_sum = p_A+ p_B;
            canonical_form s;
            s.g = log(p_sum);
            return s;
        }

        // implementation assumes scopes have same members and orders
        double det_K_A = det(form_A.K);
        double det_K_B = det(form_B.K);

        canonical_form collapsed_form;
        collapsed_form.scope = form_A.scope;
        collapsed_form.K.resize(form_A.K.size1(), form_B.K.size2());
        collapsed_form.h.resize(form_A.h.size(), 1);

        if( abs(det_K_A) < 1e-14 or abs(det_K_B) < 1e-14){
            double d = exp(form_A.g) + exp(form_B.g);
            collapsed_form.g = log(d);
        }
        else {// this functionality hasn't yet been tested-not necessary for TPD integration
            ublas::matrix<double> Sigma_A = inverse(form_A.K, true);
            ublas::matrix<double> Sigma_B = inverse(form_B.K, true);
            ublas::vector<double> mu_A = ublas::prod(Sigma_A, form_A.h);
            ublas::vector<double> mu_B = ublas::prod(Sigma_B, form_B.h);

            double n = form_A.scope.size();
            double det_Sigma_A = det(Sigma_A, true);
            double det_Sigma_B = det(Sigma_B, true);

            double temp_A = ublas::inner_prod(mu_A, form_A.h);
            double temp_B = ublas::inner_prod(mu_B, form_B.h);

            double d_A = form_A.g - 0.5*temp_A - log(pow(2*pi,n/2)*pow(det_Sigma_A, 0.5));
            double d_B = form_B.g - 0.5*temp_B - log(pow(2*pi,n/2)*pow(det_Sigma_B, 0.5));

            double w_A = log(d_A);
            double w_B = log(d_B);

            ublas::vector<double> mu_AB = w_A*mu_A + w_B*mu_B; // weighted sum of means
            ublas::vector<double> mu_diff_A = mu_A - mu_AB;
            ublas::vector<double> mu_diff_B = mu_B - mu_AB;
            ublas::matrix<double> mu_mat_A = ublas::outer_prod(mu_diff_A, mu_diff_A);
            ublas::matrix<double> mu_mat_B = ublas::outer_prod(mu_diff_B, mu_diff_B);
            ublas::matrix<double> Sigma = w_A*Sigma_A + w_B*Sigma_B + w_A*mu_mat_A + w_B*mu_mat_B;

            collapsed_form.K = inverse(Sigma);
            collapsed_form.h = ublas::prod(collapsed_form.K, mu_AB);
            // compute new g
            ublas::matrix<double> S_inv = inverse(Sigma);
            ublas::vector<double> S_inv_mu = ublas::prod(S_inv, mu_AB);
            double mu_S_inv_mu = ublas::inner_prod(mu_AB, S_inv_mu);
            double det_Sigma = det(Sigma, true);
            collapsed_form.g = -0.5*mu_S_inv_mu - log(pow(2*pi, n/2)*pow(det_Sigma, 0.5));
        }


        return collapsed_form;
    }

}
#endif //TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H
