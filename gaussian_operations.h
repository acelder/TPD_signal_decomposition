//
// Class definition for canonical Gaussian forms and associated arithmetic operations.
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

const double pi =  3.1415926535897;

/**
 * Gaussian distribution expressed as a canonical form, also stores mean and variance.
 */

namespace gml {
    // declarations
    struct canonical_form;

    canonical_form* vacuous_form(set <string> continuous_vars);
    canonical_form canonical_product(canonical_form *form1, canonical_form *form2);
    canonical_form canonical_quotient(const canonical_form *numerator, const canonical_form *denominator);
    canonical_form canonical_marginal(canonical_form *joint_form, set <string> vars);
    canonical_form canonical_reduction(canonical_form *form, set <string> vars, ublas::matrix<double> y);
    canonical_form canonical_collapse(canonical_form *formA, canonical_form *formB);

    // implementations

    struct canonical_form{
        set<string> scope;
        ublas::matrix<double> K; // inverse variance if positive definite
        ublas::vector<double> h; // inverse variance left multiplied on mean vector (Kh = mu)
        double g;

        ublas::matrix<double> Sigma;
        ublas::vector<double> mu;


        canonical_form() = default;
        explicit canonical_form(canonical_form* ptr){
            for(auto member: ptr->scope){
                scope.insert(member);
            }
            K = ptr->K;
            h = ptr->h;
            g = ptr->g;
            Sigma = ptr->Sigma;
            mu = ptr->mu;
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
        virtual double density( ublas::vector<double> &x);
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

    double canonical_form::density( ublas::vector<double> &x){
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
    canonical_form* vacuous_form(set<string> continuous_vars){
        double dimension = continuous_vars.size();
        canonical_form* vac_form = new canonical_form;
        vac_form->scope = continuous_vars;
        vac_form->K.resize(dimension, dimension);
        vac_form->h.resize(dimension, 1);

        for(int i=0; i<dimension; ++i){
            vac_form->h(i) = 0;
            for(int j=0; j<dimension; ++j){
                vac_form->K(i,j) = 0;
            }
        }

        vac_form->g = 0;
        return vac_form;
    }


    /**
 * \brief: Computes factor product of two gaussian factors in canonical form representation.
 * @param form1 : first gaussian
 * @param form2 : second gaussian
 * @return : product of gaussians
 */

    canonical_form canonical_product( canonical_form* form1, canonical_form* form2 ){
        // initialize the scope of the product as the union of the scopes of the factors
        unordered_map<string, int> map1, map2;
        set<string> scope;
        int j,i = 0;
        for( auto var : form1->scope){
            scope.insert(var);
            map1[var] = i; ++i;
        }
        i = 0;
        for( auto var : form2->scope){
            scope.insert(var);
            map2[var] = i; ++i;
        }


        // compute product h and K matrix-this implementation is somewhat naive but faster wouldn't give dramatically better speed
        ublas::matrix<double> K1 (scope.size(), scope.size());
        ublas::matrix<double> K2 (scope.size(), scope.size());
        ublas::vector<double> h1 (scope.size(), 1);
        ublas::vector<double> h2 (scope.size(), 1);

        set<string>::iterator it1, it2;
        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for(it1=scope.begin(); it1!=scope.end(); ++it1){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(form1->scope.find(*it1) == form1->scope.end()){
                h1(i) = 0;
            }
            else{
                h1(i) = form1->h(map1[*it1]);
            }
            if(form2->scope.find(*it1) == form2->scope.end()){
                h2(i) = 0;
            }
            else{
                h2(i) = form2->h(map2[*it1]);
            }

            j = 0;
            // map K matrices
            for(it2=scope.begin(); it2!=scope.end(); ++it2){
                if(form1->scope.find(*it1) == form1->scope.end() or form1->scope.find(*it2) == form1->scope.end()){ K1(i,j) = 0; }
                else{
                    K1(i,j) = form1->K(map1[*it1], map1[*it2]);
                }
                if(form2->scope.find(*it1) == form2->scope.end() or form2->scope.find(*it2) == form2->scope.end()){ K2(i,j) = 0; }
                else{
                    K2(i,j) = form2->K(map2[*it1], map2[*it2]);
                }
                ++j;
            }
            ++i;
        }

        canonical_form product;
        product.K = K1 + K2;
        product.h = h1 + h2;
        product.g = form1->g + form2->g;

        return product;
    }

    /**
     * \brief : computes the quotient of two Gaussian canonical forms.
     * @param numerator : canonical form
     * @param denominator : canonical form
     * @return : quotient canonical form
     */
    canonical_form canonical_quotient( const canonical_form* numerator, const canonical_form* denominator ){
        const canonical_form* form1 = numerator;
        const canonical_form* form2 = denominator;

        // initialize the scope of the product as the union of the scopes of the factors
        unordered_map<string, int> map1, map2;
        set<string> scope;
        int j,i = 0;
        for( auto var : form1->scope){
            scope.insert(var);
            map1[var] = i; ++i;
        }
        i = 0;
        for( auto var : form2->scope){
            scope.insert(var);
            map2[var] = i; ++i;
        }


        // compute product h and K matrix-this implementation is somewhat naive but faster wouldn't give dramatically better speed
        ublas::matrix<double> K1 (scope.size(), scope.size());
        ublas::matrix<double> K2 (scope.size(), scope.size());
        ublas::vector<double> h1 (scope.size(), 1);
        ublas::vector<double> h2 (scope.size(), 1);

        set<string>::iterator it1, it2;
        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for(it1=scope.begin(); it1!=scope.end(); ++it1){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(form1->scope.find(*it1) == form1->scope.end()){
                h1(i) = 0;
            }
            else{
                h1(i) = form1->h(map1[*it1]);
            }
            if(form2->scope.find(*it1) == form2->scope.end()){
                h2(i) = 0;
            }
            else{
                h2(i) = form2->h(map2[*it1]);
            }

            j = 0;
            // map K matrices
            for(it2=scope.begin(); it2!=scope.end(); ++it2){
                if(form1->scope.find(*it1) == form1->scope.end() or form1->scope.find(*it2) == form1->scope.end()){ K1(i,j) = 0; }
                else{
                    K1(i,j) = form1->K(map1[*it1], map1[*it2]);
                }
                if(form2->scope.find(*it1) == form2->scope.end() or form2->scope.find(*it2) == form2->scope.end()){ K2(i,j) = 0; }
                else{
                    K2(i,j) = form2->K(map2[*it1], map2[*it2]);
                }
                ++j;
            }
            ++i;
        }

        canonical_form quotient;
        quotient.K = K1 - K2;
        quotient.h = h1 - h2;
        quotient.g = form1->g - form2->g;

        return quotient;
    }

    /**
     * \brief: Computes the marginal of a canonical form over the given variables.
     * @param joint_form : form prior to marginalization
     * @param vars : set of variables to be marginalized
     * @return : marginalizd form
     */
    canonical_form canonical_marginal( canonical_form* joint_form, set<string> vars){
        // declare a random access iterator over sets of strings
        set<string>::iterator it1, it2;
        // map variables to their order in the joint matrix
        unordered_map<string, int> joint_map;
        int j,i = 0;
        for(it1=joint_form->scope.begin(); it1!=joint_form->scope.end(); ++it1){
            joint_map[(*it1)] = i;
            ++i;
        }
        // populate the marginal scopes
        set<string> marginal_scope;
        for(it1=joint_form->scope.begin(); it1!=joint_form->scope.end(); ++it1){
            if(vars.find(*it1) == vars.end()){
                marginal_scope.insert(*it1);
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
        for(it1=marginal_scope.begin(); it1!=marginal_scope.end(); ++it1){
            h_x(i) = joint_form->h(joint_map[*it1]);
            j = 0;
            for(it2=marginal_scope.begin(); it2!=marginal_scope.end(); ++it2){
                K_xx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            h_y(i) = joint_form->h(joint_map[*it1]);
            j = 0;
            for(it2=vars.begin(); it2!=vars.end(); ++it2){
                K_yy(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=marginal_scope.begin(); it1!=marginal_scope.end(); ++it1){
            j = 0;
            for(it2=vars.begin(); it2!=vars.end(); ++it2){
                K_xy(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            j = 0;
            for(it2=marginal_scope.begin(); it2!=marginal_scope.end(); ++it2){
                K_yx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
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
        double g = joint_form->g+0.5*(log(2*pi*det(K_yy_inv, true)) + hKh);

        canonical_form marginal;
        marginal.scope = marginal_scope;
        marginal.K = K;
        marginal.h = h;
        marginal.g = g;
        return marginal;
    }

    /**
     * \brief : Reduces a canonical form for known data
     * @param joint_form : form to be reduced
     * @param vars : variables to assign
     * @param y : vector of assignments
     * @return : reduced form
     */

    canonical_form canonical_reduction(canonical_form* joint_form, set<string> vars, ublas::vector<double> y){
        // declare a random access iterator over sets of strings
        set<string>::iterator it1, it2;
        // map variables to their order in the joint matrix
        unordered_map<string, int> joint_map;
        int j,i = 0;
        for(it1=joint_form->scope.begin(); it1!=joint_form->scope.end(); ++it1){
            joint_map[(*it1)] = i;
            ++i;
        }
        // populate the marginal scopes
        set<string> unassigned_scope;
        for(it1=joint_form->scope.begin(); it1!=joint_form->scope.end(); ++it1){
            if(vars.find(*it1) == vars.end()){
                unassigned_scope.insert(*it1);
            }
        }

        // populate sub-matrices and sub-vectors
        ublas::matrix<double> K_xx (unassigned_scope.size(), unassigned_scope.size());
        ublas::matrix<double> K_yy (vars.size(), vars.size());
        ublas::matrix<double> K_xy (unassigned_scope.size(), vars.size());
        ublas::matrix<double> K_yx (vars.size(), unassigned_scope.size());

        ublas::vector<double> h_x (unassigned_scope.size(), 1);
        ublas::vector<double> h_y (vars.size(), 1);
        i = 0;
        for(it1=unassigned_scope.begin(); it1!=unassigned_scope.end(); ++it1){
            h_x(i) = joint_form->h(joint_map[*it1]);
            j = 0;
            for(it2=unassigned_scope.begin(); it2!=unassigned_scope.end(); ++it2){
                K_xx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            h_y(i) = joint_form->h(joint_map[*it1]);
            j = 0;
            for(it2=vars.begin(); it2!=vars.end(); ++it2){
                K_yy(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=unassigned_scope.begin(); it1!=unassigned_scope.end(); ++it1){
            j = 0;
            for(it2=vars.begin(); it2!=vars.end(); ++it2){
                K_xy(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            j = 0;
            for(it2=unassigned_scope.begin(); it2!=unassigned_scope.end(); ++it2){
                K_yx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
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
        double g = joint_form->g + hy - 0.5*yKy;

        canonical_form reduced_form;
        reduced_form.scope = unassigned_scope;
        reduced_form.K = K_xx;
        reduced_form.h = h;
        reduced_form.g = g;

        return reduced_form;
    }

    /**
     * \brief : Collapses a mixture of two Gaussian canonical form into a single Gaussian canonical form
     * @param form_A : first form
     * @param form_B : second form
     * @return : collapsed form
     *
     * \note :
     * Necessarily an approximation of the Gaussian mixture but needed for most inference algorithms. Inaccurate
     * for Gaussians with widely separated means.
     */
    canonical_form canonical_collapse(canonical_form* form_A, canonical_form* form_B){
        // implementation assumes scopes have same members and orders
        double det_K_A = det(form_A->K);
        double det_K_B = det(form_B->K);

        canonical_form collapsed_form;
        collapsed_form.scope = form_A->scope;
        collapsed_form.K.resize(form_A->K.size1(), form_B->K.size2());
        collapsed_form.h.resize(form_A->h.size(), 1);

        if( abs(det_K_A) < 1e-14 or abs(det_K_B) < 1e-14){
            double d = exp(form_A->g) + exp(form_B->g);
            collapsed_form.g = log(d);
        }
        else {// this functionality hasn't yet been tested-not necessary for TPD integration
            ublas::matrix<double> Sigma_A = inverse(form_A->K, true);
            ublas::matrix<double> Sigma_B = inverse(form_B->K, true);
            ublas::vector<double> mu_A = ublas::prod(Sigma_A, form_A->h);
            ublas::vector<double> mu_B = ublas::prod(Sigma_B, form_B->h);

            double n = form_A->scope.size();
            double det_Sigma_A = det(Sigma_A, true);
            double det_Sigma_B = det(Sigma_B, true);

            double temp_A = ublas::inner_prod(mu_A, form_A->h);
            double temp_B = ublas::inner_prod(mu_B, form_B->h);

            double d_A = form_A->g - 0.5*temp_A - log(pow(2*pi,n/2)*pow(det_Sigma_A, 0.5));
            double d_B = form_B->g - 0.5*temp_B - log(pow(2*pi,n/2)*pow(det_Sigma_B, 0.5));

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
