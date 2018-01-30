//
// Created by Alexander Elder on 1/26/18.
//

#ifndef TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H
#define TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H

#include "lin_alg.h"
#include <string>
#include <set>
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace lin_alg;
namespace ublas = boost::numeric::ublas;

const double pi =  3.1415926535897;

namespace gml {
    // declarations
    struct canonical_form;

    canonical_form vacuous_form(set <string> continuous_vars);
    canonical_form canonical_product(canonical_form *form1, canonical_form *form2);
    canonical_form canonical_quotient(const canonical_form *numerator, const canonical_form *denominator);
    canonical_form canonical_marginal(canonical_form *joint_form, set <string> vars);
    canonical_form canonical_reduction(canonical_form *form, set <string> vars, ublas::matrix<double> y);
    canonical_form canonical_collapse(canonical_form *formA, canonical_form *formB);

    // implementations

    struct canonical_form{
        set<string> scope;
        ublas::matrix<double> K; // inverse variance if positive definite
        ublas::matrix<double> h; // inverse variance left multiplied on mean vector (Kh = mu)
        double g;

        ublas::matrix<double> Sigma;
        ublas::matrix<double> mu;


        canonical_form() = default;

        explicit canonical_form(double g) : g(g){
            K.resize(1,1);
            h.resize(1,1);
            K(0,0) = 0;
            h(0,0) = 0;
        }

        void compute_variance();
        void compute_mean();
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
     * \brief: Returns a vacuous form over the set of variables.
     * @param dimension
     * @return vacuous form
     *
     * Notes:
     * -passes manual test, automated test not necessary at this time
     */
    canonical_form vacuous_form(set<string> continuous_vars){
        double dimension = continuous_vars.size();
        canonical_form vac_form;
        vac_form.scope = continuous_vars;
        vac_form.K.resize(dimension, dimension);
        vac_form.h.resize(dimension, 1);

        for(int i=0; i<dimension; ++i){
            vac_form.h(i,0) = 0;
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
        ublas::matrix<double> h1 (scope.size(), 1);
        ublas::matrix<double> h2 (scope.size(), 1);

        set<string>::iterator it1, it2;
        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for(it1=scope.begin(); it1!=scope.end(); ++it1){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(form1->scope.find(*it1) == form1->scope.end()){
                h1(i,0) = 0;
            }
            else{
                h1(i,0) = form1->h(map1[*it1],0);
            }
            if(form2->scope.find(*it1) == form2->scope.end()){
                h2(i,0) = 0;
            }
            else{
                h2(i,0) = form2->h(map2[*it1], 0);
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
        ublas::matrix<double> h1 (scope.size(), 1);
        ublas::matrix<double> h2 (scope.size(), 1);

        set<string>::iterator it1, it2;
        i = 0;
        // loops over both matrices can be coded compactly due to matrix symmetry
        for(it1=scope.begin(); it1!=scope.end(); ++it1){
            // if the variable pointed to by it1 is not in a given form scope, set its h value to zero, otherwise maps
            if(form1->scope.find(*it1) == form1->scope.end()){
                h1(i,0) = 0;
            }
            else{
                h1(i,0) = form1->h(map1[*it1],0);
            }
            if(form2->scope.find(*it1) == form2->scope.end()){
                h2(i,0) = 0;
            }
            else{
                h2(i,0) = form2->h(map2[*it1], 0);
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

        ublas::matrix<double> h_x (marginal_scope.size(),1);
        ublas::matrix<double> h_y (vars.size(),1);
        i = 0;
        for(it1=marginal_scope.begin(); it1!=marginal_scope.end(); ++it1){
            h_x(i,0) = joint_form->h(joint_map[*it1],0);
            j = 0;
            for(it2=marginal_scope.begin(); it2!=marginal_scope.end(); ++it2){
                K_xx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            h_y(i,0) = joint_form->h(joint_map[*it1], 0);
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

        ublas::matrix<double> Kh_inv = ublas::prod(K_yy_inv, h_y);
        ublas::matrix<double> K_xy_Kh_inv = ublas::prod(K_xy, K_yy_inv);


        // compute marginal h vector
        ublas::matrix<double> h = h_x - K_xy_Kh_inv;
        // compute marginal g value
        ublas::matrix<double> hKh = ublas::prod(ublas::trans(h_y), Kh_inv);
        double g = joint_form->g+0.5*(log(2*pi*det(K_yy_inv, true)) + hKh(0,0));

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

    canonical_form canonical_reduction(canonical_form* joint_form, set<string> vars, ublas::matrix<double> y){
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

        ublas::matrix<double> h_x (unassigned_scope.size(), 1);
        ublas::matrix<double> h_y (vars.size(), 1);
        i = 0;
        for(it1=unassigned_scope.begin(); it1!=unassigned_scope.end(); ++it1){
            h_x(i,0) = joint_form->h(joint_map[*it1], 0);
            j = 0;
            for(it2=unassigned_scope.begin(); it2!=unassigned_scope.end(); ++it2){
                K_xx(i,j) = joint_form->K(joint_map[*it1], joint_map[*it2]);
                ++j;
            }
            ++i;
        }
        i = 0;
        for(it1=vars.begin(); it1!=vars.end(); ++it1){
            h_y(i,0) = joint_form->h(joint_map[*it1], 0);
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
        ublas::matrix<double> h = h_x - prod(K_xy, y);

        // assign new g value
        ublas::matrix<double> K_yyY = prod(y,K_yy);
        ublas::matrix<double> yKy   = prod( K_yyY, y);

        ublas::matrix<double> hy    = prod(h_y, trans(y));
        double g = joint_form->g + hy(0,0) - 0.5*yKy(0,0);

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
        collapsed_form.h.resize(form_A->h.size1(), 1);

        if( abs(det_K_A) < 1e-14 or abs(det_K_B) < 1e-14){
            double d = exp(form_A->g) + exp(form_B->g);
            collapsed_form.g = log(d);
        }
        else {// this functionality hasn't yet been tested-not necessary for TPD integration
            ublas::matrix<double> Sigma_A = inverse(form_A->K, true);
            ublas::matrix<double> Sigma_B = inverse(form_B->K, true);
            ublas::matrix<double> mu_A = ublas::prod(Sigma_A, form_A->h);
            ublas::matrix<double> mu_B = ublas::prod(Sigma_B, form_B->h);

            double n = form_A->scope.size();
            double det_Sigma_A = det(Sigma_A, true);
            double det_Sigma_B = det(Sigma_B, true);

            ublas::matrix<double> temp_A = ublas::prod(ublas::trans(mu_A), form_A->h);
            ublas::matrix<double> temp_B = ublas::prod(ublas::trans(mu_B), form_B->h);

            double d_A = form_A->g - 0.5*temp_A(0,0) - log(pow(2*pi,n/2)*pow(det_Sigma_A, 0.5));
            double d_B = form_B->g - 0.5*temp_B(0,0) - log(pow(2*pi,n/2)*pow(det_Sigma_B, 0.5));

            double w_A = log(d_A);
            double w_B = log(d_B);

            ublas::matrix<double> mu = w_A*mu_A + w_B*mu_B;
            ublas::matrix<double> mu_diff_A = mu_A - mu;
            ublas::matrix<double> mu_diff_B = mu_B - mu;
            ublas::matrix<double> Sigma = w_A*Sigma_A + w_B*Sigma_B + w_A*ublas::prod(mu_diff_A, ublas::trans(mu_diff_A))
                                          + w_B*ublas::prod(mu_diff_B, ublas::trans(mu_diff_B));

            collapsed_form.K = inverse(Sigma);
            collapsed_form.h = ublas::prod(collapsed_form.K, mu);
            // compute new g
            ublas::matrix<double> S_inv = inverse(Sigma);
            ublas::matrix<double> S_inv_mu = ublas::prod(S_inv, mu);
            ublas::matrix<double> mu_S_inv_mu = ublas::prod(trans(mu), S_inv_mu);
            double det_Sigma = det(Sigma, true);
            collapsed_form.g = -0.5*mu_S_inv_mu(0,0) - log(pow(2*pi, n/2)*pow(det_Sigma, 0.5));
        }


        return collapsed_form;
    }

}
#endif //TPD_SIGNAL_DECOMPOSITION_GAUSSIAN_OPERATIONS_H
