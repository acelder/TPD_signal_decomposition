//
// Class definition and basic manipulation functions for factorization of joint hybrid distributions.
// Created by Alexander Elder on 1/26/18.
//
// 

#ifndef TPD_SIGNAL_DECOMPOSITION_FACTOR_OPERATIONS_H
#define TPD_SIGNAL_DECOMPOSITION_FACTOR_OPERATIONS_H

#include "gaussian_operations.h"
#include <string>
#include <set>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace lin_alg;
namespace ublas = boost::numeric::ublas;

namespace gml{
    /*
     * Declarations.
     */

    // classes
    class variable;
    class factor;

    // factor manipulation functions
    factor* factor_product( factor* factor1, factor* factor2, const string& name );
    factor* discrete_marginal( factor* target_factor, variable* target_variable, const string& name );


    struct variable{
        // variable name
        string name;
        // table CPD parameters
        int cardinality;
        // canonical form parameters
        ublas::vector<canonical_form*> canonical_column;
        // vector of factors with variable in scope
        vector<factor*> factors;
        // discrete indicator

        // constructor
        variable(string name = "default") : name(move(name)) {}
        variable(string name = "default", int cardinality = 1) : name(move(name)), cardinality(cardinality){
            canonical_column.resize(cardinality);
            for(int i=0; i<cardinality; ++i){
                canonical_column[i] = new canonical_form;
            }
        }

        virtual ~variable() {
            for(int i=0; i<canonical_column.size(); ++i){
                delete canonical_column[i];
            }
        }

    };


    /**
     * An object of class factor represents a conditional probability distribution (CPD) as a canonical table.
     * For more information, see Koller and Friedman p. 618.
     */
    class factor{
    public:
        // factor data
        string name;

        // scope and values
        set<variable*> scope;
        set<string> continuous_scope;// !! NOT YET IMPLEMENTED
        ublas::vector<canonical_form*> canonical_table;

        // Gaussian canonical form

        // vectors of parents and children pointers
        vector<factor*> parents;
        vector<factor*> children;

        // constructors and destructor
        explicit factor(string name = "default") : name(move(name)) {}
        explicit factor(factor *ptr);
        factor(string name, variable* var);

        virtual ~factor(){
            for(int i=0; i<canonical_table.size(); ++i){
                delete canonical_table[i];
            }
        }

        // factor manipulation functions
        void print();
        vector<int> compute_stride();
        ublas::matrix<int> assignment(const ublas::vector<int>& index);
        ublas::vector<int> index(const ublas::matrix<int>& assignment);
        void observe_evidence( variable* target_variable, int var_assign );
    };

    factor::factor(factor *ptr){
        name = ptr->name;
        scope = ptr->scope;
        continuous_scope = ptr->continuous_scope;
        canonical_table = ptr->canonical_table;
        parents = ptr->parents;
        children = ptr->children;
    }

    /**
     * Constructor for single variable factors
     * @param name : name of factor
     * @param var : variable in factor scope
     */
    factor::factor(string name, variable* var) : name(name) {
        scope.insert(var); // add variable to scope
        var->factors.push_back(this); // add reference to this factor to the variable's factor list
        // assign variable values to the factor
        canonical_table.resize(var->cardinality);

        // We want factors to have their own data, so that operations on factors can happen independent of underlying
        // variables. To this end, we'll copy each entry of canonical table of the input variable
        for(int i=0; i<var->cardinality; ++i){
            auto *form_copy = new canonical_form(var->canonical_column(i));
            canonical_table(i) = form_copy;
        }

    }

    /**
     * @brief Computes the strides for each variable in the scope.
     * @return : vector of strides ordered by variable order in scope set.
     */
    vector<int> factor::compute_stride(){
        int m = scope.size();

        // compute stride as a function of vector cardinalities
        vector<int> stride (m);
        stride[0] = 1;
        set<variable*>::iterator it = scope.begin();
        for(int i=1; i<m; ++i){
            stride[i] = stride[i-1]*(*it)->cardinality;
            ++it;
        }
        return stride;
    }

    /**
     * print tabulated assignment values and associated exp(g)
     */
    void factor::print(){
        cout << name << endl;
        for(auto it = scope.begin(); it != scope.end(); ++it){
            cout << (*it)->name << " ";
        }
        cout << endl;
        double num_entries = canonical_table.size();
        for(int i=0; i<num_entries; ++i){
            ublas::vector<int> index (1);
            index(0) = i;// assignment function requires vector type, overload later
            cout << assignment(index) << " " << exp(canonical_table(i)->g) << endl;
        }
        cout << endl;
    }

    /**
     * @brief compute matrix of discrete variable assignments
     * @param index : indices of each assignment
     * @return : matrix of numbers
     */
    ublas::matrix<int> factor::assignment(const ublas::vector<int>& index){

        size_t n = index.size();
        size_t m = scope.size();

        // extract joint cardinality from factor scope
        ublas::vector<int> joint_cardinality (m);
        auto it = scope.begin();
        for(int i=0; i<m; ++i){
            joint_cardinality(i) = (*it)->cardinality;
            ++it;
        }

        vector<int> stride = compute_stride();

        // compute assignments to the input indices
        ublas::matrix<int> assignments (n,m);
        for(int i=0; i<n; ++i){
            for(int j=0; j<m; ++j){
                assignments(i,j) = (index[i]/stride[j]) % joint_cardinality[j];
            }
        }
        return assignments;
    }

    /**
     * @brief Compute table index based on variable assignment
     * @param assignment : discrete variable assignment
     * @return : table index
     */
    ublas::vector<int> factor::index(const ublas::matrix<int>& assignment) {
        size_t n = assignment.size1(); // number of indices to return
        size_t m = assignment.size2(); // should be equal to number of parents + 1

        // extract joint cardinality from factor scope
        ublas::vector<int> joint_cardinality (m);
        int i = 0;
        for(auto& it : scope){
            joint_cardinality(i) = it->cardinality;
            ++i;
        }

        // compute stride as a function of joint cardinality
        vector<int> stride = compute_stride();

        ublas::vector<int> indices (n);
        for(int i=0; i<n; ++i){
            indices(i) = 0;
            for(int j=0; j<m; j++){
                indices(i) += assignment(i,j) * stride[j];
            }
        }
        return indices;
    }


    void factor::observe_evidence( variable* target_variable, int var_assign ){
        // naive implementation, not necessary to optimize as this operation is not limiting

        // get position of variable in set
        auto it = scope.find(target_variable);
        auto i = distance(scope.begin(), it);

        // compute assignment matrix
        size_t N = canonical_table.size();
        ublas::vector<int> ind (N);
        for(int j=0; j<N; ++j){ ind(j) = j; }// fill ind with indices

        ublas::matrix<int> A = assignment(ind);

        //Asserts floating point compatibility at compile time
        static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");

        for(int j=0; j<N; ++j){
            if(A(j,i) != var_assign){
                canonical_table[j]->g = -INFINITY;
            }
        }
    }

    /**
     * @brief Compute product of two factors
     * @param factor1 : pointer to left factor
     * @param factor2 : pointer to right factor
     * @param name : name of new factor
     * @return : product of factor1 and factor2
     */
    factor* factor_product( factor* factor1, factor* factor2, const string& name ){

        // check that the scopes of each factor are not empty
        if( factor1->scope.empty() ){ return factor2; }
        if( factor2->scope.empty() ){ return factor1; }

        set<variable*> scope;

        for(auto& it : factor1->scope){
            scope.insert(it);
        }
        for(auto& it : factor2->scope){
            scope.insert(it);
        }

        // compute strides
        vector<int> factor1_stride = factor1->compute_stride();
        vector<int> factor2_stride = factor2->compute_stride();

        unordered_map<string, int> stride1;
        unordered_map<string, int> stride2;

        // map variable names to their stride values in each factor
        int i = 0;
        for(auto& it : factor1->scope){
            stride1[it->name] = factor1_stride[i];
            ++i;
        }

        i = 0;
        for(auto& it : factor2->scope){
            stride2[it->name] = factor2_stride[i];
            ++i;
        }

        // iterate over the scope, computing the total number of product values as the product of the cardinalities
        // also, map variable names to stride=zero if they are not in a given factor's scope
        int n_vals = 1;
        for(auto& it : scope){
            n_vals = n_vals * it->cardinality;

            if( factor1->scope.find(it) == factor1->scope.end()){
                stride1[it->name] = 0;
            }
            if( factor2->scope.find(it) == factor2->scope.end()){
                stride2[it->name] = 0;
            }
        }

        // main loop
        ublas::vector<canonical_form*> canonical_table (n_vals);// this needs to be updated for table of pointers
        vector<int> assignment(scope.size(), 0); // initialize assignment vector

        ublas::vector<canonical_form*> phi1 = factor1->canonical_table; // extract value arrays
        ublas::vector<canonical_form*> phi2 = factor2->canonical_table;


        int j = 0, k = 0;
        for(i=0; i<canonical_table.size(); ++i){
            // canonical form product
            canonical_form p = canonical_product(phi1[j], phi2[k]);// compute the canonical product, which is locally output to the stack
            canonical_table[i] = new canonical_form(&p);// copy the canonical form to the heap to access it out of scope
            int l = 0;
            for(auto& it : scope){
                assignment[l]++;
                if(assignment[l] == it->cardinality){
                    assignment[l] = 0;
                    j = j - (it->cardinality - 1) * stride1[it->name];
                    k = k - (it->cardinality - 1) * stride2[it->name];
                }
                else{
                    j = j + stride1[it->name];
                    k = k + stride2[it->name];
                    break;
                }
                ++l;
            }
        }

        // construct product factor
        factor* product = new factor(name);
        product->scope = scope;
        product->canonical_table = canonical_table;// consider returning a pointer

        return product;
    }

    /**
     * \brief: Computes the marginalization of a canonical tabl CPD over a single variable.
     * @param joint_factor : pointer to factor to be marginalized
     * @param target_variable : pointer to variable to be marginalized
     * @param name : name of marginal factor
     * @param normalize : indicates whether to perform a normalization after marginalization
     * @return  : a (possibly normalized) marginalization of the join_factor
     */

    factor* discrete_marginal( factor* joint_factor, variable* target_variable, const string& name ){

        vector<int> stride = joint_factor->compute_stride();
        unordered_map<string, int> stride_joint;
        unordered_map<string, int> stride_marginal;

        set<variable*> marginal_scope;

        // construct stride maps for both the joint and marginal distributions
        int target_stride = 0, i = 0;
        for(auto& it : joint_factor->scope){
            stride_joint[it->name] = stride[i];
            if( it != target_variable ){
                marginal_scope.insert(it);
            }
            else{
                target_stride = stride[i];
            }
            ++i;
        }

        factor* marginal_factor = new factor(name); // construct new factor
        marginal_factor->scope = marginal_scope;
        stride = marginal_factor->compute_stride();
        i = 0;
        for(auto& it : marginal_factor->scope){
            stride_marginal[it->name] = stride[i];
        }

        // main loop
        long n_vals = joint_factor->canonical_table.size()/target_variable->cardinality; // number of values in marginalized table
        marginal_factor->canonical_table.resize(n_vals);
        vector<int> assignment (marginal_scope.size(), 0);
        int k = 0;
        for(i=0; i<n_vals; ++i){
            // compute marginal over current assignment
            canonical_form marginal_val = *joint_factor->canonical_table(k); // first value
            for(int j=1; j<target_variable->cardinality; ++j){
                int s = stride_joint[target_variable->name];
                canonical_form next_val = *joint_factor->canonical_table(k + j*s);// subsequent values
                marginal_val = canonical_collapse(&marginal_val, &next_val);
            }
            marginal_factor->canonical_table(i) = new canonical_form(&marginal_val);// move the locally computed marginal value onto the heap

            // update assignment vector
            int l = 0;
            for(auto& it : marginal_scope){
                ++assignment[l];
                if( assignment[l] == it->cardinality ){
                    assignment[l] = 0;
                }
                else{
                    break;
                }
                ++l;
            }
            // update k as the stride-weighted sum of assignments
            k = 0, l = 0;
            for(auto& it : marginal_scope){
                k += assignment[l]*stride_joint[it->name];
            }
        }


        return marginal_factor;
    }


}


#endif //TPD_SIGNAL_DECOMPOSITION_FACTOR_OPERATIONS_H
