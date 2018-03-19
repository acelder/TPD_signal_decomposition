//
// Created by Alexander Elder on 1/26/18.
//
//////Needs to be modified-functions should return pointers, not objects to be consistent with current type implementation.

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

    // factor declaration
    class factor;

    // factor manipulation functions
    factor factor_product( const factor &factor1, const factor &factor2, const string &name );
    factor continuous_marginal( const factor &joint_factor, variable<double>* target_variable, const string &name);
    factor discrete_marginal( const factor &joint_factor, variable<int>* target_variable, const string &name );
    factor unit_factor();


    /**
     * An object of class factor represents a conditional probability distribution (CPD) as a canonical table.
     * For more information, see Koller and Friedman p. 618.
     */
    class factor{
    public:
        // factor data
        string name;

        // scope and values
        vector<variable<int>> discrete_scope;
        vector<variable<double>> continuous_scope;
        ublas::vector<canonical_form> canonical_table;

        // Gaussian canonical form

        // vectors of parents and children pointers
        vector<factor*> parents;
        vector<factor*> children;

        // constructors and destructor
        factor() = default;
        explicit factor(string name = "default") : name(move(name)) {}
        explicit factor(factor *ptr);
        factor(const factor& phi);
        explicit factor(const variable<int> &A);
        explicit factor(const variable<double> &X);
        explicit factor(const vector<variable<int>> &Phi);
//        factor(const string &name, variable<int>* var);

        ~factor()= default;

        // factor manipulation functions
        void print();
        vector<int> compute_stride() const;
        ublas::matrix<int> assignment(const ublas::vector<int>& index);
        ublas::vector<int> index(const ublas::matrix<int>& assignment);
        void observe_evidence( const variable<int> &target_variable, int var_assign );
    };

    factor::factor(factor *ptr){
        name = ptr->name;
        discrete_scope = ptr->discrete_scope;
        continuous_scope = ptr->continuous_scope;
        canonical_table = ptr->canonical_table;
        parents = ptr->parents;
        children = ptr->children;
    }

    factor::factor(const factor &phi){
        name = phi.name;
        discrete_scope = phi.discrete_scope;
        continuous_scope = phi.continuous_scope;
        canonical_table = phi.canonical_table;
        parents = phi.parents;
        children = phi.children;
    }

    factor::factor(const variable<int> &A) : name(A.name){
        discrete_scope.push_back(A); // add variable to scope
        canonical_table.resize(A.cardinality);
    }

    factor::factor(const variable<double> &X) :name(X.name){
        continuous_scope.push_back(X); // add variable to continuous scope
        canonical_table.resize(X.cardinality);
    }

    factor::factor(const vector<variable<int>> &Phi){
        name = Phi.front().name;
        for(auto &phi: Phi){
            discrete_scope.push_back(phi);
        }
    }


    /**
     * @brief Computes the strides for each variable in the scope.
     * @return : vector of strides ordered by variable order in scope set.
     */
    vector<int> factor::compute_stride() const{
        size_t m = discrete_scope.size();

        // compute stride as a function of vector cardinalities
        vector<int> stride (m);
        stride[0] = 1;
        auto it = discrete_scope.begin();
        for(int i=1; i<m; ++i){
            stride[i] = stride[i-1]*(*it).cardinality;
            ++it;
        }
        return stride;
    }

    /**
     * print tabulated assignment values and associated exp(g)
     */
    void factor::print(){
        cout << "Factor name: " << name << endl;
        for(auto var: discrete_scope){
            cout << var.name << " ";
        }
        cout << endl;
        double num_entries = canonical_table.size();
        for(int i=0; i<num_entries; ++i){
            ublas::vector<int> index (1);
            index(0) = i;// assignment function requires vector type, overload later
            cout << assignment(index) << " " << exp(canonical_table(i).g) << endl;
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
        size_t m = discrete_scope.size();

        // extract joint cardinality from factor scope
        ublas::vector<int> joint_cardinality (m);
        for(int i=0; i<m; ++i){
            joint_cardinality(i) = discrete_scope[i].cardinality;
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
        for(auto& it : discrete_scope){
            joint_cardinality(i) = it.cardinality;
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


//    void factor::observe_evidence( const variable<int> &target_variable, int var_assign ){
//        // naive implementation, not necessary to optimize as this operation is not limiting
//
//        // get position of variable in set
//
//        auto it = find(discrete_scope.begin(), discrete_scope.end(), target_variable);//discrete_scope.find(target_variable);
//        auto i = distance(discrete_scope.begin(), it);
//
//        // compute assignment matrix
//        size_t N = canonical_table.size();
//        ublas::vector<int> ind (N);
//        for(int j=0; j<N; ++j){ ind(j) = j; }// fill ind with indices
//
//        ublas::matrix<int> A = assignment(ind);
//
//        //Asserts floating point compatibility at compile time
//        static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");
//
//        for(int j=0; j<N; ++j){
//            if(A(j,i) != var_assign){
//                canonical_table[j].g = -INFINITY;
//            }
//        }
//    }

    /**
     * @brief Compute product of two factors
     * @param factor1 : pointer to left factor
     * @param factor2 : pointer to right factor
     * @param name : name of new factor
     * @return : product of factor1 and factor2
     */
    factor factor_product( const factor &factor1, const factor &factor2, const string &name ){
        // initialize return factor
        size_t factor_1_scope_size = factor1.discrete_scope.size();
        size_t factor_2_scope_size = factor2.discrete_scope.size();

        // check that the scopes of each factor are not empty
        if( !factor_1_scope_size ){
            factor p (factor2);
            return p;
        }
        if( !factor_2_scope_size ){
            factor p (factor1);
            return p;
        }

        vector<variable<int>> product_scope;

        set<string> scope_union;

        for(auto &v: factor1.discrete_scope){
            scope_union.insert(v.name);
            product_scope.push_back(v);
        }
        for(auto &v: factor2.discrete_scope){
            if(scope_union.find(v.name) == scope_union.end()){
                scope_union.insert(v.name);
                product_scope.push_back(v);
            }
        }


        // compute strides
        vector<int> factor1_stride = factor1.compute_stride();
        vector<int> factor2_stride = factor2.compute_stride();

        unordered_map<string, int> stride1;
        unordered_map<string, int> stride2;

        // map discrete variable names to their stride values in each factor
        int i = 0;
        for(auto &v : factor1.discrete_scope){
            stride1[v.name] = factor1_stride[i];
            ++i;
        }

        i = 0;
        for(auto &v: factor2.discrete_scope){
            stride2[v.name] = factor2_stride[i];
            ++i;
        }

        // iterate over the scope, computing the total number of product values as the product of the cardinalities
        // also, map variable names to stride=zero if they are not in a given factor's scope
        int n_assignments = 1;
        for(auto& v : product_scope){
            n_assignments = n_assignments * v.cardinality;

            if( stride1.find(v.name) == stride1.end()){
                stride1[v.name] = 0;
            }
            if( stride2.find(v.name) == stride2.end()){
                stride2[v.name] = 0;
            }
        }

        // main loop
        ublas::vector<canonical_form> canonical_table (n_assignments);
        vector<int> assignment(n_assignments, 0); // initialize assignment vector-equals 1 if neither factors have discrete components

        ublas::vector<canonical_form> phi1 = factor1.canonical_table; // extract value arrays
        ublas::vector<canonical_form> phi2 = factor2.canonical_table;


        int j = 0, k = 0;
        for(i=0; i<canonical_table.size(); ++i){
            // canonical form product
            canonical_table[i] = canonical_product(phi1[j], phi2[k]);
            int l = 0;
            for(auto &v: product_scope){
                assignment[l]++;
                if(assignment[l] == v.cardinality){
                    assignment[l] = 0;
                    j = j - (v.cardinality - 1) * stride1[v.name];
                    k = k - (v.cardinality - 1) * stride2[v.name];
                }
                else{
                    j = j + stride1[v.name];
                    k = k + stride2[v.name];
                    break;
                }
                ++l;
            }
        }

        // construct product factor
        factor product(name);
        product.discrete_scope = product_scope;
        product.canonical_table = canonical_table;
        return product;
    }

    /**
     * @brief compute the marginal marginal over a continuous variable
     * @param joint_factor
     * @param target_variable
     * @param name
     * @return
     */

    factor continuous_marginal( const factor &joint_factor, variable<double>* target_variable, const string &name){
        // copy all but the factor being marginalized
        factor marginal(name);
        marginal.discrete_scope = joint_factor.discrete_scope;
        for(auto v: joint_factor.continuous_scope){
            if(&v != target_variable){
                marginal.continuous_scope.push_back(v);
            }
        }
        marginal.canonical_table.resize(joint_factor.canonical_table.size());
        for(int i=0; i<marginal.canonical_table.size(); ++i){
            vector<variable<double>*> target_variables;
            target_variables.push_back(target_variable);
            marginal.canonical_table[i] = canonical_marginal(joint_factor.canonical_table(i), target_variables);
        }

        return marginal;
    }

    /**
     * \brief: Computes the discrete marginalization of a canonical table CPD over a single variable.
     * @param joint_factor : pointer to factor to be marginalized
     * @param target_variable : pointer to variable to be marginalized
     * @param name : name of marginal factor
     * @param normalize : indicates whether to perform a normalization after marginalization
     * @return  : a (possibly normalized) marginalization of the join_factor
     */

    factor discrete_marginal(const factor &joint_factor, variable<int>* target_variable, const string &name  ){

        vector<int> stride = joint_factor.compute_stride();
        unordered_map<string, int> stride_joint;
        unordered_map<string, int> stride_marginal;
        string target_name = target_variable->name;

        vector<variable<int>> marginal_scope;

        // construct stride maps for the both joint and marginal distributions
        int target_stride = 0, i = 0;
        for(auto v: joint_factor.discrete_scope){

            if( v.name.compare(target_name) ){
                marginal_scope.push_back(v);
                stride_joint[v.name] = stride[i];
            }
            else{
                target_stride = stride[i];
            }
            ++i;
        }

        factor marginal_factor (name);
        marginal_factor.discrete_scope = marginal_scope;
        stride = marginal_factor.compute_stride();

        i = 0;
        for(auto& v: marginal_factor.discrete_scope){
            stride_marginal[v.name] = stride[i];
            ++i;
        }

        // main loop
        long n_vals = joint_factor.canonical_table.size()/target_variable->cardinality; // number of values in marginalized table
        marginal_factor.canonical_table.resize(n_vals);
        vector<int> assignment (marginal_scope.size(), 0);
        cout << endl;
        int k = 0;
        for(i=0; i<n_vals; ++i){
            // compute marginal over current assignment
            canonical_form marginal_val = joint_factor.canonical_table(k); // first value
            for(int j=1; j<target_variable->cardinality; ++j){// j might need to start from 1..
                int s = target_stride;
                canonical_form next_val = joint_factor.canonical_table(k + j*s);// subsequent values
                marginal_val = canonical_collapse(marginal_val, next_val);// this may not be correct in current implementation
            }
            marginal_factor.canonical_table(i) = marginal_val;// move the locally computed marginal value onto the heap

            // update assignment vector
            int l = 0;
            for(auto& v : marginal_scope){
                ++assignment[l];
                if( assignment[l] == v.cardinality ){
                    assignment[l] = 0;
                    k = k - (v.cardinality - 1) * stride_joint[v.name];
                }
                else{
                    k = k + stride_joint[v.name];
                    break;
                }
                ++l;
            }
        }


        return marginal_factor;
    }

    factor unit_factor(){
        factor unit("1");
        return unit;
    }

}


#endif //TPD_SIGNAL_DECOMPOSITION_FACTOR_OPERATIONS_H
