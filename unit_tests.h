//
// unit_tests.cpp
// Created by Alexander Elder on 1/15/18.
//
// Contains unit tests for the TPD_signal_decomposition program. Tests are organized as follows:
// * classes
//   * constructor tests
//   * member function tests
// * functions
//   * test on appropriate data
//   * test on inappropriate data

#ifndef TPD_SIGNAL_DECOMPOSITION_UNIT_TESTS_H
#define TPD_SIGNAL_DECOMPOSITION_UNIT_TESTS_H

#include "factor_operations.h"
#include "gaussian_operations.h"
#include <random>
#include <string>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "storage_adaptors.h"

using namespace std;
using namespace gml;
using namespace lin_alg;

typedef boost::numeric::ublas::vector< double > double_vector;
typedef boost::numeric::ublas::matrix< double > double_matrix;
typedef boost::numeric::ublas::vector< int > int_vector;
typedef boost::numeric::ublas::matrix< int > int_matrix;

/*
 * Declare test functions:
 */
bool canonical_form_constructor_test();
bool variable_constructor_test();
bool factor_constructor_test1();
bool p_graph_constructor_test1();
bool compute_stride_test();
bool assignment_test();
bool index_test();
bool canonical_product_test();
bool canonical_quotient_test();
bool canonical_marginal_test();
bool canonical_reduction_test();
bool factor_product_test1();
bool factor_product_test2();
bool discrete_marginal_test1();
void canonical_form_tests();
void variable_tests();
void factor_tests();
void run_unit_tests();

// --------------------------------------- class tests -------------------------------------------

// ##### canonical_form #####
bool canonical_form_constructor_test(){
    // initialize random number generator
    uniform_real_distribution<double> unif(0,10);
    default_random_engine re;

    double_matrix K_test (3,3);
    double_matrix h_test (3,1);
    double g_test = unif(re);
    for(int i=0; i<3; ++i){
        h_test(i,0) = unif(re);
        for(int j=0; j<3; ++j){
            K_test(i,j) = unif(re);
        }
    }

    set<string> scope_test;
    scope_test.insert("A");
    scope_test.insert("B");
    scope_test.insert("C");

    canonical_form cf_test;
    cf_test.scope = scope_test;
    cf_test.K = K_test;
    cf_test.h = h_test;
    cf_test.g = g_test;

    bool pass_flag = true;
    if((*cf_test.scope.begin()) != (*cf_test.scope.begin())){
        pass_flag = false;
    }
    for(int i=0; i<3; ++i){
        if( abs(h_test(i,0) - cf_test.h(i,0)) > 1e-6 ){
            pass_flag = false;
        }
        for(int j=0; j<3; ++j){
            if( abs(K_test(i,j) - cf_test.K(i,j)) > 1e-6 ){
                pass_flag = false;
            }
        }
    }
    if( abs(g_test - g_test) > 1e-6 ){
        pass_flag = false;
    }

    return pass_flag;
}

// ##### variable #####
bool variable_constructor_test(){
    string test_name = "A";
    int test_card = 2;
    variable test_variable(test_name, test_card);

    bool pass_flag = true;
    if(test_name != test_variable.name){
        pass_flag = false;
    }
    if(test_card != test_variable.cardinality or test_card != test_variable.canonical_column.size()){
        pass_flag = false;
    }
    return pass_flag;
}

// ##### factor #####
bool factor_constructor_test1(){
    string test_name = "phi1";
    variable var1 ("A", 2);
    var1.canonical_column[0].g = log(0.4);
    var1.canonical_column[0].g = log(0.6);

    factor phi1 (test_name, &var1);

    bool pass_flag = true;

    auto i = phi1.scope.begin();
    if( (*i)->name != var1.name ){
        pass_flag = false;
    }
    for(int j=0; j<2; ++j){
        if( abs( phi1.canonical_table(j).g - var1.canonical_column[j].g ) > 1e-3 ){
            pass_flag = false;
        }
    }
    return pass_flag;

}

// ##### p_graph #####

//bool p_graph_constructor_test1(){
//    variable test_var("A", 2);
//    factor test_factor("phi1", &test_var);
//
//    p_graph test_graph;
//    test_graph.variables.push_back(&test_var);
//    test_graph.factors.push_back(&test_factor);
//
//    bool pass_flag = true;
//
//    if(test_graph.variables[0] != &test_var or test_graph.factors[0] != &test_factor){
//        pass_flag = false;
//    }
//
//    return pass_flag;
//}

// -------------------------------------- member function tests -------------------------------------------

// ###### factor member functions #######

bool compute_stride_test(){
    variable A ("A", 2);
    variable B ("B", 3);
    variable C ("C", 2);

    factor phi("phi1");
    phi.scope.insert(&A);
    phi.scope.insert(&B);
    phi.scope.insert(&C);

    vector<int> stride = phi.compute_stride();

    bool pass_flag = true;

    if(stride[0] != 1){
        pass_flag = false;
    }
    if(stride[1] != 2){
        pass_flag = false;
    }
    if(stride[2] != 6){
        pass_flag = false;
    }

    return pass_flag;
}

bool assignment_test(){
    variable A ("A", 2);
    variable B ("B", 3);
    variable C ("C", 2);

    factor phi("phi1");
    phi.scope.insert(&A);
    phi.scope.insert(&B);
    phi.scope.insert(&C);

    int_vector index (12);
    for(int j=0; j<12; ++j){ index[j] = j; }

    int_matrix assignment_test = phi.assignment(index);

    bool pass_flag = true;

    int assignment[12][3] = {
            0,0,0,
            1,0,0,
            0,1,0,
            1,1,0,
            0,2,0,
            1,2,0,
            0,0,1,
            1,0,1,
            0,1,1,
            1,1,1,
            0,2,1,
            1,2,1
    };


    for(int i=0; i<12; ++i){
        for(int j=0; j<3; ++j){
            if(assignment[i][j] != assignment_test(i,j)){
                pass_flag = false;
            }
        }
    }
    return pass_flag;
}

bool index_test(){
    variable A ("A", 2);
    variable B ("B", 3);
    variable C ("C", 2);

    factor phi("phi1");
    phi.scope.insert(&A);
    phi.scope.insert(&B);
    phi.scope.insert(&C);


    int value[12][3] = {
            0,0,0,
            1,0,0,
            0,1,0,
            1,1,0,
            0,2,0,
            1,2,0,
            0,0,1,
            1,0,1,
            0,1,1,
            1,1,1,
            0,2,1,
            1,2,1
    };

    int_matrix assignments (3,12);
    assignments = ublas::make_matrix_from_pointer(value);

    int_vector index_test = phi.index(assignments);

    bool pass_flag = true;

    for(int i=0; i<index_test.size(); ++i){
        if( index_test[i] != i){
            pass_flag = false;
            break;
        }
    }

    return pass_flag;
}

// ###### canonical form manipulation functions #######

bool canonical_product_test(){
    auto f1 = new canonical_form;

    f1->scope.insert("A");
    f1->scope.insert("B");
    f1->K = eye(2);

    f1->h.resize(2,1);
    f1->h(0,0) = 1.0;
    f1->h(1,0) = 2.0;

    f1->g = -3.0;

    auto f2 = new canonical_form;

    f2->scope.insert("B");
    f2->scope.insert("C");
    f2->K = eye(2);
    f2->K = 2.0*f2->K;

    f2->h.resize(2,1);
    f2->h(0,0) = -1.0;
    f2->h(1,0) = 1.0;

    f2->g = 1;


    canonical_form f_prod = canonical_product(f1, f2);

    double correct_K[3][3] = {
            1.0,0.0,0.0,
            0.0,3.0,0.0,
            0.0,0.0,2.0
    };

    double correct_h[3] = {1.0,1.0,1.0};

    double correct_g = -2.0;

    bool pass_flag = true;

    for(int i=0; i<3; ++i){
        if( abs(f_prod.h(i,0) - correct_h[i]) > 1e-8 ){
            pass_flag = false;
            break;
        }
        for(int j=0; j<3; ++j){
            if(abs(f_prod.K(i,j) - correct_K[i][j]) > 1e-8){
                pass_flag = false;
                break;
            }
        }
        if(!pass_flag){
            break;
        }
    }

    if( abs(f_prod.g - correct_g) > 1e-8 ){
        pass_flag = false;
    }

    delete f1;
    delete f2;
    return pass_flag;

}


bool canonical_quotient_test(){
    auto f1 = new canonical_form;

    f1->scope.insert("A");
    f1->scope.insert("B");
    f1->K = eye(2);

    f1->h.resize(2,1);
    f1->h(0,0) = 1.0;
    f1->h(1,0) = 2.0;

    f1->g = -3.0;

    auto f2 = new canonical_form;

    f2->scope.insert("B");
    f2->scope.insert("C");
    f2->K = eye(2);
    f2->K = 2.0*f2->K;

    f2->h.resize(2,1);
    f2->h(0,0) = -1.0;
    f2->h(1,0) = 1.0;

    f2->g = 1;

//    p_graph g;
    canonical_form f_quot = canonical_quotient(f1, f2);

    double correct_K[3][3] = {
            1.0, 0.0, 0.0,
            0.0,-1.0, 0.0,
            0.0, 0.0, -2.0
    };

    double correct_h[3] = {1.0,3.0,-1.0};

    double correct_g = -4.0;

    bool pass_flag = true;

    for(int i=0; i<3; ++i){
        if( abs(f_quot.h(i,0) - correct_h[i]) > 1e-8 ){
            pass_flag = false;
            break;
        }
        for(int j=0; j<3; ++j){
            if(abs(f_quot.K(i,j) - correct_K[i][j]) > 1e-8){
                pass_flag = false;
                break;
            }
        }
        if(!pass_flag){
            break;
        }
    }

    if( abs(f_quot.g - correct_g) > 1e-8 ){
        pass_flag = false;
    }

    delete f1;
    delete f2;
    return pass_flag;

}

bool canonical_marginal_test(){

    // construct test canonical form
    double test_K_values[3][3] = {
            4,1,2,
            1,2,1,
            2,1,3
    };

    double test_h_values[3][1] = {3,2,1};

    double test_g = 1;

    canonical_form test_form;
    test_form.scope.insert("A");
    test_form.scope.insert("B");
    test_form.scope.insert("C");
    test_form.K.resize(3,3);
    test_form.K = ublas::make_matrix_from_pointer(test_K_values);

    test_form.h.resize(3,1);
    test_form.h = ublas::make_matrix_from_pointer(test_h_values);

    test_form.g = 1;

    set<string> marginalized_vars;
    marginalized_vars.insert("C");

    // run canonical_marginal on test form
    canonical_form test_marginal = canonical_marginal(&test_form, marginalized_vars);

    // evaluate results
    bool pass_flag = true;
    double true_marginal_K[2][2] = {
            2.6667,    0.3333,
            0.3333,    1.6667
    };
    double true_marginal_h[2] = {2.3333, 1.6667};


    for(int i=0; i<2; ++i){
        if( abs(test_marginal.h(i,0) - true_marginal_h[i]) > 1e-3 ){ pass_flag = false; }
        for(int j=0; j<2; ++j){
            if( abs(test_marginal.K(i,j) - true_marginal_K[i][j]) > 1e-3 ){ pass_flag = false; break; }
        }
    }

    if( abs(test_marginal.g - 1.5363) > 1e-3){ pass_flag = false; }

    return pass_flag;
}

/**
 * \brief: test canonical_reduction function
 * @return : pass_flag
 */
bool canonical_reduction_test(){
    // construct test canonical form
    double test_K_values[3][3] = {
            4,1,2,
            1,2,1,
            2,2,3
    };

    double test_h_values[3][1] = {
            3,
            2,
            1
    };

    canonical_form test_form;
    test_form.scope.insert("A");
    test_form.scope.insert("B");
    test_form.scope.insert("C");
    test_form.K.resize(3,3);
    test_form.K = ublas::make_matrix_from_pointer(test_K_values);

    test_form.h.resize(3,1);
    test_form.h = ublas::make_matrix_from_pointer(test_h_values);

    test_form.g = 1;

    set<string> reduced_vars;
    reduced_vars.insert("C");

    double_matrix assignment(1,1);
    assignment(0,0) = 3;

    // run canonical_reduction on test_form
    canonical_form test_reduction = canonical_reduction(&test_form, reduced_vars, assignment);

    // evaluate test results
    double correct_reduced_K[2][2] = {
            4.0, 1.0,
            1.0, 2.0
    };
    double correct_reduced_h[2] = {-3.0, -1.0};

    bool pass_flag = true;

    for(int i=0; i<2; ++i){
        if( abs(test_reduction.h(i,0) - correct_reduced_h[i]) > 1e-3 ){ pass_flag = false; }
        for(int j=0; j<2; ++j){
            if( abs(test_reduction.K(i,j) - correct_reduced_K[i][j]) > 1e-3){ pass_flag = false; }
        }
    }

    if( abs(test_reduction.g - (-9.5)) > 1e-3 ){ pass_flag = false; }

    return pass_flag;
}

/**
 * \brief: test the product of two discrete factors, the first will be an "all one" factor
 * @return : bool value indicating test pass
 */
bool factor_product_test1(){
    // define continuous variable
    string c_var = "X";
    set<string> c_vars;
    c_vars.insert(c_var);

    // initialize identity
    variable D ("D", 2);
    D.canonical_column[0] = vacuous_form(c_vars);
    D.canonical_column[1] = vacuous_form(c_vars);

    // initialize non-identity factor
    variable C ("C", 2);
    C.canonical_column[0] = vacuous_form(c_vars);
    C.canonical_column[1] = vacuous_form(c_vars);

    C.canonical_column[0].g = log(0.1);
    C.canonical_column[1].g = log(0.9);

    factor phi_1 ("phi_1", &C);
    factor phi_2 ("phi_2", &D);

    factor phi_3 = factor_product(&phi_1, &phi_2, "phi_3");

    bool pass_flag = true;

    if( abs(phi_3.canonical_table[0].g - log(0.1)) > 1e-5 or
            abs(phi_3.canonical_table[1].g - log(0.9)) > 1e-5 or
            abs(phi_3.canonical_table[2].g - log(0.1)) > 1e-5 or
            abs(phi_3.canonical_table[3].g - log(0.9)) > 1e-5){
        pass_flag = false;
    }

    return pass_flag;
}

/**
 * \brief: second factor product test
 * @return: pass_flag
 */
bool factor_product_test2(){
    // define continuous variable
    string c_var = "X";
    set<string> c_vars;
    c_vars.insert(c_var);

    // initialize non-identity factor
    variable D ("D", 2);
    D.canonical_column[0] = vacuous_form(c_vars);
    D.canonical_column[1] = vacuous_form(c_vars);

    D.canonical_column[0].g = log(0.1);
    D.canonical_column[1].g = log(0.9);

    // initialize non-identity factor
    variable C ("C", 2);
    C.canonical_column[0] = vacuous_form(c_vars);
    C.canonical_column[1] = vacuous_form(c_vars);

    C.canonical_column[0].g = log(0.4);
    C.canonical_column[1].g = log(0.6);


    factor phi_1 ("phi_1", &C);
    factor phi_2 ("phi_2", &D);

    factor phi_3 = factor_product(&phi_1, &phi_2, "phi_3");

    bool pass_flag = true;

    if( abs(phi_3.canonical_table[0].g - log(0.04)) > 1e-5 or
        abs(phi_3.canonical_table[1].g - log(0.06)) > 1e-5 or
        abs(phi_3.canonical_table[2].g - log(0.36)) > 1e-5 or
        abs(phi_3.canonical_table[3].g - log(0.54)) > 1e-5){
        pass_flag = false;
    }

    return pass_flag;
}

/**
 * \brief: test discrete_marginal function
 * @return: pass_flag
 */

bool discrete_marginal_test1(){
    // define continuous variable
    string c_var = "X";
    set<string> c_vars;
    c_vars.insert(c_var);

    // initialize non-trivial factor
    variable C ("C", 2);
    C.canonical_column[0] = vacuous_form(c_vars);
    C.canonical_column[1] = vacuous_form(c_vars);
    C.canonical_column[0].g = log(0.4);
    C.canonical_column[1].g = log(0.6);

    // initialize all one factor
    variable D ("D", 2);
    D.canonical_column[0] = vacuous_form(c_vars);
    D.canonical_column[1] = vacuous_form(c_vars);

    D.canonical_column[0].g = log(0.1);
    D.canonical_column[1].g = log(0.9);

    factor phi_1 ("phi_1", &C);
    factor phi_2 ("phi_2", &D);

    factor phi_3 = factor_product(&phi_1, &phi_2, "phi_3");

    factor phi_4 = discrete_marginal(&phi_3, &C, "marg_4");

    bool pass_flag = true;

    if( abs(phi_4.canonical_table[0].g - log(0.1)) > 1e-5 or
            abs(phi_4.canonical_table[1].g - log(0.9)) > 1e-5){
        pass_flag = false;
    }
    return pass_flag;
}

void report( bool test_flag){
    if(test_flag){
        cout << "test passed!\n";
    }
    else{
        cout << "test failed!\n";
    }
}

/**
 * \brief run canonical_form tests, prints which tests passed/failed
 */
void canonical_form_tests(){
    cout << "Beginning canonical form tests.\n";
    cout << "Testing constructor...\n";
    bool constructor_test = canonical_form_constructor_test();
    report(constructor_test);
    cout << "Testing canonical_product...\n";
    bool product_test = canonical_product_test();
    report(product_test);
    cout << "Testing canonical_quotient...\n";
    bool quotient_test = canonical_quotient_test();
    report(quotient_test);
    cout << "Testing canonincal_marginal...\n";
    bool marginal_test = canonical_marginal_test();
    report(marginal_test);
    cout << "Testing canonical_reduction...\n";
    bool reduction_test = canonical_reduction_test();
    report(reduction_test);
    cout << "Tests complete!\n";
}
/**
 * \brief runs variable tests, prints which tests passed/failed
 */

void variable_tests(){
    cout << "Beginning variable tests.\n";
    cout << "Testing constructor...\n";
    bool constructor_test = variable_constructor_test();
    report(constructor_test);
    cout << "Tests complete!\n";
}

/**
 * \brief runs factor tests, prints which tests passed/failed
 */

void factor_tests(){
    cout << "Beginning factor tests.\n";
    cout << "Testing constructor...\n";
    bool constructor_test = factor_constructor_test1();
    report(constructor_test);
    cout << "Testing stride computation...\n";
    bool stride_result = compute_stride_test();
    report(stride_result);
    cout << "Testing assignment computation...\n";
    bool assignment_result = assignment_test();
    report(assignment_result);
    cout << "Testing index computation...\n";
    bool index_result = index_test();
    report(index_result);
    cout << "Testing factor_product: test 1...\n";
    bool product_result1 = factor_product_test1();
    report(product_result1);
    cout << "Testing factor product: test 2...\n";
    bool product_result2 = factor_product_test2();
    report(product_result2);
    cout << "Testing discrete factor marginalization...\n";
    bool marginal_result = discrete_marginal_test1();
    report(marginal_result);
    cout << "Tests complete!";
}

void run_unit_tests(){
    canonical_form_tests();
    variable_tests();
    factor_tests();
}


#endif //TPD_SIGNAL_DECOMPOSITION_UNIT_TESTS_H
