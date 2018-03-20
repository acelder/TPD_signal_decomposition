//
// Created by Alexander Elder on 1/26/18.
//
// This file contains matrix decomposition and inversion functions. 

#ifndef TPD_SIGNAL_DECOMPOSITION_LIN_ALG_H
#define TPD_SIGNAL_DECOMPOSITION_LIN_ALG_H

#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <stdexcept>


namespace ublas = boost::numeric::ublas;

namespace lin_alg{
    /*
     * declarations
     */
    // functions
    ublas::matrix<double> eye( int dimension );
    ublas::vector<double> substitution( ublas::matrix<double> A, ublas::vector<double> b, bool backward = true);
    ublas::matrix<double> inverse( ublas::matrix<double> A, bool use_cholesky = false);
    double det( ublas::matrix<double> A, bool use_cholesky = false );

    // classes
    /**
     * Data structure for LU decomposition. Constructor computes decomposition of input matrix.
     */
    class LUP_decomposition{
    public:
        ublas::matrix<double> A; // composite matrix
        ublas::matrix<double> L; // Lower triangular decomposition
        ublas::matrix<double> U; // Upper triangular decomposition
        ublas::matrix<double> P; // Permutation matrix.

        double det_A, det_L, det_U, det_P; // determinants of A,L,U, and P
        LUP_decomposition() = default;
        LUP_decomposition(const ublas::matrix<double> A);
    };


    class cholesky_decomposition : public LUP_decomposition{
    public:

        cholesky_decomposition(ublas::matrix<double> A);
    };

    /**
     * @brief Returns boost identity matrix
     * @param dimension: number of rows and columns of desired identity matrix
     * @return : identity matrix
     */

    ublas::matrix<double> eye( int dimension ){
        ublas::matrix<double> I (dimension, dimension);
        for( int i=0; i<dimension; ++i ){
            for (int j=0; j<dimension; ++j ){
                if( i == j ){
                    I(i,j) = 1;
                }
                else{
                    I(i,j) = 0;
                }
            }
        }

        return I;
    }

    /**
     * @brief Substitution subroutine of Gaussian elimination.
     * @param A : Matrix
     * @param b : Vector
     * @param backward : direction
     * @return : x = A/b
     */
    ublas::vector<double> substitution( ublas::matrix<double> A, ublas::vector<double> b, bool backward ){
        ublas::vector<double> x (A.size1(), 0);
        size_t N = b.size();
        if(backward){ // iterate over rows from bottom to top
            for(long i=N-1; i>=0; --i){ // iterate over lower columns
                double sum = 0;
                for(long j=i+1; j<N; ++j){
                    sum += A(i,j)*x(j);
                }
                x(i) = (b(i) - sum)/A(i,i);
                if( abs(x(i)) < 1e-15 ){ x(i) = 0.0;}// round to zero for numerical stability
            }
        }
        else{
            for(long i=0; i<N; ++i){// iterate over rows from top to bottom
                double sum = 0;
                for(long j=0; j<i; ++j){ // iterate over upper coolumns
                    sum += A(i,j)*x(j);
                }
                x(i) = (b(i) - sum)/A(i,i);
                if( abs(x(i)) < 1e-15 ){ x(i) = 0.0;}// round to zero for numerical stability
            }
        }
        return x;
    }

    ublas::matrix<double> inverse( ublas::matrix<double> A, bool use_cholesky){
        ublas::matrix<double> A_inv (A.size1(), A.size2());
        ublas::matrix<double> I = eye(A.size1());
        if(A_inv.size1() == 1){// handle the 1x1 case
            A_inv(0,0) = 1/A(0,0);
            return A_inv;
        }

        if(!use_cholesky){
            LUP_decomposition LUP (A);
            for(long i=0; i<A.size1(); ++i){
                ublas::vector<double> col_i = column(I,i);
                ublas::vector<double> Pb = ublas::prod(LUP.P, col_i);
                ublas::vector<double> c = substitution(LUP.L, Pb, false);
                column(A_inv, i) = substitution(LUP.U, c, true );
            }
        }
        else{
            cholesky_decomposition decomp (A);
            ublas::matrix<double> L = decomp.L;
            for(long i=0; i<A.size1(); ++i){
                ublas::vector<double> b = column(I,i);
                ublas::vector<double> c = substitution(L, b, false);
                column(A_inv,i) = substitution(trans(L), c, true);
            }
        }


        return A_inv;
    }

    /**
     * @brief Compute determinant of a square matrix
     * @param A : square matrix
     * @param use_cholesky : indicate whether to use cholesky decomposition
     * @return : determinant of A
     */
    double det( ublas::matrix<double> A, bool use_cholesky ){
        double d;
        if( use_cholesky ){
            cholesky_decomposition decomp (A);
            d = decomp.det_A;
        }
        else{
            LUP_decomposition decomp (A);
            d = decomp.det_A;
        }
        return d;
    }

    /**
     * Constructor for LUP decomposition
     * @param A : composite matrix
     */
    LUP_decomposition::LUP_decomposition(const ublas::matrix<double> A) : A(A), U(A) {
        size_t N = A.size1();
        L.resize(N, N);
        P = eye(N);

        det_P = 1;

        // perform Gaussian elimination with partial pivoting on U storing the labmdas in L
        for(int j=0; j<N-1; ++j){
            // select row i >= j such that |a_i,j| = max_[i>=j] {|a_i,j|, |a_i+1,j|,..., |a_N,j|}
            int arg_max = j;
            for(int i=j+1; i<N; ++i){
                if( U(i,j) > U(arg_max,j) ){
                    arg_max = i;
                }
            }
            if( abs(A(arg_max, j)) < 1e-14 ){
                throw std::runtime_error("Singular matrix cannot be decomposed by LU decomposition.");
            }
            // if i != j, interchange rows i and j-NOTE: encode change in permutation and L
            if( arg_max != j ) {
                ublas::vector<double> temp = row(U, j);// may need to use double_matrix data type instead
                row(U, j) = row(U, arg_max);
                row(U, arg_max) = temp;

                temp = row(L, j);// take a closer look at what needs to be permuted
                row(L, j) = row(L, arg_max);
                row(L, arg_max) = temp;

                temp = row(P, j);
                row(P, j) = row(P, arg_max);
                row(P, arg_max) = temp;

                det_P *= -1;
            }
            for(int i=j+1; i<N; ++i){
                double lambda = U(i,j)/U(j,j);
                L(i,j) = lambda;
                for(int k=j; k<N; ++k){
                    U(i,k) = U(i,k) - lambda*U(j,k);
                }
            }
        }

        // set diagonal of L to one
        for(int i=0; i<L.size1(); ++i){
            L(i,i) = 1;
        }

        det_L = 1;
        // compute the determinant of U
        det_U = 1;
        for(int i=0; i<N; ++i){
            det_U *= U(i,i);
        }

        det_A = det_L*det_U/det_P;
    }

    /**
     * Constructor for cholesky decomposition
     * @param A : composite matrix
     */
    cholesky_decomposition::cholesky_decomposition(ublas::matrix<double> composite_A) {
        A = composite_A;
        size_t N = A.size1();
        L.resize(N, N);


        for(long i=0; i<N; ++i){
            // compute diagonal
            double diagonal = A(i,i);
            for(long k=0; k<i; ++k){
                diagonal -= L(i,k)*L(i,k);
            }
            L(i,i) = sqrt(diagonal);
            // compute lower triangle
            for(long j=i+1; j<A.size1(); ++j){
                L(j,i) = A(i,j);
                double element = A(i,j);
                for(long k=0; k<i; ++k){
                    element -=  L(i,k)*L(j,k);
                }
                element = element/L(i,i);
                L(j,i) = element;

            }
        }

        U = trans(L);
        P = eye(N);

        det_L = 1;
        for(int i=0; i<N; ++i){
            det_L *= L(i,i);
        }
        det_U = det_L, det_P = 1;
        det_A = det_U * det_L;
    }

    ublas::vector<double> left_mult(ublas::matrix<double> A, ublas::vector<double> x){
        ublas::vector<double> b = ublas::prod(A, x);
        return b;
    }
}

#endif //TPD_SIGNAL_DECOMPOSITION_LIN_ALG_H
