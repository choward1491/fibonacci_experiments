//
//  fibonacci_tests.cpp
//  tinkering
//
//  Created by Christian Howard on 8/17/20.
//  Copyright © 2020 Christian Howard. All rights reserved.
//

#include "fibonacci_tests.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <chrono>
#include <iostream>
#include <armadillo>
#include <string>

namespace fibonacci {

    /* Medium post related code */
    double log_binomial(size_t n, size_t k) {
        double num = 0.0, denom = 0.0;
        for(size_t i = 0; i < (n-k); ++i){
            num += std::log(static_cast<double>(n-i));
            denom += std::log(static_cast<double>(i+1));
        }
        return (num - denom);
    }

    __uint128_t medium_approach( size_t n ) {
        __uint128_t fn = 0;
        
        std::vector<double> numbers(n+1);
        double c_max = 0.0;
        
        // compute the coefficients using the log form
        // and find the maximum value
        for(size_t i = 1; i <= n; i += 2){
            double logbinom = log_binomial(n, i);
            double ci = logbinom + (i-1) * std::log(std::sqrt(5.0)) - (n-1) * std::log(2.0);
            c_max = std::max(c_max, ci);
            numbers[i] = ci;
        }
        
        // compute the normalized sum
        double sum = 0.0;
        for(size_t i = 1; i <= n; i += 2 ){
            sum += exp(numbers[i] - c_max);
        }
        
        // compute the resulting integer
        quad exponential = std::exp(static_cast<quad>(c_max));
        quad value = std::round(exponential * static_cast<quad>(sum));
        fn = static_cast<__uint128_t>(value);
        
        // return the result
        return fn;
    }

    namespace version2 {
    
        arma::vec log_binomial2(size_t n, const arma::vec& indices){
            arma::vec r = arma::regspace(0, n-1) + 1;
            r = arma::log(r);
            double s = arma::sum(r);
            r = arma::cumsum(r);
            arma::vec z(n+1, arma::fill::zeros);
            z.elem(arma::regspace<arma::uvec>(1, n)) = r;
            arma::vec z1 = z.elem(arma::regspace<arma::uvec>(n, -1, 0));
            z = z + z1;
            z = s - z;
            return z.elem(arma::regspace<arma::uvec>(1, 2, n));
        }
    
        __uint128_t medium_approach( size_t n ) {
            __uint128_t fn = 0;
            
            arma::vec ks = arma::regspace(0, (n+1)/2 - 1);
            arma::vec odds = 2.0 * ks + 1.0;
            arma::vec coefs = log_binomial2(n, odds);
            arma::vec terms = std::log(std::sqrt(5.0)) * odds;
            arma::vec res = coefs + terms;
            res = res - (std::log(2.0)*(n-1) + std::log(std::sqrt(5.0)));
            double m = arma::max(res);
            res = res - m;
            res = arma::exp(res);
            double sum = arma::sum(res);
            
            // compute the resulting integer
            quad exponential = std::exp(static_cast<quad>(m));
            quad value = std::round(exponential * static_cast<quad>(sum));
            fn = static_cast<__uint128_t>(value);
            
            // return the result
            return fn;
        }
    }

    namespace biginput {
        big_int medium_approach( size_t n ) {
            
            big_int fn;
            
            arma::vec ks = arma::regspace(0, (n+1)/2 - 1);
            arma::vec odds = 2.0 * ks + 1.0;
            arma::vec coefs = version2::log_binomial2(n, odds);
            arma::vec terms = std::log(std::sqrt(5.0)) * odds;
            arma::vec res = coefs + terms;
            res = res - (std::log(2.0)*(n-1) + std::log(std::sqrt(5.0)));
            double m = arma::max(res);
            res = res - m;
            res = arma::exp(res);
            double sum = arma::sum(res);
            
            // do the conversion
            big_float exponential = boost::multiprecision::exp(static_cast<big_float>(m));
            big_float value = boost::multiprecision::round(exponential * static_cast<big_float>(sum));
            fn = static_cast<big_int>(value);
            
            // return the result
            return fn;
            
        }
    }

    /* Fast doubling related code */
    using mat22 = std::array<__uint128_t, 4>;
    using vec2  = std::array<__uint128_t, 2>;
    void mat_vec_inplace( const mat22& A, vec2& x){
        // computing out = A * B
        // assume row order matrices
        auto tmp0 = A[0]*x[0] + A[1]*x[1];
        auto tmp1 = A[2]*x[0] + A[3]*x[1];
        x[0] = tmp0;
        x[1] = tmp1;
    }
    void mat_square( const mat22& A, mat22& out){
        // computing out = A * A
        // assume row order matrices
        out[0] = A[0]*A[0] + A[1]*A[2];
        out[1] = A[0]*A[1] + A[1]*A[3];
        out[2] = A[2]*A[0] + A[3]*A[2];
        out[3] = A[2]*A[1] + A[3]*A[3];
    }

    std::string to_string(__uint128_t num) {
        auto tenPow19 = 10000000000000000000;
        std::string str;
        do {
            uint64_t digits = num % tenPow19;
            auto digitsStr = std::to_string(digits);
            auto leading0s = (digits != num) ? std::string(19 - digitsStr.length(), '0') : "";
            str = leading0s + digitsStr + str;
            num = (num - digits) / tenPow19;
        } while (num != 0);
        return str;
    }

    __uint128_t fast_doubling( size_t n ) {
        
        __uint128_t fn = 0;
        
        // return trivial base cases
        if( n == 0 ){ return 0; }
        if( n == 1 ){ return 1; }
        
        // initialize the list of pre-computed matrices
        --n;
        std::vector<mat22> computed_mats(std::log2(n)+1);
        size_t size = computed_mats.size();
        
        // initialize first matrix
        mat22& A = computed_mats[0];
        A[0] = 1; A[1] = 1;
        A[2] = 1; A[3] = 0;
        
        // compute matrices using previous matrix
        for(size_t i = 1; i < size; ++i){
            mat_square(computed_mats[i-1], computed_mats[i]);
        }// end for i
        
        // compute the final result
        vec2 x;
        x[0] = 1; x[1] = 0;
        size_t pow = 1;
        for(size_t i = 0; i < size; ++i){
            if( (pow & n) ){
                mat_vec_inplace(computed_mats[i], x);
            }
            pow *= 2;
        }// end for i
        
        // set the fn value and return it
        fn = x[0];
        return fn;
    }

    namespace biginput {
        using mat22 = std::array<big_int, 4>;
        using vec2  = std::array<big_int, 2>;
        void mat_vec_inplace( const mat22& A, vec2& x){
            // computing out = A * B
            // assume row order matrices
            big_int tmp0 = (A[0]*x[0]) + (A[1]*x[1]);
            big_int tmp1 = (A[2]*x[0]) + (A[3]*x[1]);
            x[0] = tmp0;
            x[1] = tmp1;
        }
        void mat_square( const mat22& A, mat22& out){
            // computing out = A * A
            // assume row order matrices
            out[0] = A[0]*A[0] + A[1]*A[2];
            out[1] = A[0]*A[1] + A[1]*A[3];
            out[2] = A[2]*A[0] + A[3]*A[2];
            out[3] = A[2]*A[1] + A[3]*A[3];
        }
    
        big_int fast_doubling( size_t n ) {
            big_int fn = 0;
            
            // return trivial base cases
            if( n == 0 ){ return 0; }
            if( n == 1 ){ return 1; }
            
            // initialize the list of pre-computed matrices
            --n;
            std::vector<mat22> computed_mats(std::log2(n)+1);
            size_t size = computed_mats.size();
            
            // initialize first matrix
            mat22& A = computed_mats[0];
            A[0] = 1; A[1] = 1;
            A[2] = 1; A[3] = 0;
            
            // compute matrices using previous matrix
            for(size_t i = 1; i < size; ++i){
                mat_square(computed_mats[i-1], computed_mats[i]);
            }// end for i
            
            // compute the final result
            vec2 x;
            x[0] = 1; x[1] = 0;
            size_t pow = 1;
            for(size_t i = 0; i < size; ++i){
                if( (pow & n) ){
                    mat_vec_inplace(computed_mats[i], x);
                }
                pow *= 2;
            }// end for i
            
            // set the fn value and return it
            fn = x[0];
            return fn;
        }
    }

    /* Simple test to run and compare approaches */
    void test_fibonacci_num() {
        
        // list of input n values to try to see how
        // the two algorithms compare
        std::array<size_t, 8> n_values { 2, 17, 35, 53, 61, 80, 100, 200 };
        
        // loop over value of n and see how the different
        // techniques perform for each input, on average
        size_t num_samples = 200;
        for(auto n: n_values){
            
            double mruntime = 0.0, fdruntime = 0.0;
            __uint128_t val1, val2;
            {// test the medium post version
                size_t j = 0;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                
                for(size_t i = 0; i < num_samples; ++i){
                    val1 = fibonacci::version2::medium_approach(n);
                    j = j + i*i; // help make sure optimizer keeps each iteration of loop
                }
                
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                mruntime = time_span.count()/static_cast<double>(num_samples);
                std::cout << "Found F(" << n << ") using medium post approach in " << mruntime << " seconds on average" << std::endl;
            }
            {// test the naive approach
                size_t j = 0;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                
                for(size_t i = 0; i < num_samples; ++i){
                    val2 = fibonacci::fast_doubling(n);
                    j = j + i*i; // help make sure optimizer keeps each iteration of loop
                }
                
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                fdruntime = time_span.count()/static_cast<double>(num_samples);
                std::cout << "Found F(" << n << ") using fast doubling approach in " << fdruntime << " seconds on average" << std::endl;
            }
            
            quad difference = 0.0;
            if( val1 >= val2 ){
                difference = static_cast<quad>( val1 - val2 );
            }else{
                difference = static_cast<quad>( val2 - val1 );
            }
            double percent_different =  difference / static_cast<quad>(val2);
            percent_different = std::abs(percent_different);
            
            std::cout << "Fast doubling is " << (mruntime / fdruntime) << " times faster than the medium post approach" << std::endl;
            std::cout << "The error in medium approach is " << percent_different << "%" << std::endl;
            std::cout << std::endl;
        }// end for loop
    }

    void test_exact_fibonacci_bignum() {
        
        // list of input n values to try to see how
        // the two algorithms compare
        std::vector<size_t> n_values { 10000, 25000, 50000, 100000, 200000, 400000, 800000, 1600000, 2000000 };
        
        size_t num_samples = 5;
        for(auto n: n_values){
            {// test the exact big num version
                big_int out;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                
                for(size_t i = 0; i < num_samples; ++i){
                    out = fibonacci::biginput::fast_doubling(n);
                }
                
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                double mruntime = time_span.count()/static_cast<double>(num_samples);
                std::cout << "In " << mruntime << " seconds on avg, found F(" << n << ")" << std::endl;
                //std::cout << "F(" << n << ") = " << ".. yeah, too big" << std::endl; //out << std::endl;
            }
        }
    }

    void test_approx_fibonacci_bignum() {
        
        // list of input n values to try to see how
        // the two algorithms compare
        std::vector<size_t> n_values { 10000, 25000, 50000, 100000, 200000, 400000, 800000, 1600000, 2000000 };
        
        size_t num_samples = 5;
        for(auto n: n_values){
            {// test the exact big num version
                big_int out;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                
                for(size_t i = 0; i < num_samples; ++i){
                    out = fibonacci::biginput::medium_approach(n);
                }
                
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                double mruntime = time_span.count()/static_cast<double>(num_samples);
                std::cout << "In " << mruntime << " seconds on avg, found F(" << n << ")" << std::endl;
                //std::cout << "F(" << n << ") = " << ".. yeah, too big" << std::endl; //out << std::endl;
            }
        }
            
    }

}
