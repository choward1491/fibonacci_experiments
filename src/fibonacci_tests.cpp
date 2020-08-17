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

    /* Simple test to run and compare approaches */
    void test_fibonacci_num() {
        
        // list of input n values to try to see how
        // the two algorithms compare
        std::array<size_t, 7> n_values { 2, 17, 35, 53, 61, 80, 170 };
        
        // loop over value of n and see how the different
        // techniques perform for each input, on average
        size_t num_samples = 100;
        for(auto n: n_values){
            
            double mruntime = 0.0, fdruntime = 0.0;
            {// test the medium post version
                size_t j = 0;
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                
                for(size_t i = 0; i < num_samples; ++i){
                    auto val1 = fibonacci::medium_approach(n);
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
                    auto val2 = fibonacci::fast_doubling(n);
                    j = j + i*i; // help make sure optimizer keeps each iteration of loop
                }
                
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                fdruntime = time_span.count()/static_cast<double>(num_samples);
                std::cout << "Found F(" << n << ") using fast doubling approach in " << fdruntime << " seconds on average" << std::endl;
            }
            
            std::cout << "Fast doubling is " << (mruntime / fdruntime) << " times faster than the medium post approach" << std::endl;
            std::cout << std::endl;
        }// end for loop
    }

}
