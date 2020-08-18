//
//  fibonacci_tests.hpp
//  tinkering
//
//  Created by Christian Howard on 8/17/20.
//  Copyright © 2020 Christian Howard. All rights reserved.
//

#ifndef fibonacci_tests_hpp
#define fibonacci_tests_hpp

#include <stddef.h>
#include <string>
#include <boost/multiprecision/gmp.hpp>

namespace fibonacci {
    using big_int = boost::multiprecision::mpz_int;
    using quad = long double;

    std::string to_string(__uint128_t num);
    __uint128_t medium_approach( size_t n );
    namespace version2 {
    __uint128_t medium_approach( size_t n );
    }
    __uint128_t fast_doubling( size_t n );

    namespace biginput {
        big_int fast_doubling( size_t n );
    }

    void test_fibonacci_num();
    void test_exact_fibonacci_bignum();

}



#endif /* fibonacci_tests_hpp */
