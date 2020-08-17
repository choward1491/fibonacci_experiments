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

namespace fibonacci {
    using quad = long double;

    __uint128_t medium_approach( size_t n );
    __uint128_t fast_doubling( size_t n );

    void test_fibonacci_num();

}



#endif /* fibonacci_tests_hpp */
