//
//  main.cpp
//  fibonacci_experiments
//
//  Created by Christian Howard on 8/17/20.
//  Copyright © 2020 Christian Howard. All rights reserved.
//

#include <iostream>
#include "fibonacci_tests.hpp"

int main(int argc, const char * argv[]) {
    //fibonacci::test_fibonacci_num();
    size_t n = 30000; // F(186) is the last number that can fit into __uint128_t
    auto fn1 = fibonacci::fast_doubling( n );
    std::cout << "F(n) = " << fibonacci::to_string(fn1) << std::endl;
    
    auto fn2 = fibonacci::version2::medium_approach( n );
    std::cout << "F(n) = " << fibonacci::to_string(fn2) << std::endl;
    
    auto fn3 = fibonacci::biginput::fast_doubling( n );
    std::cout << "F(n) = " << fn3 << std::endl;
    fibonacci::test_exact_fibonacci_bignum();
    return 0;
}
