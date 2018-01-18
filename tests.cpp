/*
Compile and run these tests in bash with:
g++ -std=c++11 tests.cpp rw_gjk.cpp -o tests && ./tests
*/

#include <cstdio>
#include "rw_gjk.h"

int main() {
	printf("\n * Running tests for rw_gjk *\n\n");
	
	auto s = rw_gjk::allocShape({rw_gjk::v2(1, 1), rw_gjk::v2(1, 2), rw_gjk::v2(-10, 1)});
	
	printf("\n");
	return 0;
}