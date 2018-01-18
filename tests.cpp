/*
Compile and run these tests in bash with:
g++ -std=c++11 tests.cpp rw_gjk.cpp -o tests && ./tests
*/

#include <cstdio>
#include "rw_gjk.h"

using namespace rw_gjk;

int main() {
	printf("\n * Running tests for rw_gjk *\n\n");
	
	printf("allocShape():\n");
	{
		printf("Valid shape with clockwise winding - ");
		vector<v2> corners = {v2(0, 0), v2(0, 1), v2(1, 1)};
		Shape *shape = allocShape(corners);
		shape == nullptr ? printf("FAILED\n") : printf("success\n");
	}
	
	{ // valid shape with anti-clockwise winding
		vector<v2> corners = {v2(0, 0), v2(1, 0), v2(1, 1)};
		Shape *shape = allocShape(corners);
		// assert(shape != nullptr);
	}
	
	{ // concave shape
		vector<v2> corners = {
			v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9)};
		Shape *shape = allocShape(corners);
		// assert(shape == nullptr);	
	}
	
	{ // corners form a straight line
		vector<v2> corners = {
			v2(0, 0), v2(1, 0), v2(2, 0), v2(1, 1)
		};
		Shape *shape = allocShape(corners);
		// assert(shape == nullptr);
	}
	
	{ // duplicate corners
		vector<v2> corners = {
			v2(0, 0), v2(0, 0), v2(1, 0), v2(0, 1)
		};
		Shape *shape = allocShape(corners);
		// assert(shape == nullptr);
	}
	
	printf("\n");
	return 0;
}