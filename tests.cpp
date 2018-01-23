/*
Compile and run these tests in bash with:
g++ -std=c++11 tests.cpp -o tests && ./tests

I personally like to use:
clear && echo Compiling... && g++ -std=c++11 tests.cpp -o tests && ./tests && rm tests
to clear the terminal beforehand and delete the executable after I'm done with it.
*/

#include <cstdio>
#include <string>

#include "rw_gjk.h"
#define RW_GJK_IMPLEMENTATION
#include "rw_gjk.h" // Included twice to test the include guards.

using namespace rw_gjk;

void printTestName(string name) {
	while (name.size() < 60) {
		name.push_back('_');
	}
	printf("%s", name.c_str());
}

void printTestResult(bool success) {
	printf(success ? "success\n" : "FAIL\n");
}

int main() {
	printf("\n * Running tests for rw_gjk *\n\n");
	
	Shape shape;
	
	printf("tryMakeShape():\n");
	{
		printTestName("Valid shape with clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1) };
		printTestResult(tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Valid shape with anti-clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(1, 0), v2(1, 1) };
		printTestResult(tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Null pointer argument");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		printTestResult(!tryMakeShape(corners, nullptr));
	}
	
	{
		printTestName("Invalid concave shape");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Invalid shape, only two corners");
		vector<v2> corners = { v2(0, 0), v2(0, 1) };
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Invalid shape, only one corner");
		vector<v2> corners = { v2(0, 1) };
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Invalid shape, no corners");
		vector<v2> corners = {};
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Invalid shape, corners form a straight line");
		vector<v2> corners = {
			v2(0, 0), v2(1, 0), v2(2, 0), v2(1, 1)
		};
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	{
		printTestName("Invalid shape, duplicate corners");
		vector<v2> corners = {
			v2(0, 0), v2(0, 0), v2(1, 0), v2(0, 1)
		};
		printTestResult(!tryMakeShape(corners, &shape));
	}
	
	printf("\n");
	return 0;
}