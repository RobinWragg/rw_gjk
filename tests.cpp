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

void print_test_name(string name) {
	while (name.size() < 60) {
		name.push_back('_');
	}
	printf("%s", name.c_str());
}

void print_test_result(bool success) {
	printf(success ? "success\n" : "FAIL\n");
}

int main() {
	printf("\n * Running tests for rw_gjk *\n");
	
	Shape shape;
	
	printf("\ntry_make_shape():\n");
	{
		print_test_name("Valid shape with clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1) };
		print_test_result(try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Valid shape with anti-clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(1, 0), v2(1, 1) };
		print_test_result(try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Null pointer argument");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_shape(corners, nullptr));
	}
	
	{
		print_test_name("Invalid concave shape");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, only two corners");
		vector<v2> corners = { v2(0, 0), v2(0, 1) };
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, only one corner");
		vector<v2> corners = { v2(0, 1) };
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, no corners");
		vector<v2> corners = {};
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, corners form a straight line");
		vector<v2> corners = {
			v2(0, 0), v2(1, 0), v2(2, 0), v2(1, 1)
		};
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, duplicate corners");
		vector<v2> corners = {
			v2(0, 0), v2(0, 0), v2(1, 0), v2(0, 1)
		};
		print_test_result(!try_make_shape(corners, &shape));
	}
	
	printf("\nget_overlap_info():\n");
	Shape shape_a, shape_b;
	try_make_shape({v2(-1, -0.5), v2(1, -0.5), v2(0, 1)}, &shape_a);
	try_make_shape({v2(-1, -0.5), v2(1, -0.5), v2(0, 1)}, &shape_b);
	{
		print_test_name("Shapes overlap");
		shape_a.pos = v2(42, 42);
		shape_b.pos = v2(42, 42);
		print_test_result(get_overlap_info(&shape_a, &shape_b).overlapping);
	}
	
	{
		print_test_name("Shapes don't overlap");
		shape_a.pos = v2(-10, 3);
		shape_b.pos = v2(10, 3);
		print_test_result(!get_overlap_info(&shape_a, &shape_b).overlapping);
	}
	
	printf("\n");
	return 0;
}





