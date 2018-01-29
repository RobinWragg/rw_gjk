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
	while (name.size() < 70) {
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
	
	printf("\ntry_make_polygon():\n");
	{
		print_test_name("Valid shape with clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1) };
		print_test_result(try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Valid shape with anti-clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(1, 0), v2(1, 1) };
		print_test_result(try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Null pointer argument");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_polygon(corners, nullptr));
	}
	
	{
		print_test_name("Invalid concave shape");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, only two corners");
		vector<v2> corners = { v2(0, 0), v2(0, 1) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, only one corner");
		vector<v2> corners = { v2(0, 1) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, no corners");
		vector<v2> corners = {};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, corners form a straight line");
		vector<v2> corners = {
			v2(0, 0), v2(1, 0), v2(2, 0), v2(1, 1)
		};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid shape, duplicate corners");
		vector<v2> corners = {
			v2(0, 0), v2(0, 0), v2(1, 0), v2(0, 1)
		};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	printf("\nget_overlap_amount():\n");
	const double AMOUNT_TOLERANCE = 0.000001;
	Shape shape_a, shape_b;
	const double identical_polygons_width = 0.2;
	const double identical_polygons_half_width = identical_polygons_width / 2;
	vector<rw_gjk::v2> corners = {
		rw_gjk::v2(-identical_polygons_half_width, -identical_polygons_half_width),
		rw_gjk::v2(identical_polygons_half_width, -identical_polygons_half_width),
		rw_gjk::v2(identical_polygons_half_width, identical_polygons_half_width),
		rw_gjk::v2(-identical_polygons_half_width, identical_polygons_half_width)
	};
	try_make_polygon(corners, &shape_a);
	try_make_polygon(corners, &shape_b);
	
	{
		print_test_name("Identical polygons overlap when both at origin");
		shape_a.pos = v2(0, 0);
		shape_b.pos = v2(0, 0);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		print_test_result(amount.x != 0 || amount.y != 0);
	}
	
	{
		print_test_name("Identical polygons overlap when in same location");
		shape_a.pos = v2(124.32, 74.428);
		shape_b.pos = v2(124.32, 74.428);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		print_test_result(amount.x != 0 || amount.y != 0);
	}
	
	{
		print_test_name("Identical polygons with one below the origin overlap correctly");
		double offset = -0.00198573451;
		shape_a.pos = v2(0, offset);
		shape_b.pos = v2(0, 0);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		
		double offset_amount_difference = fabs(amount.y) - (identical_polygons_width - fabs(offset));
		print_test_result(amount.x == 0 && amount.y > 0
			&& offset_amount_difference > 0
			&& offset_amount_difference < AMOUNT_TOLERANCE);
	}
	
	{
		print_test_name("Identical polygons with one above the origin overlap correctly");
		double offset = 0.0012375095;
		shape_a.pos = v2(0, offset);
		shape_b.pos = v2(0, 0);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		
		double offset_amount_difference = fabs(amount.y) - (identical_polygons_width - fabs(offset));
		print_test_result(amount.x == 0 && amount.y < 0
			&& offset_amount_difference > 0
			&& offset_amount_difference < AMOUNT_TOLERANCE);
	}
	
	{
		print_test_name("Identical polygons with one left of the origin overlap correctly");
		double offset = -0.00198573451;
		shape_a.pos = v2(offset, 0);
		shape_b.pos = v2(0, 0);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		
		double offset_amount_difference = fabs(amount.x) - (identical_polygons_width - fabs(offset));
		print_test_result(amount.x > 0 && amount.y == 0
			&& offset_amount_difference > 0
			&& offset_amount_difference < AMOUNT_TOLERANCE);
	}
	
	{
		print_test_name("Identical polygons with one right of the origin overlap correctly");
		double offset = 0.0025823875955451;
		shape_a.pos = v2(offset, 0);
		shape_b.pos = v2(0, 0);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		
		double offset_amount_difference = fabs(amount.x) - (identical_polygons_width - fabs(offset));
		print_test_result(amount.x < 0 && amount.y == 0
			&& offset_amount_difference > 0
			&& offset_amount_difference < AMOUNT_TOLERANCE);
	}
	
	{
		print_test_name("Polygons don't overlap");
		shape_a.pos = v2(-10, 3);
		shape_b.pos = v2(10, 3);
		v2 amount = get_overlap_amount(&shape_a, &shape_b);
		print_test_result(amount.x == 0 && amount.y == 0);
	}
	
	printf("\n");
	return 0;
}





