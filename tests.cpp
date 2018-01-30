/*
Compile and run these tests in bash with:
g++ -std=c++11 tests.cpp -o tests && ./tests

I personally like to use:
clear && echo Compiling... && g++ -std=c++11 tests.cpp -o tests && ./tests && rm tests
to clear the terminal beforehand and delete the executable after I'm done with it.
*/

#include <cstdio>
#include <string>

#include "rw_gjk.cpp"

using namespace rw_gjk;

void print_test_name(string name) {
	while (name.size() < 70) {
		name.push_back('_');
	}
	printf("%s", name.c_str());
}

int num_failed_tests = 0;

void print_test_result(bool success) {
	printf(success ? "success\n" : "FAIL\n");
	if (!success) num_failed_tests++;
}

int main() {
	printf("\n * Running tests for rw_gjk *\n");
	
	Shape shape;
	
	printf("\nis_convex_and_ordered():\n");
	{
		print_test_name("Valid test A");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1) };
		print_test_result(is_convex_and_ordered(corners));
	}
	
	{
		print_test_name("Valid test B");
		vector<v2> corners = { v2(1, 1), v2(0, 0), v2(0, 1) };
		print_test_result(is_convex_and_ordered(corners));
	}
	
	{
		print_test_name("Valid test C");
		vector<v2> corners = {
			v2(0.38129108817537805, 0.0073923092139486363),
			v2(-0.48871174908274423, 0.034026436793289747),
			v2(-0.078922328003752942, -0.41132716476704667)
		};
	  print_test_result(is_convex_and_ordered(corners));
	}
	
	{
		print_test_name("Valid test D");
		vector<v2> corners = {
			v2(-0.48871174908274423, 0.034026436793289747),
			v2(-0.078922328003752942, -0.41132716476704667),
			v2(0.38129108817537805, 0.0073923092139486363)
		};
	  print_test_result(is_convex_and_ordered(corners));
	}
	
	{
		print_test_name("Valid test E");
		vector<v2> corners = {
			v2(-0.078922328003752942, -0.41132716476704667),
			v2(-0.48871174908274423, 0.034026436793289747),
			v2(0.38129108817537805, 0.0073923092139486363)
		};
	  print_test_result(is_convex_and_ordered(corners));
	}
	
	{
		print_test_name("Valid test F");
		vector<v2> corners = { v2(0, 0), v2(1, 1), v2(0, 1) };
		print_test_result(is_convex_and_ordered(corners));
	}
	
	printf("\ntry_make_polygon():\n");
	{
		print_test_name("Valid polygon with clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1) };
		print_test_result(try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Valid polygon with anti-clockwise winding");
		vector<v2> corners = { v2(0, 0), v2(1, 0), v2(1, 1) };
		print_test_result(try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Null pointer argument");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_polygon(corners, nullptr));
	}
	
	{
		print_test_name("Invalid concave polygon");
		vector<v2> corners = { v2(0, 0), v2(0, 1), v2(1, 1), v2(0.1, 0.9) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid polygon, only two corners");
		vector<v2> corners = { v2(0, 0), v2(0, 1) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid polygon, only one corner");
		vector<v2> corners = { v2(0, 1) };
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid polygon, no corners");
		vector<v2> corners = {};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid polygon, corners form a straight line");
		vector<v2> corners = {
			v2(0, 0), v2(1, 0), v2(2, 0), v2(1, 1)
		};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		print_test_name("Invalid polygon, duplicate corners");
		vector<v2> corners = {
			v2(0, 0), v2(0, 0), v2(1, 0), v2(0, 1)
		};
		print_test_result(!try_make_polygon(corners, &shape));
	}
	
	{
		printf("\nshapes_are_overlapping():\n");
		// const double AMOUNT_TOLERANCE = 0.000001;
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
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Identical polygons overlap when in same location");
			v2 location = v2(124.32, 74.428);
			shape_a.pos = location;
			shape_b.pos = location;
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Identical polygons with one below the origin overlap");
			shape_a.pos = v2(0, -0.00198573451);
			shape_b.pos = v2(0, 0);
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Identical polygons with one above the origin overlap");
			shape_a.pos = v2(0, 0.0012375095);
			shape_b.pos = v2(0, 0);
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Identical polygons with one left of the origin overlap");
			double offset = -0.00198573451;
			shape_a.pos = v2(-0.00198573451, 0);
			shape_b.pos = v2(0, 0);
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Identical polygons with one right of the origin overlap");
			shape_a.pos = v2(0.0025823875955451, 0);
			shape_b.pos = v2(0, 0);
			print_test_result(shapes_are_overlapping(&shape_a, &shape_b));
		}
		
		{
			print_test_name("Polygons don't overlap");
			shape_a.pos = v2(-10, 3);
			shape_b.pos = v2(10, 3);
			print_test_result(!shapes_are_overlapping(&shape_a, &shape_b));
		}
	} // end shapes_are_overlapping()
	
	{
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
	} // end get_overlap_amount()
	
	printf("\nNumber of failed tests: %i\n\n", num_failed_tests);
	return 0;
}





