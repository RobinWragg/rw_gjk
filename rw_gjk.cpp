/*
rw_gjk

A 2D collision detection and resolution library.
Based on the Gilbert–Johnson–Keerthi distance algorithm (GJK).
Written by Robin Wragg.

See end of file for license.
*/

// TODO: check for NaN?
// TODO: remove occurrences of while (true).
// TODO: remove asserts and handle those states in a graceful way.
// TODO: add caching of rotations.
// TODO: add error strings.
// TODO: calculate TINY_NUMBER appropriately. It has something to do with adding largest_radius in shapes_are_overlapping(). TINY_NUMBER is basically used as half the thickness of a line. It should be very small but never small enough to cause IEEE-float-related problems.

#include <vector>
#include <cmath>
#include <cassert>
#include <string>

namespace rw_gjk {
	#include "vectors.cpp"
	
	using namespace std;
	
	const double TINY_NUMBER = 0.0000001;
	const int LOOP_LIMIT = 1024;
	v2 ORIGIN = v2(0, 0);
	
	struct Shape {
		v2 pos;
		double radius;
		bool is_circle;
		
		float angle;
		vector<v2> corners;
	};
	
	bool contains_duplicates(vector<v2> vertices) {
		if (vertices.size() == 1) return false;
		for (int s0 = 0; s0 < vertices.size(); s0++) {
			for (int s1 = s0+1; s1 < vertices.size(); s1++) {
				if (vertices[s0] == vertices[s1]) return true;
			}
		}
		return false;
	}
	
	vector<v2> get_ordered_convex_corners(vector<v2> corners) {
		assert(!contains_duplicates(corners));
		vector<v2> convex_corners;
		
		// add the leftmost corner
		v2 leftmost_corner = corners.front();
		for (auto &corner: corners) {
			if (corner.x < leftmost_corner.x) leftmost_corner = corner;
		}
		convex_corners.push_back(leftmost_corner);
		
		// build a convex shape out of the corners, in any order.
		v2 search_dir = v2(0, 1);
		while (true) {
			// find the next corner
			double highest_corner_dot = -INFINITY;
			v2 best_corner;
			for (auto &corner: corners) {
				if (corner == convex_corners.back()) continue;
				
				v2 corner_difference = corner - convex_corners.back();
				double corner_dot = dot(search_dir, corner_difference.normalised_or_zero());
				
				if (corner_dot >= 0 && corner_dot > highest_corner_dot) {
					highest_corner_dot = corner_dot;
					best_corner = corner;
				}
			}
			
			// if a valid corner was found
			if (highest_corner_dot >= 0) {
				// check if it's the same as the first corner
				if (best_corner == convex_corners[0]) break; // the convex shape is complete
				else {
					// add the corner to the shape and update the search direction
					search_dir = (best_corner - convex_corners.back()).normalised_or_zero();
					convex_corners.push_back(best_corner);
				}
			} else {
				// no valid corner was found in the current search direction, so rotate it 90 degrees
				search_dir = search_dir.right_normal_or_zero();
			}
		}
		
		assert(convex_corners.size() <= corners.size());
		assert(!contains_duplicates(convex_corners));
		return convex_corners;
	}
	
	bool is_convex_and_ordered(vector<v2> corners) {
		if (corners.size() < 3) return true;
		
		auto ordered_corners = get_ordered_convex_corners(corners);
		if (ordered_corners.size() != corners.size()) return false;
		
		int start_index = 0;
		while (corners[start_index] != ordered_corners[0]) start_index++;
		
		bool first_winding_failed = false;
		
		for (int i = 0; i < ordered_corners.size(); i++) {
			int corner_index = (start_index+i) % corners.size();
			if (corners[corner_index] != ordered_corners[i]) {
				first_winding_failed = true;
				break;
			}
		}
		
		if (first_winding_failed) {
			for (int i = 0; i < ordered_corners.size(); i++) {
				int corner_index = (start_index-i) % corners.size();
				if (corners[corner_index] != ordered_corners[i]) {
					first_winding_failed = true;
					break;
				}
			}
		}
		
		return true;
	}
	
	bool is_valid(vector<v2> corners) {
		bool dup = contains_duplicates(corners);
		bool conv_ordered = is_convex_and_ordered(corners);
		
		assert(!dup);
		assert(conv_ordered);
		
		return (!dup) && conv_ordered;
	}
	
	void make_circle(double radius, Shape *shape_out) {
		shape_out->radius = radius;
		shape_out->pos = ORIGIN;
		shape_out->angle = 0;
		shape_out->is_circle = true;
	}
	
	bool try_make_polygon(vector<v2> corners, Shape *shape_out) {
		if (corners.size() < 3) {
			return false;
		}
		
		if (contains_duplicates(corners)) {
			return false;
		}
		
		Shape shape;
		
		// check that the corners are convex and also order them
		{
			auto ordered_convex_corners = get_ordered_convex_corners(corners);
			assert(is_valid(ordered_convex_corners));
			
			if (corners.size() != ordered_convex_corners.size()) {
				return false;
			} else {
				shape.corners = ordered_convex_corners;
			}
		}
		
		// set radius
		for (auto &corner: shape.corners) {
			if (corner.length() > shape.radius) shape.radius = corner.length();
		}
		
		*shape_out = shape;
		shape_out->pos = ORIGIN;
		shape_out->angle = 0;
		shape_out->is_circle = false;
		
		return true;
	}
	
	v2 get_minkowski_diffed_corner(Shape *shape, Shape *other_shape, v2 direction) {
		auto get_corner_of_rotated_shape = [](Shape *shape, v2 direction) {
			if (shape->is_circle) {
				return direction.normalised_or_zero() * shape->radius;
			} else {
				v2 best_corner;
				double best_dot = -INFINITY;
				
				for (const auto &corner: shape->corners) {
					v2 rotated_corner = corner.rotated(shape->angle);
					double new_dot = dot(rotated_corner, direction);
					
					if (new_dot > best_dot) {
						best_corner = rotated_corner;
						best_dot = new_dot;
					}
				}
				
				return best_corner;
			}
		};
		
		v2 shape_world_corner = shape->pos + get_corner_of_rotated_shape(shape, direction);
		v2 other_shape_world_corner = other_shape->pos + get_corner_of_rotated_shape(other_shape, -direction);
		return shape_world_corner - other_shape_world_corner;
	}
	
	bool origin_is_between_points(v2 a, v2 b) {
		v2 ao = ORIGIN - a;
		v2 bo = ORIGIN - b;
		v2 ab = b - a;
		v2 ba = a - b;
		return (dot(ao, ab) >= 0) && (dot(bo, ba) >= 0);
	}
	
	enum Origin_Location {
		OL_ON_LINE,
		OL_NEAR_LINE,
		OL_NEAR_CORNER_0,
		OL_NEAR_CORNER_1
	};
	
	bool improve_2_simplex(vector<v2> &simplex, v2 &search_dir) {
		assert(simplex.size() == 2);
		
		Origin_Location origin_location;
		{
			if (origin_is_between_points(simplex[0], simplex[1])) {
				v2 line_normal = (simplex[1] - simplex[0]).right_normal_or_zero();
				double origin_distance_from_line = dot(line_normal, ORIGIN - simplex[0]);
				
				if (fabs(origin_distance_from_line) <= TINY_NUMBER) origin_location = OL_ON_LINE;
				else origin_location = OL_NEAR_LINE;
				
			} else if (dot(simplex[1] - simplex[0], ORIGIN - simplex[0]) <= 0) {
				origin_location = OL_NEAR_CORNER_0;
			} else {
				assert(dot(simplex[0] - simplex[1], ORIGIN - simplex[1]) <= 0);
				origin_location = OL_NEAR_CORNER_1;
			}
		}
		
		switch (origin_location) {
			case OL_ON_LINE: {
				return true;
				break;
			}
			case OL_NEAR_LINE: {
				// The simplex is correct. Search on the side of the 2-simplex that contains the origin.
				v2 simplex_line = simplex[1] - simplex[0];
				v2 direction_to_origin = ORIGIN - simplex[0];
				search_dir = simplex_line.normal_in_direction_or_zero(direction_to_origin);
				break;
			}
			case OL_NEAR_CORNER_0: {
				simplex = {simplex[0]};
				search_dir = (ORIGIN - simplex[0]).normalised_or_zero();
				break;
			}
			case OL_NEAR_CORNER_1: {
				simplex = {simplex[1]};
				search_dir = (ORIGIN - simplex[1]).normalised_or_zero();
				break;
			}
			default: {
				assert(false);
				break;
			}
		}
		
		assert(search_dir.length() > 0.9999 && search_dir.length() < 1.0001);
		return false;
	}
	
	// returns true when the simplex contains the origin.
	bool improve_simplex(vector<v2> &simplex, v2 &search_dir) {
		if (simplex.size() == 3) {
			// cache some basic vectors
			v2 ab = simplex[1] - simplex[0];
			v2 bc = simplex[2] - simplex[1];
			v2 ca = simplex[0] - simplex[2];
			
			v2 ab_normal_away_from_c = ab.normal_in_direction_or_zero(ca);
			v2 bc_normal_away_from_a = bc.normal_in_direction_or_zero(ab);
			v2 ca_normal_away_from_b = ca.normal_in_direction_or_zero(bc);
			
			// find which side of the triangle the origin is on, or if it's inside it.
			if (dot(ab_normal_away_from_c, ORIGIN - simplex[0]) > 0) {
				simplex = {simplex[0], simplex[1]};
				return improve_2_simplex(simplex, search_dir);
			} else if (dot(bc_normal_away_from_a, ORIGIN - simplex[1]) > 0) {
				simplex = {simplex[1], simplex[2]};
				return improve_2_simplex(simplex, search_dir);
			} else if (dot(ca_normal_away_from_b, ORIGIN - simplex[2]) > 0) {
				simplex = {simplex[2], simplex[0]};
				return improve_2_simplex(simplex, search_dir);
			} else {
				assert(is_valid(simplex));
				return true; // the origin is inside the simplex.
			}
		}
		
		assert(simplex.size() == 2);
		return improve_2_simplex(simplex, search_dir);
	}
	
	bool shapes_are_overlapping(
		Shape *shape_a, Shape *shape_b,
		vector<v2> *simplex_out = nullptr // This is only used internally.
		) {
		
		v2 search_direction = v2(1, 0); // this can be anything except zero, I think? TODO
		vector<v2> simplex = {get_minkowski_diffed_corner(shape_a, shape_b, search_direction)};
		
		search_direction = ORIGIN - simplex[0]; // search toward the origin
		
		while (true) {
			assert(simplex.size() < 3);
			simplex.push_back(get_minkowski_diffed_corner(shape_a, shape_b, search_direction));
			
			double dot_result = dot(simplex.back(), search_direction);
			if (dot_result <= TINY_NUMBER) {
				return false;
			}
			
			if (improve_simplex(simplex, search_direction)) {
				if (simplex_out != nullptr) *simplex_out = simplex;
				return true;	
			}
		}
	}
	
	v2 get_overlap_amount(Shape *shape_a, Shape *shape_b) {
		return v2(0, 0); // NO-OP
		
		// NOTE NOTE NOTE NOTE: These comments don't take into account origins directly on lines.
		// begin infinite loop
			// get simplex line closest to origin
			
			// get the outer normal of that line and get minkowski diffed corner in that direction
			
			// if the new corner is almost identical to one of the points that made the simplex line,
				// find the point on the line that is closest to the origin
				// the difference between the origin and that point is the overlap amount.
			// else add the new corner to the simplex, turning the existing line into two.
		// (back to beginning of loop)
		
		
		
	}
}

/*
MIT License

Copyright (c) 2018 Robin Wragg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
