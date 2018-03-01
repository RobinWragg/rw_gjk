/*
rw_gjk

A detection and resolution library for 2D shapes.
Based on the Gilbert–Johnson–Keerthi distance algorithm (GJK).
Written by Robin Wragg.

Many GJK libraries and most popular GJK tutorials assume that a line is infinitely thin, i.e. points
that fall directly on the line are not handled correctly. This fails robustness tests because GJK
requires a lot of testing of whether a point is on one side of a line or the other. rw_gjk avoids
this issue by treating lines as thin strips that have area, and points that land in those areas are
explicitly handled. In the case of the origin on a simplex line, this is treated as the origin being
inside the simplex. This adds negligible running cost for the vast majority of situations, and the
accuracy of overlap detection/resolution is not affected because the line thickness is intelligently
set based on a combination of the size of the shapes being tested and IEEE float error margins.

See end of file for license.
*/

// TODO: check for NaN?
// TODO: remove asserts and handle those states in a graceful way.
// TODO: calculate TINY_NUMBER appropriately. It has something to do with adding largest_radius in shapes_are_overlapping(). TINY_NUMBER is basically used as half the thickness of a line. It should be very small but never small enough to cause IEEE-float-related problems.
// TODO: get_minkowski_diffed_corner() can return in-line corners. Investigate.

#include <vector>
#include <cmath>
#include <cassert>
#include <string>

namespace rw_gjk {
	#include "vectors.cpp"
	
	using namespace std;
	
	const double LINE_THICKNESS = 0.0000001;
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
	
	bool is_convex(vector<v2> corners) {
		assert(corners.size() > 2);
		
		// return false if any three points are collinear, i.e. form a straight line.
		for (int c0 = 0; c0 < corners.size(); c0++) {
			for (int c1 = 0; c1 < corners.size(); c1++) {
				if (c1 == c0) continue;
				for (int c2 = 0; c2 < corners.size(); c2++) {
					if (c2 == c0 || c2 == c1) continue;
					
					v2 a = (corners[c1] - corners[c0]).normalised_or_0();
					v2 b = (corners[c2] - corners[c1]).normalised_or_0();
					
					if (dot(a, b) == 1) return false;
				}
			}
		}
		
		/*
		The rest of this function detects concavity by finding the convex hull of all corners using the
		gift wrapping algorithm (https://en.wikipedia.org/wiki/Gift_wrapping_algorithm). If there are
		more corners than those that constructed the convex hull, that means the remaining corners are
		concave.
		*/
		vector<v2> convex_hull;
		
		// start with the leftmost corner. If two corners are equally leftmost, choose the upper one.
		{
			v2 leftmost_corner = corners.front();
			for (auto &corner: corners) {
				if (corner.x < leftmost_corner.x
					|| (corner.x == leftmost_corner.x && corner.y > leftmost_corner.y)) {
					leftmost_corner = corner;
				}
			}
			convex_hull.push_back(leftmost_corner);
		}
		
		// build a convex shape out of the corners, in any order.
		v2 search_direction = v2(0, 1);
		while (true) {
			// find the corner that is closest to parallel with the search direction.
			int next_corner_index = -1;
			{
				double highest_corner_dot = -INFINITY;
				for (int c = 0; c < corners.size(); c++) {
					if (corners[c] == convex_hull.back()) continue;
					
					v2 corner_direction = (corners[c] - convex_hull.back()).normalised_or_0();
					double corner_dot = dot(search_direction, corner_direction);
					
					if (corner_dot > highest_corner_dot) {
						highest_corner_dot = corner_dot;
						next_corner_index = c;
					}
				}
			}
			
			// check if the new corner is the same as the first corner.
			if (corners[next_corner_index] == convex_hull[0]) break; // the hull is complete.
			else {
				// add the corner to the shape and update the search direction
				search_direction = (corners[next_corner_index] - convex_hull.back()).normalised_or_0();
				convex_hull.push_back(corners[next_corner_index]);
			}
		}
		
		return corners.size() == convex_hull.size();
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
		
		if (!is_convex(corners)) {
			return false;
		}
		
		if (shape_out) {
			*shape_out = shape;
			shape_out->corners = corners;
			// set radius
			for (auto &corner: shape.corners) {
				if (corner.length() > shape.radius) shape.radius = corner.length();
			}
			shape_out->pos = ORIGIN;
			shape_out->angle = 0;
			shape_out->is_circle = false;
			return true;
		} else {
			return false;
		}
	}
	
	v2 get_minkowski_diffed_corner(Shape *shape, Shape *other_shape, v2 direction) {
		auto get_baked_corner_of_shape = [](Shape *shape, v2 corner_direction) {
			if (shape->is_circle) {
				return shape->pos + (corner_direction.normalised_or_0() * shape->radius);
			} else {
				v2 best_rotated_corner;
				double best_dot = -INFINITY;
				
				for (const auto &corner: shape->corners) {
					v2 rotated_corner = corner.rotated(shape->angle);
					double new_dot = dot(rotated_corner, corner_direction);
					
					if (new_dot > best_dot) {
						best_rotated_corner = rotated_corner;
						best_dot = new_dot;
					}
				}
				
				return shape->pos + best_rotated_corner;
			}
		};
		
		v2 baked_corner = get_baked_corner_of_shape(shape, direction);
		v2 other_baked_corner = get_baked_corner_of_shape(other_shape, -direction);
		
		return baked_corner - other_baked_corner;
	}
	
	bool origin_is_between_points(v2 a, v2 b) {
		v2 ao = ORIGIN - a;
		v2 bo = ORIGIN - b;
		v2 ab = b - a;
		v2 ba = a - b;
		return (dot(ao, ab) >= 0) && (dot(bo, ba) >= 0);
	}
	
	bool improve_2_simplex(vector<v2> &simplex, v2 &search_direction) {
		/*
		Find which simplex component the origin is closest
		to, or whether it is on the simplex line itself.
		*/
		if (origin_is_between_points(simplex[0], simplex[1])) {
			v2 line_normal = (simplex[1] - simplex[0]).right_normal_or_0();
			double origin_distance_from_line = dot(line_normal, ORIGIN - simplex[0]);
			
			if (fabs(origin_distance_from_line) <= LINE_THICKNESS) {
				return true; // The simplex contains the origin.
			} else {
				// The simplex is correct. Search on the side of the 2-simplex that contains the origin.
				search_direction = (simplex[1] - simplex[0]).normal_in_direction_or_0(ORIGIN - simplex[0]);
			}
		} else if (dot(simplex[1] - simplex[0], ORIGIN - simplex[0]) <= 0) {
			simplex = {simplex[0]}; // The origin is closest to point 0.
			search_direction = (ORIGIN - simplex[0]).normalised_or_0();
		} else {
			assert(dot(simplex[0] - simplex[1], ORIGIN - simplex[1]) <= 0);
			simplex = {simplex[1]}; // The origin is closest to point 1.
			search_direction = (ORIGIN - simplex[1]).normalised_or_0();
		}
		
		return false;
	}
	
	// returns true when the simplex contains the origin.
	bool improve_simplex(vector<v2> &simplex, v2 &search_direction) {
		if (simplex.size() == 3) {
			// cache some basic vectors
			v2 ab = simplex[1] - simplex[0];
			v2 bc = simplex[2] - simplex[1];
			v2 ca = simplex[0] - simplex[2];
			
			v2 ab_normal_away_from_c = ab.normal_in_direction_or_0(ca);
			v2 bc_normal_away_from_a = bc.normal_in_direction_or_0(ab);
			v2 ca_normal_away_from_b = ca.normal_in_direction_or_0(bc);
			
			// find which side of the triangle the origin is on, or if it's inside it.
			if (dot(ab_normal_away_from_c, ORIGIN - simplex[0]) > 0) {
				simplex = {simplex[0], simplex[1]};
			} else if (dot(bc_normal_away_from_a, ORIGIN - simplex[1]) > 0) {
				simplex = {simplex[1], simplex[2]};
			} else if (dot(ca_normal_away_from_b, ORIGIN - simplex[2]) > 0) {
				simplex = {simplex[2], simplex[0]};
			} else {
				return true; // the origin is inside the simplex.
			}
		}
		
		assert(simplex.size() == 2);
		return improve_2_simplex(simplex, search_direction);
	}
	
	bool shapes_are_overlapping(
		Shape *shape_a, Shape *shape_b,
		vector<v2> *simplex_out = nullptr // This is only used internally.
		) {
		
		// setting the initial direction like this maximises the
		// chance of the simplex covering the origin early.
		v2 search_direction = (shape_b->pos - shape_a->pos).right_normal_or_0();
		if (search_direction.is_0()) search_direction = v2(1, 0);
		
		vector<v2> simplex = { get_minkowski_diffed_corner(shape_a, shape_b, search_direction) };
		search_direction = ORIGIN - simplex[0]; // search toward the origin
		
		while (true) {
			simplex.push_back(get_minkowski_diffed_corner(shape_a, shape_b, search_direction));
			
			if (dot(simplex.back() - ORIGIN, search_direction) <= LINE_THICKNESS) {
				return false;
			}
			
			if (improve_simplex(simplex, search_direction)) {
				if (simplex_out != nullptr) *simplex_out = simplex;
				return true;
			}
		}
	}
	
	// Returns the amount that a is overlapping b.
	// Negating this amount from a->pos will resolve the overlap.
	v2 get_overlap_amount(Shape *shape_a, Shape *shape_b) {
		vector<v2> simplex;
		
		if (!shapes_are_overlapping(shape_a, shape_b, &simplex)) {
			return v2(0, 0);
		}
		
		if (simplex.size() < 3) {
			v2 pos_vector = (shape_b->pos - shape_a->pos).normalised_or_0();
			
			if (pos_vector.is_0()) {
				pos_vector.x = 1;
			}
			
			return pos_vector * LINE_THICKNESS;
		}
		
		int overlap_line_index = -1;
		
		while (true) {
			const double CORNER_SIMILARITY_TOLERANCE = LINE_THICKNESS;
			
			// get simplex line closest to origin
			double closest_line_distance = INFINITY;
			int closest_line_index = -1;
			for (int s0 = 0; s0 < simplex.size(); s0++) {
				int s1 = (s0+1) % simplex.size();
				v2 simplex_line_normal = (simplex[s1] - simplex[s0]).right_normal_or_0();
				
				double line_distance = fabs(dot(simplex_line_normal, ORIGIN - simplex[s0]));
				
				if (line_distance < closest_line_distance) {
					closest_line_index = s0;
					closest_line_distance = line_distance;
				}
			}
			
			// get the outer normal of that line and get the minkowski diffed corner in that direction
			v2 new_corner;
			int s0 = closest_line_index;
			int s1 = (s0+1) % simplex.size();
			v2 outer_normal = (simplex[s1] - simplex[s0]).normal_in_direction_or_0(simplex[s0] - ORIGIN);
			assert(!outer_normal.is_0());
			new_corner = get_minkowski_diffed_corner(shape_a, shape_b, outer_normal);
			
			// if the new corner is almost identical to one of the points that made the simplex,
			// break and handle it outside of the loop.
			{
				bool found_match = false;
				for (auto &simplex_corner : simplex) {
					if (simplex_corner.distance(new_corner) <= CORNER_SIMILARITY_TOLERANCE) {
						found_match = true;
						break;
					}
				}
				
				if (found_match) {
					overlap_line_index = s0;
					break;
				}
			}
			
			// else add the new corner to the simplex, turning the existing line into two.
			simplex.insert(simplex.begin()+s1, new_corner);
			assert(!contains_duplicates(simplex));
		}
		
		assert(overlap_line_index >= 0);
		
		// find the point on the line that is closest to the origin
		v2 overlap_line = simplex[overlap_line_index+1] - simplex[overlap_line_index];
		v2 overlap_line_unit = overlap_line.normalised_or_0();
		double len = dot(overlap_line_unit, ORIGIN - simplex[overlap_line_index]);
		v2 point_of_overlap = simplex[overlap_line_index] + overlap_line_unit * len;
		
		// the difference between the origin and that point is the overlap amount.
		v2 overlap_vector = point_of_overlap - ORIGIN;
		v2 overlap_direction_unit = overlap_vector.normalised_or_0();
		return overlap_direction_unit * (overlap_vector.length() + LINE_THICKNESS);
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
