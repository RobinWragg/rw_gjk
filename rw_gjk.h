/*
rw_gjk

A 2D collision detection and resolution library.
Based on the Gilbert–Johnson–Keerthi distance algorithm (GJK).
Written by Robin Wragg.

See end of file for license information.
*/

// TODO: check for NaN?
// TODO: remove occurrences of while (true).
// TODO: remove asserts and handle those states in a graceful way.
// TODO: add caching of rotations.
// TODO: add error strings.
// TODO: calculate TINY_NUMBER appropriately. It has something to do with largest_radius in shapes_are_overlapping(). TINY_NUMBER is basically used as half the thickness of a line. It should be very small but never small enough to cause IEEE-float-related problems.

#ifndef RW_GJK_HEADER
#define RW_GJK_HEADER

#include <vector>

namespace rw_gjk {
	using namespace std;
	
	struct v2 {
		double x, y;
		
		v2();
		v2(double x_, double y_);
		
		double length() const;
		double distance(const v2 &) const;
		bool is_zero() const;
		v2 normalised_or_zero() const;
		v2 right_normal_or_zero() const;
		v2 normal_in_direction_or_0(v2 direction) const;
		v2 rotated(double radians) const;
		
		bool operator==(const v2 &) const;
		v2 operator+(const v2 &) const;
		v2 operator-(const v2 &) const;
		v2 operator-() const;
		v2 operator*(const double &) const;
		v2 operator/(const double &) const;
	};
	
	struct Shape {
		v2 pos;
		float angle;
		double radius;
		
		vector<v2> corners;
		
		float cached_angle;
		vector<v2> cached_corners;
	};
	
	bool try_make_shape(vector<v2> corners, Shape *shape_out);
	
	v2 get_overlap_amount(Shape *a, Shape *b);
}

#endif

#ifdef RW_GJK_IMPLEMENTATION

#include <cmath>
#include <cassert>
#include <string>

namespace rw_gjk {
	const double TINY_NUMBER = 0.0000001;
	
	double dot(const v2 &a, const v2 &b) {
		return a.x*b.x + a.y*b.y;
	}
	
	v2::v2() {
		// These are NANs at the moment to catch uninitialised-variable bugs
		x = NAN;
		y = NAN;
	}
	
	v2::v2(double x_, double y_) {
		x = x_;
		y = y_;
	}
	
	double v2::length() const {
		return hypot(x, y);
	}
	
	double v2::distance(const v2 &rh) const {
		return hypot(x - rh.x, y - rh.y);
	}
	
	bool v2::is_zero() const {
		return x == 0 && y == 0;
	}
	
	v2 v2::normalised_or_zero() const {
		if (is_zero()) {
			return v2(0, 0);
		}
		
		double l = length();
		assert(l > 0);
		return *this / l;
	}
	
	v2 v2::right_normal_or_zero() const {
		return v2(y, -x).normalised_or_zero();
	}
	
	v2 v2::normal_in_direction_or_0(v2 direction) const {
		v2 normal_a = right_normal_or_zero(); // TODO: I don't think this needs to be normalised until we return.
		double dot_result = dot(normal_a, direction);
		
		if (dot_result > 0) {
			assert(fabs(dot(normal_a, *this)) < TINY_NUMBER);
			return normal_a;
		} else if (dot_result < 0) {
			v2 normal_b = normal_a * -1;
			assert(fabs(dot(normal_b, *this)) < TINY_NUMBER);
			return normal_b;
		} else {
			return v2(0, 0);
		}
	}
	
	v2 v2::rotated(double radians) const {
		v2 oldv = *this;
		radians *= -1; // flip the sign so that a positive number rotates the vector clockwise
		return v2(oldv.x*cos(radians) - oldv.y*sin(radians),
			oldv.x*sin(radians) + oldv.y*cos(radians));
	}
	
	bool v2::operator==(const v2 &rh) const {
		return x == rh.x && y == rh.y;
	}
	
	v2 v2::operator+(const v2 &rh) const {
		return v2(x + rh.x, y + rh.y);
	}
	
	v2 v2::operator-(const v2 &rh) const {
		return v2(x - rh.x, y - rh.y);
	}
	
	v2 v2::operator-() const {
		return v2(-x, -y);
	}
	
	v2 v2::operator*(const double &rh) const {
		return v2(x * rh, y * rh);
	}
	
	v2 v2::operator/(const double &rh) const {
		return v2(x / rh, y / rh);
	}
	
	v2 origin = v2(0, 0);
	
	bool contains_duplicates(vector<v2> vertices) {
		for (int s0 = 0; s0 < vertices.size(); s0++) {
			for (int s1 = s0+1; s1 < vertices.size(); s1++) {
				if (vertices[s0] == vertices[s1]) return true;
			}
		}
		return false;
	}
	
	bool try_make_shape(vector<v2> corners, Shape *shape_out) {
		if (corners.size() < 3) {
			return false;
		}
		
		if (contains_duplicates(corners)) {
			return false;
		}
		
		Shape shape;
		
		// check if the corners form a convex shape
		{
			// add the leftmost corner to the shape
			v2 leftmost_corner = corners.front();
			for (auto &corner: corners) {
				if (corner.x < leftmost_corner.x) leftmost_corner = corner;
			}
			shape.corners.push_back(leftmost_corner);
			
			// build a convex shape out of the corners, in any order.
			v2 search_dir = v2(0, 1);
			while (true) {
				// find the next corner
				double highest_corner_dot = -INFINITY;
				v2 best_corner;
				for (auto &corner: corners) {
					if (corner == shape.corners.back()) continue;
					
					v2 corner_difference = corner - shape.corners.back();
					double corner_dot = dot(search_dir, corner_difference.normalised_or_zero());
					
					if (corner_dot >= 0 && corner_dot > highest_corner_dot) {
						highest_corner_dot = corner_dot;
						best_corner = corner;
					}
				}
				
				// if a valid corner was found
				if (highest_corner_dot >= 0) {
					// check if it's the same as the first corner
					if (best_corner == shape.corners[0]) break; // the shape is complete
					else {
						// add the corner to the shape and update the search direction
						search_dir = (best_corner - shape.corners.back()).normalised_or_zero();
						shape.corners.push_back(best_corner);
					}
				} else {
					// no valid corner was found in the current search direction, so rotate it 90 degrees
					search_dir = search_dir.right_normal_or_zero();
				}
			}
			
			// compare the shape against the original corners.
			assert(shape.corners.size() <= corners.size());
			
			if (shape.corners.size() < corners.size()) {
				return false; // at least one of the corners forms a concave angle.
			}
		}
		
		// set radius
		for (auto &corner: shape.corners) {
			if (corner.length() > shape.radius) shape.radius = corner.length();
		}
		
		*shape_out = shape;
		shape_out->pos = origin;
		shape_out->angle = 0;
		shape_out->cached_angle = NAN;
		shape_out->cached_corners.resize(shape_out->corners.size());
		
		return true;
	}
	
	void update_shape_cache(Shape *shape) {
		if (shape->angle == shape->cached_angle) return;
		
		for (int i = 0; i < shape->corners.size(); i++) {
			shape->cached_corners[i] = shape->corners[i].rotated(shape->angle);
		}
		
		shape->cached_angle = shape->angle;
	}
	
	v2 get_minkowski_diffed_edge(Shape *shape, Shape *other_shape, v2 direction) {
		auto getEdgeOfRotatedShape = [](Shape *shape, v2 direction) {
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
		};
		
		v2 shape_world_edge = shape->pos + getEdgeOfRotatedShape(shape, direction);
		v2 other_shape_world_edge = other_shape->pos + getEdgeOfRotatedShape(other_shape, -direction);
		return shape_world_edge - other_shape_world_edge;
	}
	
	bool origin_is_between_points(v2 a, v2 b) {
		v2 ao = origin - a;
		v2 bo = origin - b;
		v2 ab = b - a;
		v2 ba = a - b;
		return (dot(ao, ab) > 0) && (dot(bo, ba) > 0);
	}
	
	bool improve_2_simplex(vector<v2> &simplex, v2 &search_dir) {
		assert(simplex.size() == 2);
		assert(!contains_duplicates(simplex));
		
		if (origin_is_between_points(simplex[0], simplex[1])) {
			// the simplex is correct.
			// search on the side of the 2-simplex that contains the origin.
			v2 simplex_line = simplex[1] - simplex[0];
			v2 direction_to_origin = origin - simplex[0];
			search_dir = simplex_line.normal_in_direction_or_0(direction_to_origin);
		} else if (dot(simplex[1] - simplex[0], origin - simplex[0]) <= 0) {
			// the closest part of the simplex is the point simplex[0], so update the simplex.
			simplex = {simplex[0]};
			search_dir = (origin - simplex[0]).normalised_or_zero();
		} else {
			// the closest part of the simplex is the point simplex[1], so update the simplex.
			assert(dot(simplex[0] - simplex[1], origin - simplex[1]) <= 0);
			simplex = {simplex[1]};
			search_dir = (origin - simplex[1]).normalised_or_zero();
		}
		
		if (search_dir.is_zero()) {
			return true; // the origin is precisely on the simplex edge.
		}
		
		assert(simplex.size() == 1 || simplex.size() == 2);
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
			
			v2 ab_normal_away_from_c = ab.normal_in_direction_or_0(ca);
			v2 bc_normal_away_from_a = bc.normal_in_direction_or_0(ab);
			v2 ca_normal_away_from_b = ca.normal_in_direction_or_0(bc);
			
			// find which side of the triangle the origin is on, or if it's inside it.
			if (dot(ab_normal_away_from_c, origin - simplex[0]) > 0) {
				simplex = {simplex[0], simplex[1]};
				return improve_2_simplex(simplex, search_dir);
			} else if (dot(bc_normal_away_from_a, origin - simplex[1]) > 0) {
				simplex = {simplex[1], simplex[2]};
				return improve_2_simplex(simplex, search_dir);
			} else if (dot(ca_normal_away_from_b, origin - simplex[2]) > 0) {
				simplex = {simplex[2], simplex[0]};
				return improve_2_simplex(simplex, search_dir);
			} else {
				return true; // the origin is inside the simplex.
			}
		}
		
		assert(simplex.size() == 2);
		return improve_2_simplex(simplex, search_dir);
	}
	
	bool shapes_are_close(Shape *a, Shape *b) {
		return a->pos.x+a->radius >= b->pos.x-b->radius
			&& a->pos.x-a->radius <= b->pos.x+b->radius
			&& a->pos.y+a->radius >= b->pos.y-b->radius
			&& a->pos.y-a->radius <= b->pos.y+b->radius;
	}
	
	bool shapes_are_overlapping(
		Shape *shape_a, Shape *shape_b,
		vector<v2> *simplexOut = nullptr // This is only used internally.
		) {
		
		// if (!shapes_are_close(shape_a, shape_b)) return false; // TODO: uncomment this.
		
		{ // Unfinished code for correctly setting TINY_NUMBER. See todo list at top of file.
			double largest_radius = max(shape_a->radius, shape_b->radius);
			
		}
		
		v2 search_direction = v2(0, 1); // this can be anything except zero
		vector<v2> simplex = {get_minkowski_diffed_edge(shape_a, shape_b, search_direction)};
		
		search_direction = origin - simplex[0]; // search toward the origin
		
		while (true) {
			assert(simplex.size() < 3);
			simplex.push_back(get_minkowski_diffed_edge(shape_a, shape_b, search_direction));
			
			double dot_result = dot(simplex.back(), search_direction);
			if (dot_result <= TINY_NUMBER) {
				return false;
			}
			
			assert(!contains_duplicates(simplex));
			
			if (improve_simplex(simplex, search_direction)) {
				if (simplexOut != nullptr) *simplexOut = simplex;
				return true;	
			}
		}
	}
	
	v2 get_overlap_amount(Shape *shape_a, Shape *shape_b) {
		
		auto get_amount_for_origin_on_edge = [shape_a, shape_b]() {
			v2 pos_difference = shape_b->pos - shape_a->pos;
			
			if (pos_difference.is_zero()) {
				pos_difference.x = TINY_NUMBER;
			}
			
			return pos_difference.normalised_or_zero() * TINY_NUMBER;
		};
		
		vector<v2> simplex;
		if (!shapes_are_overlapping(shape_a, shape_b, &simplex)) {
			return v2(0, 0);
		}
		
		if (simplex.size() < 2) {
			v2 amount = get_amount_for_origin_on_edge();
			assert(!amount.is_zero());
			return amount;
		}
		
		while (true) {
			// find simplex line closest to origin
			int closest_edge_index = 0;
			{
				double closest_distance = INFINITY;
				
				for (int s0 = 0; s0 < simplex.size(); s0++) {
					int s1 = (s0 + 1) % simplex.size();
					int s2 = (s0 + 2) % simplex.size();
					
					v2 simplex_line = simplex[s1] - simplex[s0];
					v2 approx_outer_normal_dir = simplex[s0] - simplex[s2];
					v2 outer_normal = simplex_line.normal_in_direction_or_0(approx_outer_normal_dir);
					
					if (outer_normal.is_zero()) {
						v2 amount = get_amount_for_origin_on_edge();
						assert(!amount.is_zero());
						return amount;
					}
					
					double distance = dot(outer_normal, simplex[s0]);
					if (distance < closest_distance) {
						closest_distance = distance;
						closest_edge_index = s0;
					}
				}
			}
			
			// get outer normal of the line and get minkowski diffed edge in that direction
			v2 new_simplex_corner;
			int s0 = closest_edge_index;
			int s1 = (closest_edge_index + 1) % simplex.size();
			int s2 = (closest_edge_index + 2) % simplex.size();
			{
				v2 outer_normal =
					(simplex[s1] - simplex[s0]).normal_in_direction_or_0(simplex[s0] - simplex[s2]);
				assert(!outer_normal.is_zero());
				new_simplex_corner = get_minkowski_diffed_edge(shape_a, shape_b, outer_normal);
			}
			
			// if the new edge is almost identical to one of the points that made the simplex line,
			if (new_simplex_corner.distance(simplex[s0]) < TINY_NUMBER
				|| new_simplex_corner.distance(simplex[s1]) < TINY_NUMBER) {
				// find the point on the line that is closest to the origin
				v2 outer_normal =
					(simplex[s1] - simplex[s0]).normal_in_direction_or_0(simplex[s0] - simplex[s2]);
				assert(!outer_normal.is_zero());
				double length_of_point = dot(outer_normal, simplex[s0]);
				
				// the difference between the origin and that point is the overlap amount.
				v2 amount = outer_normal.normalised_or_zero() * (length_of_point + TINY_NUMBER);
				assert(!amount.is_zero());
				return amount;
			} else {
				// add the new corner to the simplex, turning the existing line into two.
				simplex.insert(simplex.begin()+s1, new_simplex_corner);
				assert(!contains_duplicates(simplex));
			}
		}
		
		// Control should never reach here, unless we break out of the 'while (true)'' loop to prevent
		// an infinite loop, at which point we should write an error description.
		// see the todo list.
		assert(false);
		return v2(NAN, NAN);
	}
}

#endif

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