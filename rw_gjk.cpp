/*
rw_gjk

A 2D collision detection and resolution library.
Based on the Gilbert–Johnson–Keerthi distance algorithm (GJK).
Written by Robin Wragg.

See end of file for license information.
*/

// TODO: check for NaN
// TODO: look at better ways to iterate over vectors.
// TODO: remove asserts and handle those states in a graceful way
// TODO: write tests for allocShape():
	// will return nullptr if given less than three corners
	// will return nullptr if the collection of corners contains duplicates
	// will return nullptr if the corners form a concave shape
	// will return nullptr if any three of the corners form a straight line

#include "rw_gjk.h"
#include <cmath>
#include <cassert>
#include <string>

namespace rw_gjk {
	using namespace std;
	
	typedef float f32;
	typedef double f64;
	
	const f64 TINY_NUMBER = 0.0000001; // TODO: calculate how tiny this should be in a smarter way. this is basically half the thickness of a line.
	
	f64 dot(const v2 &a, const v2 &b) {
		return a.x*b.x + a.y*b.y;
	}
	
	v2::v2() {}
	
	v2::v2(f64 x_, f64 y_) {
		x = x_;
		y = y_;
	}
	
	f64 v2::length() const {
		return hypot(x, y);
	}
	
	f64 v2::distance(const v2 &rh) const {
		return hypot(x - rh.x, y - rh.y);
	}
	
	bool v2::isZero() const {
		return x == 0 && y == 0;
	}
	
	v2 v2::normalisedOrZero() const {
		if (isZero()) {
			return v2(0, 0);
		}
		
		f64 l = length();
		assert(l > 0);
		return *this / l;
	}
	
	v2 v2::rightNormalOrZero() const {
		return v2(y, -x).normalisedOrZero();
	}
	
	v2 v2::normalInDirectionOrZero(v2 direction) const {
		v2 normal0 = rightNormalOrZero(); // TODO: I don't think this needs to be normalised until we return.
		f64 dotResult = dot(normal0, direction);
		
		if (dotResult > 0) {
			assert(fabs(dot(normal0, *this)) < 0.00000001);
			return normal0;
		} else if (dotResult < 0) {
			v2 normal1 = normal0 * -1;
			assert(fabs(dot(normal1, *this)) < 0.00000001);
			return normal1;
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
	
	v2 v2::operator*(const f64 &rh) const {
		return v2(x * rh, y * rh);
	}
	
	v2 v2::operator/(const f64 &rh) const {
		return v2(x / rh, y / rh);
	}
	
	v2 origin = v2(0, 0);
	string errorDescription = "";
	
	bool containsDuplicates(vector<v2> vertices) {
		for (int s0 = 0; s0 < vertices.size(); s0++) {
			for (int s1 = s0+1; s1 < vertices.size(); s1++) {
				if (vertices[s0] == vertices[s1]) return true;
			}
		}
		return false;
	}
	
	Shape *allocShape(vector<v2> corners) {
		if (corners.size() < 3) {
			errorDescription = "Less than 3 corners were given to allocShape().";
			return nullptr;
		}
		
		if (containsDuplicates(corners)) {
			return nullptr;
		}
		
		Shape *shape = new Shape;
		
		// check if the corners form a convex shape
		{
			// add the leftmost corner to the shape
			v2 leftmostCorner = corners.front();
			for (auto &corner: corners) {
				if (corner.x < leftmostCorner.x) leftmostCorner = corner;
			}
			shape->corners.push_back(leftmostCorner);
			
			// build a convex shape out of the corners, in any order.
			v2 searchDir = v2(0, 1);
			while (true) {
				// find the next corner
				f64 highestCornerDot = -INFINITY;
				v2 bestCorner;
				for (int c = 0; c < corners.size(); c++) {
					if (corners[c] == shape->corners.back()) continue;
					
					v2 cornerDifference = corners[c] - shape->corners.back();
					f64 cornerDot = dot(searchDir, cornerDifference.normalisedOrZero());
					
					if (cornerDot >= 0 && cornerDot > highestCornerDot) {
						highestCornerDot = cornerDot;
						bestCorner = corners[c];
					}
				}
				
				// if a valid corner was found
				if (highestCornerDot >= 0) {
					// check if it's the same as the first corner
					if (bestCorner == shape->corners[0]) break; // the shape is complete
					else {
						// add the corner to the shape and update the search direction
						searchDir = (bestCorner - shape->corners.back()).normalisedOrZero();
						shape->corners.push_back(bestCorner);
					}
				} else {
					// no valid corner was found in the current search direction, so rotate it 90 degrees
					searchDir = searchDir.rightNormalOrZero();
				}
			}
			
			// compare the shape against the original corners.
			assert(shape->corners.size() <= corners.size());
			
			if (shape->corners.size() < corners.size()) {
				delete shape;
				return nullptr; // at least one of the corners forms a concave angle.
			}
		}
		
		shape->pos = origin;
		shape->angle = 0;
		
		return shape;
	}
	
	v2 getMinkowskiDiffedEdge(Shape *shape, Shape *otherShape, v2 direction) {
		auto getEdgeOfRotatedShape = [](Shape *shape, v2 direction) {
			v2 bestCrnr;
			f64 bestDot = -INFINITY;
			
			for (const auto &corner: shape->corners) {
				v2 rotatedCrnr = corner.rotated(shape->angle);
				f64 newDot = dot(rotatedCrnr, direction);
				
				if (newDot > bestDot) {
					bestCrnr = rotatedCrnr;
					bestDot = newDot;
				}
			}
			
			return bestCrnr;
		};
		
		v2 shapeWorldEdge = shape->pos + getEdgeOfRotatedShape(shape, direction);
		v2 otherShapeWorldEdge = otherShape->pos + getEdgeOfRotatedShape(otherShape, -direction);
		return shapeWorldEdge - otherShapeWorldEdge;
	}
	
	bool originIsBetweenPoints(v2 a, v2 b) {
		v2 ao = origin - a;
		v2 bo = origin - b;
		v2 ab = b - a;
		v2 ba = a - b;
		return (dot(ao, ab) > 0) && (dot(bo, ba) > 0);
	}
	
	bool improve2Simplex(vector<v2> &simplex, v2 &searchDir) {
		assert(simplex.size() == 2);
		assert(!containsDuplicates(simplex));
		
		if (originIsBetweenPoints(simplex[0], simplex[1])) {
			// the simplex is correct.
			// search on the side of the 2-simplex that contains the origin.
			v2 simplexLine = simplex[1] - simplex[0];
			v2 directionToOrigin = origin - simplex[0];
			searchDir = simplexLine.normalInDirectionOrZero(directionToOrigin);
		} else if (dot(simplex[1] - simplex[0], origin - simplex[0]) <= 0) {
			// the closest part of the simplex is the point simplex[0], so update the simplex.
			simplex = {simplex[0]};
			searchDir = (origin - simplex[0]).normalisedOrZero();
		} else {
			// the closest part of the simplex is the point simplex[1], so update the simplex.
			assert(dot(simplex[0] - simplex[1], origin - simplex[1]) <= 0);
			simplex = {simplex[1]};
			searchDir = (origin - simplex[1]).normalisedOrZero();
		}
		
		if (searchDir.isZero()) {
			return true; // the origin is precisely on the simplex edge.
		}
		
		assert(simplex.size() == 1 || simplex.size() == 2);
		assert(searchDir.length() > 0.9999 && searchDir.length() < 1.0001);
		return false;
	}
	
	// returns true when the simplex contains the origin.
	bool improveSimplex(vector<v2> &simplex, v2 &searchDir) {
		if (simplex.size() == 3) {
			v2 abNormalAwayFromC = (simplex[1] - simplex[0]).normalInDirectionOrZero(simplex[0] - simplex[2]);
			assert(!abNormalAwayFromC.isZero());
			v2 bcNormalAwayFromA = (simplex[2] - simplex[1]).normalInDirectionOrZero(simplex[1] - simplex[0]);
			assert(!bcNormalAwayFromA.isZero());
			v2 caNormalAwayFromB = (simplex[0] - simplex[2]).normalInDirectionOrZero(simplex[2] - simplex[1]);
			assert(!caNormalAwayFromB.isZero());
			
			// find which side of the triangle the origin is on, or if it's inside it.
			if (dot(abNormalAwayFromC, origin - simplex[0]) > 0) {
				simplex = {simplex[0], simplex[1]};
				return improve2Simplex(simplex, searchDir);
			} else if (dot(bcNormalAwayFromA, origin - simplex[1]) > 0) {
				simplex = {simplex[1], simplex[2]};
				return improve2Simplex(simplex, searchDir);
			} else if (dot(caNormalAwayFromB, origin - simplex[2]) > 0) {
				simplex = {simplex[2], simplex[0]};
				return improve2Simplex(simplex, searchDir);
			} else {
				return true; // the origin is inside the simplex.
			}
		}
		
		assert(simplex.size() == 2);
		return improve2Simplex(simplex, searchDir);
	}
	
	bool shapesAreOverlapping(
		Shape *shapeA, Shape *shapeB,
		vector<v2> *simplexOut = nullptr // This is only used internally.
		) {
		
		v2 searchDirection = v2(0, 1); // this can be anything except zero
		vector<v2> simplex = {getMinkowskiDiffedEdge(shapeA, shapeB, searchDirection)};
		
		searchDirection = origin - simplex[0]; // search toward the origin
		
		while (true) { // TODO: remove this infinite loop.
			assert(simplex.size() < 3);
			simplex.push_back(getMinkowskiDiffedEdge(shapeA, shapeB, searchDirection));
			
			f64 dotResult = dot(simplex.back(), searchDirection);
			if (dotResult <= TINY_NUMBER) {
				return false;
			}
			
			assert(!containsDuplicates(simplex));
			
			if (improveSimplex(simplex, searchDirection)) {
				if (simplexOut != nullptr) *simplexOut = simplex;
				return true;	
			}
		}
	}
	
	struct OverlapInfo {
		bool overlapping;
		v2 amount; // TODO: rename?
		bool failure;
	};
	
	vector<vector<v2>> sss; // TODO: this is temporary, for debugging issues in getOverlapInfo().
	vector<v2> mink;
	OverlapInfo getOverlapInfo(Shape *shapeA, Shape *shapeB) {
		
		auto getAmountForOriginOnEdge = [shapeA, shapeB]() {
			v2 posDifference = shapeB->pos - shapeA->pos;
			
			if (posDifference.isZero()) {
				posDifference.x = TINY_NUMBER;
			}
			
			return posDifference.normalisedOrZero() * TINY_NUMBER;
		};
		
		OverlapInfo info;
		vector<v2> simplex;
		
		info.overlapping = shapesAreOverlapping(shapeA, shapeB, &simplex);
		info.failure = false;
		sss = {simplex};
		
		if (info.overlapping) {
			if (simplex.size() < 2) {
				info.amount = getAmountForOriginOnEdge();
				return info;
			}
			
			while (true) {
				// find simplex line closest to origin
				int closestEdgeIndex = 0;
				{
					f64 closestDistance = INFINITY;
					
					for (int s0 = 0; s0 < simplex.size(); s0++) {
						int s1 = (s0 + 1) % simplex.size();
						int s2 = (s0 + 2) % simplex.size();
						
						v2 simplexLine = simplex[s1] - simplex[s0];
						v2 approxOuterNormalDir = simplex[s0] - simplex[s2];
						v2 outerNormal = simplexLine.normalInDirectionOrZero(approxOuterNormalDir);
						
						if (outerNormal.isZero()) {
							info.amount = getAmountForOriginOnEdge();
							return info;
						}
						
						f64 distance = dot(outerNormal, simplex[s0]);
						if (distance < closestDistance) {
							closestDistance = distance;
							closestEdgeIndex = s0;
						}
					}
				}
				
				// get outer normal of the line and get minkowski diffed edge in that direction
				v2 newSimplexCrnr;
				int s0 = closestEdgeIndex;
				int s1 = (closestEdgeIndex + 1) % simplex.size();
				int s2 = (closestEdgeIndex + 2) % simplex.size();
				{
					v2 outerNormal = (simplex[s1] - simplex[s0]).normalInDirectionOrZero(simplex[s0] - simplex[s2]);
					assert(!outerNormal.isZero());
					newSimplexCrnr = getMinkowskiDiffedEdge(shapeA, shapeB, outerNormal);
				}
				
				// if the new edge is almost identical to one of the points that made the simplex line,
				if (newSimplexCrnr.distance(simplex[s0]) < TINY_NUMBER || newSimplexCrnr.distance(simplex[s1]) < TINY_NUMBER) {
					// find the point on the line that is closest to the origin
					v2 outerNormal = (simplex[s1] - simplex[s0]).normalInDirectionOrZero(simplex[s0] - simplex[s2]);
					assert(!outerNormal.isZero());
					f64 lengthOfPoint = dot(outerNormal, simplex[s0]);
					
					// the difference between the origin and that point is the overlap amount.
					info.amount = outerNormal.normalisedOrZero() * (lengthOfPoint + TINY_NUMBER);
				
					// break out of the loop
					break;
				} else {
					// add the new corner to the simplex, turning the existing line into two.
					simplex.insert(simplex.begin()+s1, newSimplexCrnr);
					assert(!containsDuplicates(simplex));
					sss.push_back(simplex);
				}
			}
		}
		
		return info;
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