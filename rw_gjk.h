/*
rw_gjk

A 2D collision detection and resolution library.
Based on the Gilbert–Johnson–Keerthi distance algorithm (GJK).
Written by Robin Wragg.

See end of file for license information.
*/

#ifndef RW_GJK_H
#define RW_GJK_H

#include <vector>

namespace rw_gjk {
	using namespace std;
	
	typedef double RW_FLOAT; // NOTE: ablility to use single-precision floats has not been proven.
	
	struct v2 {
		RW_FLOAT x, y;
		
		v2();
		v2(RW_FLOAT x_, RW_FLOAT y_);
		
		RW_FLOAT length() const;
		RW_FLOAT distance(const v2 &) const;
		bool isZero() const;
		v2 normalisedOrZero() const;
		v2 rightNormalOrZero() const;
		v2 normalInDirectionOrZero(v2 direction) const;
		v2 rotated(RW_FLOAT radians) const;
		
		bool operator==(const v2 &) const;
		v2 operator+(const v2 &) const;
		v2 operator-(const v2 &) const;
		v2 operator-() const;
		v2 operator*(const RW_FLOAT &) const;
		v2 operator/(const RW_FLOAT &) const;
	};
	
	struct Shape {
		v2 pos;
		float angle;
		
		vector<v2> corners; // TODO: Make private?
	};
	
	struct OverlapInfo {
		bool overlapping;
		v2 amount; // TODO: rename?
	};
	
	Shape *allocShape(vector<v2> corners);
	
	OverlapInfo getOverlapInfo(Shape *shapeA, Shape *shapeB);
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