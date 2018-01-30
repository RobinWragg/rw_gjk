struct v2 {
	double x, y;
	
	v2();
	v2(double x_, double y_);
	
	double length() const;
	double distance(const v2 &) const;
	bool is_zero() const;
	v2 normalised_or_zero() const;
	v2 right_normal_or_zero() const;
	v2 normal_in_direction_or_zero(v2 direction) const;
	v2 rotated(double radians) const;
	
	bool operator==(const v2 &) const;
	bool operator!=(const v2 &rh) const;
	v2 operator+(const v2 &) const;
	v2 operator-(const v2 &) const;
	v2 operator-() const;
	v2 operator*(const double &) const;
	v2 operator/(const double &) const;
};

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

v2 v2::normal_in_direction_or_zero(v2 direction) const {
	v2 normal_a = right_normal_or_zero(); // TODO: I don't think this needs to be normalised until we return.
	double dot_result = dot(normal_a, direction);
	
	if (dot_result > 0) {
		return normal_a;
	} else if (dot_result < 0) {
		v2 normal_b = normal_a * -1;
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

bool v2::operator!=(const v2 &rh) const {
	return !(*this == rh);
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