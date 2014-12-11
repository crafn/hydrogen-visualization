#ifndef QM_UTIL_HPP
#define QM_UTIL_HPP

#include <cmath>

namespace qm {

template <typename T>
struct Vec2 {
	T x, y;

	Vec2(T x= 0, T y= 0): x(x), y(y) {}
	T lengthSqr() const { return x*x + y*y; }
	T length() const { return std::sqrt(lengthSqr()); }

	Vec2 operator*(T scalar) const { return Vec2(x*scalar, y*scalar); }
	Vec2 operator*(Vec2 other) const { return Vec2(x*other.x, y*other.y); }
	Vec2 operator+(Vec2 other) const { return Vec2(x+other.x, y+other.y); }
};

typedef Vec2<float> Vec2f;
typedef Vec2<int> Vec2i;

} // qm

#endif // QM_UTIL_HPP
