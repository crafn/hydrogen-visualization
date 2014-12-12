#ifndef QM_UTIL_HPP
#define QM_UTIL_HPP

#include <cmath>

namespace qm {

const float radToDeg= 57.2957795f;
const float tau= 6.28318530f;

template <typename T>
struct Vec2 {
	T x, y;

	typedef T Value;
	Vec2(): x(0), y(0) {}
	Vec2(T x, T y): x(x), y(y) {}

	T lengthSqr() const { return x*x + y*y; }
	T length() const { return std::sqrt(lengthSqr()); }

	Vec2 operator*(T scalar) const { return Vec2(x*scalar, y*scalar); }
	Vec2 operator*(Vec2 other) const { return Vec2(x*other.x, y*other.y); }
	Vec2 operator/(Vec2 other) const { return Vec2(x/other.x, y/other.y); }
	Vec2 operator+(Vec2 other) const { return Vec2(x+other.x, y+other.y); }
	Vec2 operator-(Vec2 other) const { return Vec2(x-other.x, y-other.y); }
	
	Vec2& operator*=(T scalar) { return *this= *this*scalar; }
	Vec2& operator*=(Vec2 other) { return *this= *this*other; }
	Vec2& operator/=(Vec2 other) { return *this= *this/other; }
	Vec2& operator+=(Vec2 other) { return *this= *this+other; }
	Vec2& operator-=(Vec2 other) { return *this= *this-other; }
};

typedef Vec2<float> Vec2f;
typedef Vec2<int> Vec2i;

template <typename T, typename U>
T casted(Vec2<U> v) { return T((typename T::Value)v.x, (typename T::Value)v.y); }

inline
Vec2f fitToGrid(Vec2f v, Vec2i reso)
{
	Vec2f hreso= casted<Vec2f>(reso)*0.5f;
	Vec2i grid_v= casted<Vec2i>((v + Vec2f(1, 1))*hreso);
	return casted<Vec2f>(grid_v)/hreso- Vec2f(1, 1);
}

inline
float clamp(float v, float min, float max)
{ return v < min ? min : v > max ? max : v; }

} // qm

#endif // QM_UTIL_HPP
