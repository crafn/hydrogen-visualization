#ifndef QM_UTIL_HPP
#define QM_UTIL_HPP

#include <cassert>
#include <cmath>
#include <cstdarg>

namespace qm {

const double radToDeg= 57.2957795;

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

	bool operator==(Vec2 other) const { return x == other.x && y == other.y; }
	bool operator!=(Vec2 other) const { return !(*this == other); }
};

typedef Vec2<float> Vec2f;
typedef Vec2<int> Vec2i;

template <typename T, typename U>
T cast(Vec2<U> v) { return T((typename T::Value)v.x, (typename T::Value)v.y); }

inline
Vec2f fitToGrid(Vec2f v, Vec2i reso)
{
	Vec2f hreso= cast<Vec2f>(reso)*0.5f;
	Vec2i grid_v= cast<Vec2i>((v + Vec2f(1, 1))*hreso);
	return cast<Vec2f>(grid_v)/hreso- Vec2f(1, 1);
}

template <std::size_t Size>
struct StackString {
	char str[Size];
	std::size_t length;
};

template <std::size_t Size>
void append(StackString<Size>& s, const char* format, ...)
{
	va_list args;
	va_start(args, format);

	std::size_t added= std::vsnprintf(s.str + s.length, Size - s.length, format, args);
	assert(added >= 0);

	s.length += added;
	assert(s.length < Size);

	va_end(args);
}

struct String {
	char* str;
	std::size_t length;
};

inline
String createString()
{
	String s= {};
	s.str= (char*)std::calloc(1, 1);
	return s;
}

inline
void destroyString(String& s)
{
	std::free(s.str);
	s.str= NULL;
}

inline
void append(String& s, const char* format, ...)
{
	assert(s.str);

	va_list args;
	va_list args2;
	va_start(args, format);
	va_copy(args2, args);

	std::size_t new_len= s.length + std::vsnprintf(NULL, 0, format, args);
	s.str= (char*)std::realloc((void*)s.str, new_len + 1);
	assert(s.str);
	std::vsnprintf(s.str + s.length, new_len + 1, format, args2);
	s.length= new_len;

	va_end(args);
}

template <typename T, std::size_t Size>
struct StackArray {
	T data[Size];
	std::size_t size;
};

template <typename T, std::size_t Size>
void push(StackArray<T, Size>& a, T t)
{
	assert(a.size < Size);
	a.data[a.size]= t;
	++a.size;
	assert(a.size <= Size);
}

template <typename T, std::size_t Size>
void pop(StackArray<T, Size>& a)
{
	assert(a.size > 0);
	--a.size;
}

template <typename T, std::size_t Size>
T& last(StackArray<T, Size>& a)
{
	assert(a.size > 0);
	return a.data[a.size - 1];
}

} // qm

#endif // QM_UTIL_HPP
