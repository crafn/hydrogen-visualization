#ifndef QM_UTIL_HPP
#define QM_UTIL_HPP

#include <cassert>
#include <cmath>
#include <cstdarg>

namespace qm {

const double radToDeg= 57.2957795;
const double tau= 6.28318530;
const double pi= tau/2.0;

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

template <int Size>
struct StackString {
	char str[Size];
	int length;

	StackString(): length(0)
	{ str[0]= 0; }
	void append(const char *format, ...)
	{
		va_list args;
		va_start(args, format);

		std::size_t added= std::vsnprintf(str + length, Size - length, format, args);
		assert(added >= 0);

		va_end(args);

		length += added;
		assert(length < Size);
	}
};

template <typename T>
T clamp(T v, T min, T max)
{ return v < min ? min : v > max ? max : v; }

template <typename T>
T round(T v, int decimals= 0)
{
	double mul= std::pow(10, decimals);
	return std::floor(v*mul + 0.5)/mul;
}

inline
double fact(double n)
{ return n > 0 ? n*fact(n - 1) : 1; }

double binomial(double n, int r)
{
	double result= 1;
	for (int i= 1; i <= r; ++i)
		result *= n + 1 - i;
	for (int i= 1; i <= r; ++i)
		result /= i;
	return result;
}

/// Generalized Laguerre Polynomials
/// @param coeff should be size of n + 1, as coeff[n] will contain nth power
inline
void laguerre(double* coeff, int n, int alpha)
{
	int sign= 1;
	for (int i= 0; i <= n; ++i) {
		coeff[i]= 1.0*sign*binomial(n + alpha, n - i)/fact(i);
		sign *= -1;
	}
}

/// Legendre polynomials
/// @param coeff should be size of n + 1, as coeff[n] will contain nth power
inline
void legendre(double* coeff, int n)
{
	double mul= std::pow(2, n);
	for (int i= 0; i <= n; ++i)
		coeff[i]= mul*binomial(n, i)*binomial((i + n - 1)/2.0, n);
}

inline
void differentiate(double* coeff, int coeff_size, int diff_count)
{
	for (int i= 0; i < coeff_size; ++i) {
		if (diff_count > i)
			coeff[i]= 0;
		coeff[i - diff_count] += fact(diff_count)/fact(i - diff_count)*coeff[i];
	}
}

/// Coefficients for cosines in spherical harmonics
/// @param cos_coeff should be size of l + 1
inline
void sphericalHarmonics(double* cos_coeff, int l, int m)
{
	// Calculate coefficients for associated legendre polynomials
	// P_lm(cos(theta)) = (-1)^m * sin(theta)^m * D^m P_l(cos(theta))
	legendre(cos_coeff, l);
	differentiate(cos_coeff, l + 1, m);
	for (int i= 0; i <= l; ++i)
		cos_coeff[i] *= m % 2 ? -1 : 1;

	/// @todo Normalization factor
}

inline
void testMath()
{
	assert(fact(0) == 1);
	assert(fact(1) == 1);
	assert(fact(4) == 24);

	assert(binomial(10, 4) == 210);
	assert(binomial(1, 77) == 0);

	double lag[4]= {};
	double lag_correct[4]= {10, -10, 2.5, -1.0/6};
	laguerre(lag, 3, 2);
	for (int i= 0; i < 4; ++i)
		assert(std::abs(lag[i] - lag_correct[i]) < 0.0001);

}

} // qm

#endif // QM_UTIL_HPP
