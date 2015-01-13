#ifndef QM_MATH_HPP
#define QM_MATH_HPP

#include <cassert>

namespace qm {

const double tau= 6.28318530;
const double pi= tau/2.0;

#define CLAMP(v, min, max) (v < (min) ? (min) : v > (max) ? (max) : v)

template <typename T>
T round(T v, int decimals= 0)
{
	double mul= std::pow(10, decimals);
	return std::floor(v*mul + 0.5)/mul;
}

inline
double fact(double n)
{ return n > 0 ? n*fact(n - 1) : 1; }

inline
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
		if (i + diff_count >= coeff_size)
			coeff[i]= 0.0;
		else
			coeff[i]= fact(i + diff_count)/fact(i)*coeff[i + diff_count];
	}
}

/// Coefficients for cosines in spherical harmonics
/// Embeds plus or minus from the front of Y to the coefficients depending on m
/// @param cos_coeff should be size of l + 1
inline
void sphericalHarmonics(double* cos_coeff, int l, int m)
{
	assert(l >= 0);

	int m_sign; 
	if (m >= 0)
		m_sign= m % 2 ? -1 : 1;
	else
		m_sign= 1;
	m= std::abs(m);

	// Calculate coefficients for associated legendre polynomials
	// P_lm(cos(theta)) = (-1)^m * sin(theta)^m * D^m P_l(cos(theta))
	legendre(cos_coeff, l);
	differentiate(cos_coeff, l + 1, m);
	double normalization= 
		std::sqrt( (2*l + 1.0)/(4*pi)*fact(l - m)/fact(l + m) );
	for (int i= 0; i <= l; ++i)
		cos_coeff[i] *= m_sign*normalization;
}

inline
void testMath()
{
	assert(fact(0) == 1);
	assert(fact(1) == 1);
	assert(fact(4) == 24);

	assert(binomial(10, 4) == 210);
	assert(binomial(1, 77) == 0);

	{ // Laguerre coefficients
		const int size= 4;
		double lag[size]= {};
		double lag_correct[size]= {10, -10, 2.5, -1.0/6};
		laguerre(lag, size - 1, 2);
		for (int i= 0; i < size; ++i) {
			assert(std::abs(lag[i] - lag_correct[i]) < 0.0001);
		}
	}

	{ // Differentiation
		const int size= 5;
		double diff[size]= {5, 4, 3, 2, 1};
		double diff_correct[size]= {6, 12, 12, 0, 0};
		differentiate(diff, size, 2);
		for (int i= 0; i < size; ++i) {
			assert(std::abs(diff[i] - diff_correct[i]) < 0.0001);
		}
	}

	{ // Legendre coefficients
		const int size= 5;
		double leg[size];
		double leg_correct[size]= {3.0/8, 0, -30.0/8, 0, 35.0/8};
		legendre(leg, size - 1);
		for (int i= 0; i < size; ++i) {
			assert(std::abs(leg[i] - leg_correct[i]) < 0.0001);
		}
	}

	{ // Spherical harmonic coefficients
		{
			double sphe[1];
			sphericalHarmonics(sphe, 0, 0);
			assert(std::abs(sphe[0] - 1.0/std::sqrt(4*pi)) < 0.001);
		}

		{
			const int size= 3;
			double sphe[size];
			double sphe_correct[size]= { 0, -std::sqrt(15.0/(8*pi)), 0};
			sphericalHarmonics(sphe, size - 1, 1);
			for (int i= 0; i < size; ++i) {
				assert(std::abs(sphe[i] - sphe_correct[i]) < 0.0001);
			}
		}

		{
			const int size= 10;
			double sphe[size];
			double mul= 1.0/128*std::sqrt(40755.0/pi);
			double sphe_correct[size]= { 0, -3*mul, 0, 17*mul, 0, 0, 0, 0, 0, 0};
			sphericalHarmonics(sphe, size - 1, 6);
			for (int i= 0; i < size; ++i) {
				assert(std::abs(sphe[i] - sphe_correct[i]) < 0.0001);
			}
		}
	}
}

} // qm

#endif // QM_MATH_HPP
