#pragma once

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <iostream>
#include <complex>

/*
GNU Scientific Library C++ wrapper
*/

namespace gxx_math {

#define _RE 0
#define _IM 1

template <typename T, typename _complex_T>
class GxxComplexBase : public _complex_T
{
public:
	typedef T value_type;

	explicit
	GxxComplexBase(const T& re_ = T(), const T& im_ = T())
		: _complex_T({re_, im_}) 
	{
	}
	GxxComplexBase(const std::complex<T>& z)
		: _complex_T({ z.real(), z.imag() })
	{
	}
	GxxComplexBase(const _complex_T& z)
		: _complex_T(z)
	{
	}

	T real(const T& re) {
		return (this->dat[_RE] = re);
	}

	T real() const {
		return this->dat[_RE];
	}

	T imag(const T& im) {
		return (this->dat[_IM] = im);
	}

	T imag() const {
		return this->dat[_IM];
	}
};

/*
//float specialize template
template<>
class GxxComplex <float> : public GxxComplexBase<float, gsl_complex_float>
{
};

//double specialize template
template<>
class GxxComplex <double> : public GxxComplexBase<double, gsl_complex>
{
};

//long double specialize template
template<>
class GxxComplex <long double> : public GxxComplexBase<long double, gsl_complex_long_double>
{
};
*/

//general complex
template <typename T>
class GxxComplex : public GxxComplexBase<T, gsl_complex>
{
public:
	typedef GxxComplexBase<T, gsl_complex> _Base;

	//ctor
	explicit
	GxxComplex<T>(const T& re_ = T(), const T& im_ = T())
		: _Base(re_, im_)
	{
	}

	GxxComplex<T>& operator+=(const T& other) 
	{
		gsl_complex z = { this->real(), this->imag() };
		z = gsl_complex_add_real(z, other);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}
	GxxComplex<T>& operator+=(const GxxComplex<T>& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		gsl_complex rhs = { other.real(), other.imag() };
		z = gsl_complex_add(z, rhs);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}

	GxxComplex<T>& operator-=(const T& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		z = gsl_complex_sub_real(z, other);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}
	GxxComplex<T>& operator-=(const GxxComplex<T>& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		gsl_complex rhs = { other.real(), other.imag() };
		z = gsl_complex_sub(z, rhs);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}

	GxxComplex<T>& operator*=(const T& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		z = gsl_complex_mul_real(z, other);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}
	GxxComplex<T>& operator*=(const GxxComplex<T>& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		gsl_complex rhs = { other.real(), other.imag() };
		z = gsl_complex_mul(z, rhs);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}

	GxxComplex<T>& operator/=(const T& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		z = gsl_complex_div_real(z, other);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}
	GxxComplex<T>& operator/=(const GxxComplex<T>& other)
	{
		gsl_complex z = { this->real(), this->imag() };
		gsl_complex rhs = { other.real(), other.imag() };
		z = gsl_complex_div(z, rhs);
		this->real(GSL_REAL(z));
		this->imag(GSL_IMAG(z));
		return *this;
	}
};

//construct func
template<typename T>
GxxComplex<T>
Polar(const T& mag_, const T& rad_)
{
	gsl_complex z = gsl_complex_polar(mag_, rad_);
	return GxxComplex<T>(GSL_REAL(z), GSL_IMAG(z));
}

/* Properties of complex numbers */
template<typename T>
T Arg(const GxxComplex<T>& z)
{
	return gsl_complex_arg(z);
}
template<typename T>
T Abs(const GxxComplex<T>& z)
{
	return gsl_complex_abs(z);
}
template<typename T>
T Norm(const GxxComplex<T>& z)
{
	return gsl_complex_abs2(z);
}
template<typename T>
GxxComplex<T>
Conj(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_conjugate(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

/* Elementary Complex Functions */
template<typename T>
GxxComplex<T>
Sqrt(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_sqrt(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Sqrt(const T& z)
{
	gsl_complex c = gsl_complex_sqrt_real(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Pow(const GxxComplex<T>& a, const GxxComplex<T>& b)
{
	gsl_complex c = gsl_complex_pow(a, b);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Pow(const GxxComplex<T>& a, const T& b)
{
	gsl_complex c = gsl_complex_pow_real(a, b);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

// add parentheses when using ^ for precedence
// (a ^ b) == Pow(a, b)
template<typename T>
GxxComplex<T>
operator^(const GxxComplex<T>& a, const GxxComplex<T>& b)
{
	return Pow(a, b);
}
template<typename T>
GxxComplex<T>
operator^(const GxxComplex<T>& a, const T& b)
{
	return Pow(a, b);
}

template<typename T>
GxxComplex<T>
Exp(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_exp(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Log(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_log(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Log10(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_log10(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
LogB(const GxxComplex<T>& a, const GxxComplex<T>& b)
{
	gsl_complex c = gsl_complex_log_b(a, b);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
LogB(const GxxComplex<T>& a, const T& b)
{
	gsl_complex c = gsl_complex_log_b(a, gsl_complex_rect(b, 0.0));
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

/* Complex Trigonometric Functions */
template<typename T>
GxxComplex<T>
Sin(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_sin(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Cos(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_cos(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Sec(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_sec(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Csc(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_csc(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Tan(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_tan(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Cot(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_cot(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

/* Inverse Complex Trigonometric Functions */
template<typename T>
GxxComplex<T>
ArcSin(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arcsin(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcSin(const T& a)
{
	gsl_complex c = gsl_complex_arcsin_real(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCos(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arccos(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCos(const T& a)
{
	gsl_complex c = gsl_complex_arccos_real(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcSec(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arcsec(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcSec(const T& a)
{
	gsl_complex c = gsl_complex_arcsec_real(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCsc(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arccsc(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCsc(const T& a)
{
	gsl_complex c = gsl_complex_arccsc_real(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcTan(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arctan(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCot(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arccot(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

/* Complex Hyperbolic Functions */
template<typename T>
GxxComplex<T>
Sinh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_sinh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Cosh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_cosh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Sech(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_sech(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Csch(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_csch(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Tanh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_tanh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
Coth(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_coth(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

/* Inverse Complex Hyperbolic Functions */
#if 0
gsl_complex gsl_complex_arcsinh (gsl_complex a);  /* r=arcsinh(a) */
gsl_complex gsl_complex_arccosh (gsl_complex a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arccosh_real (double a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arcsech (gsl_complex a);  /* r=arcsech(a) */
gsl_complex gsl_complex_arccsch (gsl_complex a);  /* r=arccsch(a) */
gsl_complex gsl_complex_arctanh (gsl_complex a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arctanh_real (double a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arccoth (gsl_complex a);  /* r=arccoth(a) */
#endif

template<typename T>
GxxComplex<T>
ArcSinh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arcsinh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCosh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arccosh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCosh(const T& a)
{
	gsl_complex c = gsl_complex_arccosh_real (a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcSech(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arcsech(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCsch(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arccsch(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcTanh(const GxxComplex<T>& a)
{
	gsl_complex c = gsl_complex_arctanh(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcTanh(const T& a)
{
	gsl_complex c = gsl_complex_arctanh_real(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}
template<typename T>
GxxComplex<T>
ArcCoth(const T& a)
{
	gsl_complex c = gsl_complex_arccoth(a);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

//unary opts
template<typename T>
GxxComplex<T>
operator+(const GxxComplex<T>& z)
{
	return z;
}
template<typename T>
GxxComplex<T>
operator-(const GxxComplex<T>& z)
{
	gsl_complex c = gsl_complex_negative(z);
	return GxxComplex<T>(GSL_REAL(c), GSL_IMAG(c));
}

//binary opts
template<typename T>
GxxComplex<T> 
operator+(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs);
	z += rhs;
	return z;
}
template<typename T>
GxxComplex<T> 
operator+(const GxxComplex<T>& lhs, const T& rhs)
{
	GxxComplex<T> z(lhs);
	z += rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator+(const T& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs, 0.);
	z += rhs;
	return z;
}

template<typename T>
GxxComplex<T>
operator-(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs);
	z -= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator-(const GxxComplex<T>& lhs, const T& rhs)
{
	GxxComplex<T> z(lhs);
	z -= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator-(const T& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs, 0.);
	z -= rhs;
	return z;
}

template<typename T>
GxxComplex<T>
operator*(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs);
	z *= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator*(const GxxComplex<T>& lhs, const T& rhs)
{
	GxxComplex<T> z(lhs);
	z *= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator*(const T& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs, 0.);
	z *= rhs;
	return z;
}

template<typename T>
GxxComplex<T>
operator/(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs);
	z /= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator/(const GxxComplex<T>& lhs, const T& rhs)
{
	GxxComplex<T> z(lhs);
	z /= rhs;
	return z;
}
template<typename T>
GxxComplex<T>
operator/(const T& lhs, const GxxComplex<T>& rhs)
{
	GxxComplex<T> z(lhs, 0.);
	z /= rhs;
	return z;
}

template<typename T>
std::ostream& operator<<(std::ostream& o, const GxxComplex<T>& z)
{
	o << "(" << z.real() << "," << z.imag() << ")";
	return o;
}

//relation opts
template<typename T>
bool operator==(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	return GSL_COMPLEX_EQ(lhs, rhs);
}

template<typename T>
bool operator!=(const GxxComplex<T>& lhs, const GxxComplex<T>& rhs)
{
	return !GSL_COMPLEX_EQ(lhs, rhs);
}

template<typename T>
bool operator==(const std::complex<T>& lhs, const GxxComplex<T>& rhs)
{
	return (lhs.real() == rhs.real()) && (lhs.imag() == rhs.imag());
}
template<typename T>
bool operator!=(const std::complex<T>& lhs, const GxxComplex<T>& rhs)
{
	return !((lhs.real() == rhs.real()) && (lhs.imag() == rhs.imag()));
}


typedef GxxComplex<float> FloatComplex;
typedef GxxComplex<double> DoubleComplex;
typedef GxxComplex<long double> LDoubleComplex;

//helpers
template<typename T>
GxxComplex<T>
_I_(const T& im_)
{
	return GxxComplex<T>(0., im_);
}


}


