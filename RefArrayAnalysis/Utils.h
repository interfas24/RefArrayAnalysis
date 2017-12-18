#pragma once

#include "GxxComplex.h"
#include "gsl/gsl_math.h"
#include <vector>

class PhysicsConst {
public:
	static double LightSpeed;	//unit : m/s
	static double VacuumPermit;	//unit : F/m
	static double VacuumPermea;	//unit : H/m
};

class SpecialFunc
{
public:
	static gxx_math::DoubleComplex Fresnel(double);
};

template <typename T, std::size_t N>
size_t ArrayLength(const T (&arr)[N])
{
	return N;
}

bool DiffLTPrecision(const gxx_math::DoubleComplex &lhs, const gxx_math::DoubleComplex &rhs, double pre);

inline
double RadToDeg(double rad) {
	return rad * 180. / M_PI;
}

inline
double DegToRad(double deg) {
	return deg * M_PI / 180.;
}

inline double Hypot(double x, double y) {
	return sqrt(x*x + y*y);
}

inline double Hypot(double x, double y, double z) {
	return sqrt(x*x + y*y + z*z);
}

std::vector<double> Linspace(double start, double stop, size_t points, bool con_end = true);

/*
template<typename T>
double dB(const T&, v, double coeff = 20.)
{
	return v <= 0.0 ? 0.0 : coeff * log(v);
}

template<typename T>
double dB(const gxx_math::DoubleComplex& v, double coeff)
{
	double mag = Abs(v);
	return dB(mag, coeff);
}
*/

double dB(double v, double coeff = 20.);
double dB(const gxx_math::DoubleComplex &v, double coeff = 20.);

/*
template<typename T>
bool DiffLTPrecision(T lhs, T rhs, double pre) {
	return (abs(lhs - rhs) < pre);
}

template<typename T>
bool 
DiffLTPrecision(const gxx_math::DoubleComplex &lhs, const gxx_math::DoubleComplex &rhs, double pre)
{
	return (abs(lhs.real() - rhs.real()) < pre) && (abs(lhs.imag() - rhs.imag()) < pre);
}
*/