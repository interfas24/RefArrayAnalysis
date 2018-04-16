#pragma once

#include "GxxComplex.h"
#include "gsl/gsl_math.h"
#include <vector>

class NoCopyable {
public:
	NoCopyable() = default;
	virtual ~NoCopyable() = default;
	NoCopyable(const NoCopyable&) = delete;
	NoCopyable& operator=(const NoCopyable&) = delete;
};

class PhysicsConst {
public:
	static const double LightSpeed;	//unit : m/s
	static const double VacuumPermit;	//unit : F/m
	static const double VacuumPermea;	//unit : H/m
	static const double VacuumImped;	//unit : ohm
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

// not very suitable
struct SourcePosition
{
	double Alpha;
	double Beta;
	double Gamma;
	double Distance;
};

// E or H field, in either cartesian or spherical system
// Cartesian : X Y Z
// Spherical : R T P
struct EOrHField
{
	gxx_math::DoubleComplex A;
	gxx_math::DoubleComplex B;
	gxx_math::DoubleComplex C;
};

struct CellResponse
{
	gxx_math::DoubleComplex TETE;
	gxx_math::DoubleComplex TETM;
	gxx_math::DoubleComplex TMTE;
	gxx_math::DoubleComplex TMTM;
};

typedef size_t(*PhaseFuzzifier)(const std::vector<double>&, double);
size_t PhaseInterpFuzzifier(const std::vector<double>&, double);
size_t PhaseStepFuzzifier(const std::vector<double>&, double);

bool DiffLTPrecision(const gxx_math::DoubleComplex &lhs, const gxx_math::DoubleComplex &rhs, double pre);

inline double RadToDeg(double rad) { return rad * 180. / M_PI; }

inline double DegToRad(double deg) { return deg * M_PI / 180.; }

inline double mm2m(double mm) { return mm / 1000.0; }

inline double Hypot(double x, double y) { return sqrt(x*x + y*y); }

inline double Hypot(double x, double y, double z) { return sqrt(x*x + y*y + z*z); }

gxx_math::DoubleComplex Hypot(gxx_math::DoubleComplex x,
								gxx_math::DoubleComplex y,
								gxx_math::DoubleComplex z);

std::vector<double> Linspace(double start, double stop, size_t points, bool con_end = true);

double dB(double v, double coeff = 20.);
double dB(const gxx_math::DoubleComplex &v, double coeff = 20.);

double Sinc(double x);

