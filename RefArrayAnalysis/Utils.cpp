#include "Utils.h"

using namespace gxx_math;
using namespace std;

double PhysicsConst::LightSpeed = 299792458.0000;
double PhysicsConst::VacuumPermit = 8.8541878176 * 1e-12;
double PhysicsConst::VacuumPermea = 4 * M_PI * 1e-7;
double PhysicsConst::VacuumImped = sqrt(PhysicsConst::VacuumPermea / PhysicsConst::VacuumPermit);

DoubleComplex SpecialFunc::Fresnel(double x)
{
	const double A[] = {
		1.595769140,
		-0.000001702,
		-6.808508854,
		-0.000576361,
		6.920691902,
		-0.016898657,
		-3.050485660,
		-0.075752419,
		0.850663781,
		-0.025639041,
		-0.150230960,
		0.034404779
	};

	const double B[] = {
		-0.000000033,
		4.255387524,
		-0.000092810,
		-7.780020400,
		-0.009520895,
		5.075161298,
		-0.138341947,
		-1.363729124,
		-0.403349276,
		0.702222016,
		-0.216195929,
		0.019547031
	};

	const double C[] = {
		0,
		-0.024933975,
		0.000003936,
		0.005770956,
		0.000689892,
		-0.009497136,
		0.011948809,
		-0.006748873,
		0.000246420,
		0.002102967,
		-0.001217930,
		0.000233939
	};

	const double D[] = {
		0.199471140,
		0.000000023,
		-0.009351341,
		0.000023006,
		0.004851466,
		0.001903218,
		-0.017122914,
		0.029064067,
		-0.027928955,
		0.016497308,
		-0.005598515,
		0.000838386
	};

	const size_t len = ArrayLength(A);

	if (x == 0.0)
		return DoubleComplex();
	else if(x < 0.0)
	{
		x = abs(x);
		x = M_PI_2 * pow(x, 2.);
		DoubleComplex F;
		if (x < 4.) {
			for (size_t k = 0; k < len; k++) {
				F += (DoubleComplex(A[k], B[k]) * pow(x / 4., k));
			}
			return -(F * sqrt(x / 4.) * Exp(- _I_(x)));
		}
		else {
			for (size_t k = 0; k < len; k++) {
				F += (DoubleComplex(C[k], D[k]) * pow(4. / x, k));
			}
			return -(F * sqrt(4. / x) * Exp(- _I_(x)) + DoubleComplex(1, -1) / 2.);
		}
	}
	else
	{
		x = M_PI_2 * pow(x, 2.);
		DoubleComplex F;
		if (x < 4.) {
			for (size_t k = 0; k < len; k++) {
				F += (DoubleComplex(A[k], B[k]) * pow(x / 4., k));
			}
			return (F * sqrt(x / 4.) * Exp(- _I_(x)));
		}
		else {
			for (size_t k = 0; k < len; k++) {
				F += (DoubleComplex(C[k], D[k]) * pow(4. / x, k));
			}
			return (F * sqrt(4. / x) * Exp(- _I_(x)) + DoubleComplex(1, -1) / 2.);
		}
	}

}

bool DiffLTPrecision(const gxx_math::DoubleComplex &lhs, const gxx_math::DoubleComplex &rhs, double pre)
{
	return (abs(lhs.real() - rhs.real()) < pre) && (abs(lhs.imag() - rhs.imag()) < pre);
}

vector<double> Linspace(double start, double stop, size_t points, bool con_end)
{
	vector<double> ret;
	ret.reserve(points);

	double step;
	if (con_end) {
		step = (stop - start) / (points - 1);
	} else {
		step = (stop - start) / points;
	}

	for (size_t i = 0; i < points; i++) {
		ret.push_back(start + i * step);
	}

	return ret;
}

double dB(double v, double coeff)
{
	return v <= 0.0 ? 0.0 : coeff * log10(v);
}

double dB(const gxx_math::DoubleComplex & v, double coeff)
{
	double mag = Abs(v);
	return dB(mag, coeff);
}

double Sinc(double x)
{
	if (x == 0.0)
		return 1.0;
	else
		return sin(x) / x;
}
