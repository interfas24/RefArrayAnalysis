#include "Source.h"
#include "gsl/gsl_math.h"
#include "GxxComplex.h"
#include "Utils.h"
#include <cassert>

using namespace std;
using namespace gxx_math;

namespace source_internal {
	
	vector<DoubleComplex> mtrx_mul_vector(const vector<vector<double>>& mtx,
		const vector<DoubleComplex>& v)
	{
		assert(mtx.size() == 3 &&
			mtx[0].size() == 3 &&
			v.size() == 3);
		vector<DoubleComplex> ret;
		ret.reserve(3);
		for (size_t i = 0; i < 3; i++) {
			DoubleComplex sum(0., 0.);
			for (size_t j = 0; j < 3; j++) {
				sum += mtx[i][j] * v[j];
			}
			ret.push_back(sum);
		}
		return ret;
	}

}

PyramidalHorn::PyramidalHorn(double r1, double r2, double a, double b, double a1, double b1,
		double E0, double k0)
		: _r1(r1), _r2(r2), _a(a), _b(b), _a1(a1), _b1(b1), _E0(E0), _k0(k0) 
{
}

vector<DoubleComplex> PyramidalHorn::RadiPatternAt(const SphericalCS &s)
{
	pair<DoubleComplex, DoubleComplex> e_eta_phi = _radiation_pattern_at(s);
	return { DoubleComplex(0.0, 0.0), e_eta_phi.first, e_eta_phi.second };
}

vector<DoubleComplex> PyramidalHorn::RadiPatternAt(const CartesianCS &c)
{
	SphericalCS s = c.toSpherical();
	vector<DoubleComplex> rErtp = RadiPatternAt(s);

	//from spherical CS to cartesian CS
	//Ref. IEEE Trans. AP. VOL. AP-27, No.4,JULY 1979
	//Useful Coordinate Transformations for Antenna Applications
	double t = s.theta(), p = s.phi();
	vector<vector<double>> convMtrx = {
		{ sin(t)*cos(p), cos(t)*cos(p), -sin(p) },
		{ sin(t)*sin(p), cos(t)*sin(p), cos(p) },
		{ cos(t), -sin(t), 0.0 }
	};
	return source_internal::mtrx_mul_vector(convMtrx, rErtp);
}

std::pair<DoubleComplex, DoubleComplex>
PyramidalHorn::_radiation_pattern_at(const SphericalCS &s)
{
	double k = 2 * M_PI;
	double ky = k * sin(s.theta()) * sin(s.phi());
	double t1 = sqrt(1. / (M_PI * k * _r1)) * (-k * _b1 / 2. - ky * _r1);
	double t2 = sqrt(1. / (M_PI * k * _r1)) * (k * _b1 / 2. - ky * _r1);

	double kxp = k * sin(s.theta()) * cos(s.phi()) + M_PI / _a1;
	double kxdp = k * sin(s.theta()) * cos(s.phi()) - M_PI / _a1;
	double t1p = sqrt(1. / (M_PI * k * _r2)) * (-k * _a1 / 2. - kxp * _r2);
	double t2p = sqrt(1. / (M_PI * k * _r2)) * (k * _a1 / 2. - kxp * _r2);

	double t1dp = sqrt(1. / (M_PI * k * _r2)) * (-k * _a1 / 2. - kxdp * _r2);
	double t2dp = sqrt(1. / (M_PI * k * _r2)) * (k * _a1 / 2. - kxdp * _r2);

	DoubleComplex I1 = .5 * sqrt(M_PI * _r2 / k) * (
			Exp(_I_(kxp*kxp * _r2 / (2*k))) * (SpecialFunc::Fresnel(t2p) - SpecialFunc::Fresnel(t1p))
			+ Exp(_I_(kxdp*kxdp *_r2 / (2*k))) * (SpecialFunc::Fresnel(t2dp) - SpecialFunc::Fresnel(t1dp))
		);

	DoubleComplex I2 = sqrt(M_PI * _r1 / k) * Exp(_I_(ky*ky * _r1 / (2 * k))) *
					(SpecialFunc::Fresnel(t2) - SpecialFunc::Fresnel(t1));

	DoubleComplex Etheta, Ephi;

	k = _k0;

	Etheta = _I_(k * _E0) * Exp(_I_(-k * s.r())) / (4 * M_PI * s.r()) *
		(sin(s.phi()) * (1. + cos(s.theta())) * I1 * I2);

	Ephi = _I_(k * _E0) * Exp(_I_(-k * s.r())) / (4 * M_PI * s.r()) *
		(cos(s.phi()) * (1. + cos(s.theta())) * I1 * I2);

	return {Etheta, Ephi};
}
