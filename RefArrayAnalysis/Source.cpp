#include "Source.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_integration.h"
#include "GxxComplex.h"
#include "Utils.h"
#include <cassert>

using namespace std;
using namespace gxx_math;

namespace source_internal {
	
	const double IntegrationPre = 1e-3;

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

	double horn_radi_x_func(double x, void *p)
	{
		PyramidalHorn *horn = static_cast<PyramidalHorn*>(p);
		return horn->hornIntegralFunc(x, horn->getReserveNumber());
	}

	double horn_radi_y_func(double y, void *p)
	{
		PyramidalHorn *horn = static_cast<PyramidalHorn*>(p);
		horn->setReserveNumber(y);
		gsl_function F;
		F.function = horn_radi_x_func;
		F.params = p;
		double result, abserr;
		size_t eval;
		gsl_integration_qng(&F, horn->getThetaRange().first,
			horn->getThetaRange().second,
			IntegrationPre, IntegrationPre, &result, &abserr, &eval);
		return result;
	}

}

PyramidalHorn::PyramidalHorn(double r1, double r2, double a, double b, double a1, double b1, double freq, double E0)
	: SourceBase(freq),
	_r1_real(r1),
	_r2_real(r2),
	_a_real(a),
	_b_real(b),
	_a1_real(a1),
	_b1_real(b1),
	_theta_range({0, M_PI_2}),
	_phi_range({0, 2 * M_PI}),
	_E0(E0)
{
	resetHorn(freq);
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

double PyramidalHorn::hornIntegralFunc(double t, double p)
{
	double R = 100.0;
	vector<DoubleComplex> rtp = RadiPatternAt(SphericalCS(R, t, p));
	DoubleComplex ee = Sqrt(rtp[1] * rtp[1] + rtp[2] * rtp[2]);
	return Abs(ee * ee * R * R * sin(t));
}

void PyramidalHorn::resetHorn(double freq)
{
	resetFreq(freq);
	_r1 = _r1_real / getLambda();
	_r2 = _r2_real / getLambda();
	_a = _a_real / getLambda();
	_b = _b_real / getLambda();
	_a1 = _a1_real / getLambda();
	_b1 = _b1_real / getLambda();
	_totalPowerRadi = _compute_total_power_radiation();
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

double PyramidalHorn::_compute_total_power_radiation()
{
	gsl_function F;
	F.function = source_internal::horn_radi_y_func;
	F.params = this;
	double result, abserr;
	size_t eval;
	gsl_integration_qng(&F, this->getPhiRange().first,
		this->getPhiRange().second,
		source_internal::IntegrationPre, 
		source_internal::IntegrationPre, &result, &abserr, &eval);
	return result;
}
