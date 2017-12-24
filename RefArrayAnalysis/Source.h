#pragma once

#include "CoordinateSystem.h"
#include <utility>
#include "GxxComplex.h"
#include <vector>
#include "Utils.h"

class SourceBase {
	
public:
	//unit : Hz;
	SourceBase(double freq) {
		resetFreq(freq);
	}

	double getFreq() const {
		return _freq;
	}

	double getLambda() const {
		return _lambda;
	}

	double getWaveNumber0() const {
		return _k0;
	}

	void resetFreq(double freq) {
		_freq = freq;
		//_lambda = PhysicsConst::LightSpeed / _freq;
		_lambda = 3e8 / _freq;
		_k0 = 2 * M_PI / _lambda;
	}

protected:
	double _freq;
	double _lambda;
	double _k0;
};

//polarization along Y axis in hron axes
class PyramidalHorn : public SourceBase
{
public:
	//r1, r2, a, b, a1, b1 real size
	PyramidalHorn(double r1, double r2, double a, double b, double a1, double b1, double freq, double E0);

	//Spherical CS get rE_r_theta_phi
	std::vector<gxx_math::DoubleComplex> RadiPatternAt(const SphericalCS&);
	//Cartesian CS get rE_x_y_z
	std::vector<gxx_math::DoubleComplex> RadiPatternAt(const CartesianCS&);

	//R set to 100meter
	double hornIntegralFunc(double t, double p);
	double getTotalPowerRadi() const {
		return _totalPowerRadi;
	}

	void setReserveNumber(double tmp) {
		_reserve_num = tmp;
	}

	double getReserveNumber() const {
		return _reserve_num;
	}

	std::pair<double, double> getThetaRange() const {
		return _theta_range;
	}

	std::pair<double, double> getPhiRange() const {
		return _phi_range;
	}

	void resetHorn(double freq);

private:
	double _r1;
	double _r2;
	double _a;
	double _b;
	double _a1;
	double _b1;

	const double _r1_real;
	const double _r2_real;
	const double _a_real;
	const double _b_real;
	const double _a1_real;
	const double _b1_real;
	
	double _E0;

	double _totalPowerRadi;
	double _reserve_num;

	const std::pair<double, double> _theta_range;
	const std::pair<double, double> _phi_range;

	//internal radiation pattern method
	//return Etheta and Ephi pair, far field NO Er
	std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> _radiation_pattern_at(const SphericalCS&);

	double _compute_total_power_radiation();
};
