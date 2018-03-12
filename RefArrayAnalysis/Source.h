#pragma once

#include <vector>
#include <utility>
#include "CoordinateSystem.h"
#include "GxxComplex.h"
#include "Utils.h"

//unit : Hz;
class Source : public NoCopyable {
public:
	Source(double freq) { _reset_freq(freq); }

	virtual std::vector<gxx_math::DoubleComplex>
		RadiationPatternAt(const SphericalCS&) = 0;
	virtual std::vector<gxx_math::DoubleComplex>
		RadiationPatternAt(const CartesianCS&) = 0;
	virtual double GetTotalRadiationPower() const { return _total_radiation_power; }
	virtual void ResetSource(double freq) { _reset_freq(freq); }

	double GetFreq() const { return _freq; }
	double GetLambda() const { return _lambda; }
	double GetWaveNumber0() const { return _k0; }
	void Place(const SourcePosition&);
	SourcePosition GetPosition() const { return _source_position; }

protected:
	double _freq;
	double _lambda;
	double _k0;
	double _total_radiation_power;
	SourcePosition _source_position;

	void _reset_freq(double freq)
	{
		_freq = freq;
		//_lambda = PhysicsConst::LightSpeed / _freq;
		_lambda = 3e8 / _freq;
		_k0 = 2 * M_PI / _lambda;
	}
};

//polarization along Y axis in hron CS
class PyramidalHorn : public Source
{
public:
	//r1, r2, a, b, a1, b1 real size
	PyramidalHorn(double r1, double r2, double a, double b, double a1, double b1, double freq, double E0);

	//Spherical CS get rE_r_theta_phi
	std::vector<gxx_math::DoubleComplex> 
	RadiationPatternAt(const SphericalCS&) override;
	//Cartesian CS get rE_x_y_z
	std::vector<gxx_math::DoubleComplex>
	RadiationPatternAt(const CartesianCS&) override;

	//double GetTotalRadiationPower() const { return _total_radiation_power; }

	//R set to 100meter
	double HornIntegralFunc(double t, double p);
	std::pair<double, double> GetThetaRange() const { return _theta_range; }
	std::pair<double, double> GetPhiRange() const { return _phi_range; }

	void ResetSource(double freq) override;

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

	const std::pair<double, double> _theta_range;
	const std::pair<double, double> _phi_range;

	//internal radiation pattern method
	//return Etheta and Ephi pair, far field NO Er
	void _reset_horn(double freq);
	std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> _radiation_pattern_at(const SphericalCS&);

	double _compute_total_radiation_power();
};
