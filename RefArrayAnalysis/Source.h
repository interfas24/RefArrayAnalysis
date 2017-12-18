#pragma once

#include "CoordinateSystem.h"
#include <utility>
#include "GxxComplex.h"
#include <vector>

class SourceBase {
};

//r1, r2, a, b, a1, b1 normallize to wavelength
//polarization along Y axis in hron axes
class PyramidalHorn : public SourceBase
{
public:
	PyramidalHorn(double r1, double r2, double a, double b, double a1, double b1,
		double E0, double k0);
	
	//Spherical CS get rE_r_theta_phi
	std::vector<gxx_math::DoubleComplex> RadiPatternAt(const SphericalCS&);
	//Cartesian CS get rE_x_y_z
	std::vector<gxx_math::DoubleComplex> RadiPatternAt(const CartesianCS&);

private:
	double _r1;
	double _r2;
	double _a;
	double _b;
	double _a1;
	double _b1;
	double _E0;
	double _k0;

	//internal radiation pattern method
	//return Etheta and Ephi pair, far field NO Er
	std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> _radiation_pattern_at(const SphericalCS&);
};
