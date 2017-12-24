#include "CoordinateSystem.h"
#include "Utils.h"
#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace std;

ThreeDCoordSysBase OriginalPoint;

namespace cs_internal {
		
	vector<double> car2sph(const vector<double>& sph) {
		assert(sph.size() == 3);
		vector<double> ret(sph.size());
		ret[0] = sph[0] * sin(sph[1]) * cos(sph[2]);
		ret[1] = sph[0] * sin(sph[1]) * sin(sph[2]);
		ret[2] = sph[0] * cos(sph[1]);
		return ret;
	}

	vector<double> sph2car(const vector<double>& car) {
		assert(car.size() == 3);
		vector<double> ret(car.size());
		ret[0] = Hypot(car[0], car[1], car[2]);
		ret[1] = acos(car[2] / ret[0]);
		ret[2] = atan2(car[1], car[0]);
		return ret;
	}

}

ThreeDCoordSysBase::ThreeDCoordSysBase(initializer_list<double> li)
{
	setQuantity(li);
}

void ThreeDCoordSysBase::setQuantity(double q1, double q2, double q3)
{
	setQuantity({ q1, q2, q3 });
}

void ThreeDCoordSysBase::setQuantity(const std::vector<double>& dat)
{
	assert(dat.size() == 3);
	setQuantity({ dat[0], dat[1], dat[2] });
}

void ThreeDCoordSysBase::setQuantity(initializer_list<double> li)
{
	assert(li.size() == 3);
	copy(li.begin(), li.end(), _quantity.begin());
}

CartesianCS::CartesianCS(const SphericalCS &sph)
{
	setQuantity(cs_internal::sph2car(sph.getQuantity()));
}

SphericalCS CartesianCS::toSpherical() const
{
	return SphericalCS(*this);
}

SphericalCS::SphericalCS(const CartesianCS &car)
{
	setQuantity(cs_internal::car2sph(car.getQuantity()));
}

CartesianCS SphericalCS::toCartesian() const
{
	return CartesianCS(*this);
}