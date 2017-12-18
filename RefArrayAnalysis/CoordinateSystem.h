#pragma once

#include <array>
#include <vector>
#include <initializer_list>

//three dimensional coordinate system


class CartesianCS;
class SphericalCS;

class ThreeDCoordSysBase
{
protected:
	std::array<double, 3> _quantity;

	ThreeDCoordSysBase()
		: _quantity({0.0, 0.0, 0.0}) {}
	ThreeDCoordSysBase(double q1, double q2, double q3)
		: _quantity({q1, q2, q3}) {}

	ThreeDCoordSysBase(std::initializer_list<double> li);

public:
	std::vector<double> getQuantity() const {
		return { _quantity[0], _quantity[1], _quantity[2] };
	}

	void setQuantity(double q1, double q2, double q3);
	void setQuantity(const std::vector<double>&);
	void setQuantity(std::initializer_list<double> li);
};

class CartesianCS : public ThreeDCoordSysBase
{
public:
	CartesianCS(double x, double y, double z)
		: ThreeDCoordSysBase(x, y, z) {}

	CartesianCS(const CartesianCS&) = default;
	CartesianCS(const SphericalCS&);

	SphericalCS toSpherical() const;

	double x() const {
		return _quantity[0];
	}

	CartesianCS& x(double _x) {
		_quantity[0] = _x;
		return *this;
	}

	double y() const {
		return _quantity[1];
	}

	CartesianCS& y(double _y) {
		_quantity[1] = _y;
		return *this;
	}

	double z() const {
		return _quantity[2];
	}

	CartesianCS& z(double _z) {
		_quantity[2] = _z;
		return *this;
	}
};

// theta phi : unit rad
class SphericalCS : public ThreeDCoordSysBase
{
public:
	SphericalCS(double r, double t, double p)
		: ThreeDCoordSysBase(r, t, p) {}

	SphericalCS(const SphericalCS&) = default;
	SphericalCS(const CartesianCS&);

	CartesianCS toCartesian() const;

	double r() const {
		return _quantity[0];
	}

	SphericalCS& r(double r_) {
		_quantity[0] = r_;
		return *this;
	}

	double theta() const {
		return _quantity[1];
	}

	SphericalCS& theta(double t_) {
		_quantity[1] = t_;
		return *this;
	}

	double phi() const {
		return _quantity[2];
	}

	SphericalCS& phi(double p_) {
		_quantity[2] = p_;
		return *this;
	}
};


