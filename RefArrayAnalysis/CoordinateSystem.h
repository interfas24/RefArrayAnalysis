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

public:
	ThreeDCoordSysBase()
		: _quantity({0.0, 0.0, 0.0}) {}
	ThreeDCoordSysBase(double q1, double q2, double q3)
		: _quantity({q1, q2, q3}) {}

	ThreeDCoordSysBase(std::initializer_list<double> li);
	virtual ~ThreeDCoordSysBase() {}

	static ThreeDCoordSysBase OriginalPoint;

	std::vector<double> GetQuantity() const { return { _quantity[0], _quantity[1], _quantity[2] }; }

	void SetQuantity(double q1, double q2, double q3);
	void SetQuantity(const std::vector<double>&);
	void SetQuantity(std::initializer_list<double> li);

	bool operator==(const ThreeDCoordSysBase& other) { return _quantity == other._quantity; }
	bool operator!=(const ThreeDCoordSysBase& other) { return !(*this == other); }
};

class CartesianCS : public ThreeDCoordSysBase
{
public:
	CartesianCS(double x = 0.0, double y = 0.0, double z = 0.0)
		: ThreeDCoordSysBase(x, y, z) {}

	CartesianCS(const CartesianCS&) = default;
	CartesianCS(CartesianCS&&) = default;
	CartesianCS& operator=(const CartesianCS&) = default;
	CartesianCS(const SphericalCS&);

	SphericalCS ToSpherical() const;

	double X() const { return _quantity[0]; }
	CartesianCS& X(double _x) { _quantity[0] = _x; return *this; }

	double Y() const { return _quantity[1]; }
	CartesianCS& Y(double _y) { _quantity[1] = _y; return *this; }

	double Z() const { return _quantity[2]; }
	CartesianCS& Z(double _z) { _quantity[2] = _z; return *this; } 
};

// theta phi : unit rad
class SphericalCS : public ThreeDCoordSysBase
{
public:
	SphericalCS(double r = 0.0, double t = 0.0, double p = 0.0)
		: ThreeDCoordSysBase(r, t, p) {}

	SphericalCS(const SphericalCS&) = default;
	SphericalCS(SphericalCS&&) = default;
	SphericalCS& operator=(const SphericalCS&) = default;
	SphericalCS(const CartesianCS&);

	CartesianCS ToCartesian() const;

	double R() const { return _quantity[0]; }
	SphericalCS& R(double r_) { _quantity[0] = r_; return *this; }

	double Theta() const { return _quantity[1]; }
	SphericalCS& Theta(double t_) { _quantity[1] = t_; return *this; }

	double Phi() const { return _quantity[2]; }
	SphericalCS& Phi(double p_) { _quantity[2] = p_; return *this; }
};

//extern ThreeDCoordSysBase OriginalPoint;