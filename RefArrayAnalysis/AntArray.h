#pragma once

#include "CoordinateSystem.h"
#include <vector>

class AntArrayBase {
};

//outer vector stands for column along Y axis from -y to +y
//inner vector stands for column along X axis from -x to +x
typedef std::vector<CartesianCS> ArrayDistro;

//solver iterate on array
//units: m
class Reflectarray : public AntArrayBase 
{
public:
	Reflectarray(size_t xscale, size_t yscale, double cell_sz, const CartesianCS & fp);
	Reflectarray(const Reflectarray &) = delete;

	virtual size_t totalCells() const {
		return _xscale * _yscale;
	}

	std::vector<CartesianCS>::iterator begin() {
		return _array_info.begin();
	}
	std::vector<CartesianCS>::iterator end() {
		return _array_info.end();
	}

	std::vector<CartesianCS>::const_iterator begin() const {
		return _array_info.cbegin();
	}
	std::vector<CartesianCS>::const_iterator end() const {
		return _array_info.cend();
	}

protected:
	size_t _xscale;
	size_t _yscale;
	double _cell_sz;
	CartesianCS _feed_pos;
	ArrayDistro _array_info;

	virtual ArrayDistro initArray();
};

class RectRefArray : public Reflectarray 
{
public:

	RectRefArray(size_t xscale, size_t yscale, double cell_sz, const CartesianCS &fp)
		: Reflectarray(xscale, yscale, cell_sz, fp) {}

protected:
	ArrayDistro initArray();
};

class SquareRefArray : public RectRefArray
{
public:
	SquareRefArray(size_t scale, double cell_sz, const CartesianCS &fp)
		: RectRefArray(scale, scale, cell_sz, fp) {}
};
