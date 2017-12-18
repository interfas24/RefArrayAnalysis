#pragma once

#include "CoordinateSystem.h"
#include <vector>

class AntArrayBase {
};

//outer vector stands for column along Y axis from -y to +y
//inner vector stands for column along X axis from -x to +x
typedef std::vector<std::vector<CartesianCS>> ArrayDistro;

//solver iterate on array
//units: mm
class Reflectarray : public AntArrayBase 
{
public:
	Reflectarray(size_t xscale, size_t yscale, double cell_sz, const CartesianCS & fp);
	Reflectarray(const Reflectarray &) = delete;

protected:
	size_t _xscale;
	size_t _yscale;
	double _cell_sz;
	CartesianCS _feed_pos;
	ArrayDistro _array_info;

	virtual ArrayDistro initArray() = 0;
};

class RectRefArray : public Reflectarray 
{
public:

	RectRefArray(size_t xscale, size_t yscale, double cell_sz, const CartesianCS &fp)
		: Reflectarray(xscale, yscale, cell_sz, fp) {}

	ArrayDistro initArray();

private:
};
