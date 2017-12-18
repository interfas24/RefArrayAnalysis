#include "AntArray.h"

Reflectarray::Reflectarray(size_t xscale, size_t yscale, double cell_sz, const CartesianCS & fp = CartesianCS(0., 0., 0.))
	: _xscale(xscale), _yscale(yscale), _cell_sz(cell_sz), _feed_pos(fp)
{
	_array_info = initArray();
}

ArrayDistro RectRefArray::initArray()
{
	return ArrayDistro();
}
