#include "AntArray.h"
#include <stdexcept>

using namespace std;

namespace antarray_internal {
	
}

ArrayDistro RectRefArray::initArray()
{
	if (_xscale < 2 || _yscale < 2) {
		throw logic_error("Too small array");
	}

	ArrayDistro ret;
	ret.reserve(_xscale * _yscale);

	vector<double> xlist(_xscale);
	vector<double> ylist(_yscale);
	for (size_t i = 0; i < _xscale; i++) {
		int left = - static_cast<int>(_xscale / 2);
		xlist[i] = left * _cell_sz + _cell_sz / 2. + i * _cell_sz;
	}
	for (size_t i = 0; i < _yscale; i++) {
		int left = - static_cast<int>(_yscale / 2);
		ylist[i] = left * _cell_sz + _cell_sz / 2. + i * _cell_sz;
	}

	for (size_t i = 0; i < _yscale; i++) {
		for (size_t j = 0; j < _xscale; j++) {
			ret.push_back(CartesianCS(xlist[j], ylist[i], 0.0));
		}
	}

	return ret;
}

void Reflectarray::setupHorn(PyramidalHorn * horn, double alpha, double beta, double gamma, double fdr)
{
	_py_horn = horn;
	_fdr = fdr;
	_compute_feed_pos(alpha, beta, gamma);
}

void Reflectarray::updateEulerianAngle(double alpha, double beta, double gamma)
{
	_compute_feed_pos(alpha, beta, gamma);
}

void Reflectarray::_compute_feed_pos(double alpha, double beta, double gamma)
{
}


