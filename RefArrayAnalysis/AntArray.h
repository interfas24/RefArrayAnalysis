#pragma once

#include "CoordinateSystem.h"
#include "Source.h"
#include "Utils.h"
#include <algorithm>
#include <vector>

class AntArray : public NoCopyable
{
public:
	virtual size_t TotalCells() const = 0;
	virtual double MaxScale() const = 0;

protected:
	std::vector<Source*> _sources;
	std::vector<gxx_math::DoubleComplex> _Emn;
};

//outer vector stands for column along Y axis from -y to +y
//inner vector stands for column along X axis from -x to +x
typedef std::vector<CartesianCS> ArrayDistro;

//solver iterate on array
//units: m

//set reflectarray step
//step1 : setup horn
//step2 : pre-calculate Emn (given phase distribution and csv data file)
//step3 : 
class Reflectarray : public AntArray
{
public:
	Reflectarray(size_t xscale, size_t yscale, double cell_sz)
		: _xscale(xscale), _yscale(yscale), _cell_sz(cell_sz)
	{
	}

	size_t TotalCells() const override { return _xscale * _yscale; }
	double MaxScale() const override { return std::max(_xscale, _yscale) * _cell_sz; }

	//set feed horn and position
	//Feed coordinate and Array coordinate conversion defined by Eulerian Angles (Unit : rad)
	//Ref. IEEE Trans. ON AP. Useful Coordinate Transformation for Antenna Application
	void SetupHorn(PyramidalHorn *horn,
					double alpha, double beta, double gamma, double fdr);
	void UpdateFDR(double fdr) { _fdr = fdr; }
	void UpdateEulerianAngle(double alpha, double beta, double gamma);

	void AddSource(Source*);
	void ResetSource();

	//set TeTm data and phase distribution method(Freq band)
	std::vector<CartesianCS>::iterator Begin() { return _array_info.begin(); }
	std::vector<CartesianCS>::iterator End() { return _array_info.end(); }
	std::vector<CartesianCS>::const_iterator Cbegin() const { return _array_info.cbegin(); }
	std::vector<CartesianCS>::const_iterator Cend() const { return _array_info.cend(); }

protected:
	size_t _xscale;
	size_t _yscale;
	double _cell_sz;
	ArrayDistro _array_info;

	PyramidalHorn *_py_horn;
	double _fdr;
	CartesianCS _feed_pos;

private:
	//step1
	void _compute_feed_pos(double alpha, double beta, double gamma);

	//step2
	void _recompute_Emn();
};

class RectRefArray : public Reflectarray 
{
public:

	RectRefArray(size_t xscale, size_t yscale, double cell_sz)
		: Reflectarray(xscale, yscale, cell_sz) 
	{
		_array_info = initArray();
	}

protected:
	ArrayDistro initArray();
};

class SquareRefArray : public RectRefArray
{
public:
	SquareRefArray(size_t scale, double cell_sz)
		: RectRefArray(scale, scale, cell_sz) {}
};
