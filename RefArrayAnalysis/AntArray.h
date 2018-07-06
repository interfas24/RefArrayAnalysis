#pragma once

#include "CoordinateSystem.h"
#include "Source.h"
#include "Utils.h"
#include <algorithm>
#include <vector>
#include <memory>

class AntArray : public NoCopyable
{
public:
	virtual size_t TotalCells() const = 0;
	virtual double MaxScale() const = 0;

protected:
	std::vector<std::shared_ptr<Source>> _sources;
	std::vector<std::vector<gxx_math::DoubleComplex>> _incident_field;
};

//outer vector stands for column along Y axis from -y to +y
//inner vector stands for column along X axis from -x to +x
typedef std::vector<CartesianCS> ArrayDistro;

//solver iterate on array
//units: m

//set reflectarray step
//step1 : setup horn
//step2 : pre-calculate Emn (given phase distribution and csv data file)
class Reflectarray : public AntArray
{
public:
	Reflectarray(size_t xscale, size_t yscale, double cell_sz)
		: _xscale(xscale), _yscale(yscale), _cell_sz(cell_sz)
	{
		_incident_field.resize(TotalCells());
		for (size_t i = 0; i < TotalCells(); i++)
			_incident_field[i].resize(3);
		_array_info.resize(TotalCells());
	}

	std::vector<double> Tests();
	std::vector<std::vector<gxx_math::DoubleComplex>>
		GetIncidentField();
	ArrayDistro GetArrayDistro() { return _array_info; }

	size_t TotalCells() const override { return _xscale * _yscale; }
	double MaxScale() const override { return std::max(_xscale, _yscale) * _cell_sz; }
	size_t XScale() const { return _xscale; }
	size_t YScale() const { return _yscale; }
	double CellSize() const { return _cell_sz; }

	void AddSource(std::shared_ptr<Source> src);
	std::shared_ptr<Source> GetSource() const { return _sources[0]; }
	void ResetSource();
	std::vector<CartesianCS> GetSourcesPos() const;

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

private:
	void _recompute_incident_field();
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
