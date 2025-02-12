#include "AntArray.h"
#include <stdexcept>
#include <cassert>

using namespace std;

namespace antarray_internal 
{
	vector<vector<double>> R2F(double alpha, double beta, double gamma)
	{
		vector<vector<double>> ret = {
			{0, 0, 0},
			{0, 0, 0},
			{0, 0, 0}
		};
		vector<vector<double>> tmp = ret;
		vector<vector<double>> mtx1 = {
			{cos(gamma), sin(gamma), 0},
			{-sin(gamma), cos(gamma), 0},
			{0, 0, 1}
		};
		vector<vector<double>> mtx2 = {
			{1, 0, 0},
			{0, cos(beta), sin(beta)},
			{0, -sin(beta), cos(beta)}
		};
		vector<vector<double>> mtx3 = {
			{cos(alpha), sin(alpha), 0},
			{-sin(alpha), cos(alpha), 0},
			{0, 0, 1}
		};

		for (size_t i = 0; i < 3; i++) {
			for (size_t j = 0; j < 3; j++) {
				for (size_t k = 0; k < 3; k++) {
					tmp[i][j] += mtx1[i][k] * mtx2[k][j];
				}
			}
		}
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = 0; j < 3; j++) {
				for (size_t k = 0; k < 3; k++) {
					ret[i][j] += tmp[i][k] * mtx3[k][j];
				}
			}
		}

		return ret;
	}

	vector<vector<double>> F2R(double alpha, double beta, double gamma)
	{
		vector<vector<double>> ret = R2F(alpha, beta, gamma);
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = 0; j < 3; j++) {
				if (i != j) {
					double tmp = ret[i][j];
					ret[i][j] = ret[j][i];
					ret[j][i] = tmp;
				}
			}
		}
		return ret;
	}
}

vector<vector<gxx_math::DoubleComplex>> 
Reflectarray::GetIncidentField() { return _incident_field; }

void Reflectarray::AddSource(shared_ptr<Source> src)
{
	_sources.push_back(src);
}

void Reflectarray::ResetSource()
{
	_recompute_incident_field();
}

std::vector<CartesianCS> Reflectarray::GetSourcesPos() const
{
	vector<CartesianCS> res;

	for (auto pSrc : _sources) {
		SourcePosition sp = pSrc->GetPosition();

		vector<double> from{ 0., 0., 0. }, to{ 0., 0., 0. };
		from[2] -= sp.Distance;

		auto f2r = antarray_internal::F2R(sp.Alpha, sp.Beta, sp.Gamma);
		for (size_t m = 0; m < 3; m++) {
			for (size_t n = 0; n < 3; n++) {
				to[m] += f2r[m][n] * from[n];
			}
		}
		res.push_back(CartesianCS(to[0], to[1], to[2]));
	}

	return res;
}

void Reflectarray::_recompute_incident_field()
{
	assert(_sources.size() > 0);
	_incident_field.clear();
	_incident_field.resize(TotalCells());
	for (size_t i = 0; i < TotalCells(); i++)
		_incident_field[i].resize(3);
	for (auto src : _sources) {
		for (size_t i = 0; i < _array_info.size(); i++) {
			vector<double> pt = { _array_info[i].X(),
								  _array_info[i].Y(),
								  0.0 };
			double a = src->GetPosition().Alpha;
			double b = src->GetPosition().Beta;
			double g = src->GetPosition().Gamma;
			double dis = src->GetPosition().Distance;
			vector<double> fpt = { 0., 0., 0. };
			auto r2f = antarray_internal::R2F(a, b, g);
			for (size_t m = 0; m < 3; m++) {
				for (size_t n = 0; n < 3; n++) {
					fpt[m] += r2f[m][n] * pt[n];
				}
			}
			//add z
			fpt[2] += dis;
			auto field = src->RadiationPatternAt(CartesianCS(fpt[0], fpt[1], fpt[2]));
			auto f2r = antarray_internal::F2R(a, b, g);
			vector<gxx_math::DoubleComplex> rField = { gxx_math::DoubleComplex(),
													   gxx_math::DoubleComplex(),
													   gxx_math::DoubleComplex() };
			for (size_t m = 0; m < 3; m++) {
				for (size_t n = 0; n < 3; n++) {
					rField[m] += f2r[m][n] * field[n];
				}
			}

			for (size_t m = 0; m < 3; m++) {
				//_incident_field[i][m] += rField[m];
				_incident_field[i][m] = rField[m];
			}
		}
	}
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

	for (size_t row = 0; row < _yscale; row++) {
		for (size_t col = 0; col < _xscale; col++) {
			ret.push_back(CartesianCS(xlist[col], ylist[row], 0.0));
		}
	}

	return ret;
}



