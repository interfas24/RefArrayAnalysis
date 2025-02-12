#include "Task.h"

using namespace std;
using namespace gxx_math;

void Task::SetField(size_t i, const EOrHField & field)
{
	field_[i] = field;
}

CartesianCS Task::operator[](size_t i)
{
	return calc_points_[i];
}

std::vector<double> Gain::PostProcess()
{
	vector<double> res;

	//E field R Theta Phi
	for (const EOrHField &field : field_)
	{
		DoubleComplex total = Sqrt(field.B * field.B + field.C * field.C);
		double magGain = Abs(total) * Abs(total) * 4 * M_PI * r_ * r_ / input_power_;
		res.push_back(dB(magGain, 10));
	}

	return res;
}

std::vector<double> Gain2D::PostProcess()
{
	return Gain::PostProcess();
}

std::vector<std::vector<double>> Gain3D::PostProcess()
{
	vector<double> all_res = Gain::PostProcess();
	vector<vector<double>> res;
	res.resize(phi_vec_.size());
	for (size_t i = 0; i < phi_vec_.size(); i++) {
		vector<double> row;
		row.resize(theta_vec_.size());
		for (size_t j = 0; j < theta_vec_.size(); j++) {
			row[j] = all_res[i*phi_vec_.size() + j];
		}
		res[i] = row;
	}
	return res;
}

std::vector<std::vector<EOrHField>> Observe2D::PostProcess()
{
	vector<vector<EOrHField>> res;

	for (size_t i = 0; i < plane_.NXY.first; i++)
	{
		vector<EOrHField> row;
		for (size_t j = 0; j < plane_.NXY.second; j++)
		{
			SphericalCS s = calc_points_[i*plane_.NXY.first+j].ToSpherical();
			EOrHField rtp = field_[i*plane_.NXY.first+j];

			DoubleComplex x = cos(s.Theta())*cos(s.Phi())*rtp.B - sin(s.Phi()) * rtp.C;
			DoubleComplex y = cos(s.Theta())*sin(s.Phi())*rtp.B + cos(s.Phi()) * rtp.C;
			DoubleComplex z = -sin(s.Theta()) * rtp.B;

			EOrHField xyz = { x, y, z };
			row.push_back(xyz);
		}
		res.push_back(row);
	}

	return res;
}
