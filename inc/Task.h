#pragma once

#include <vector>
#include "Utils.h"
#include "GxxComplex.h"
#include "CoordinateSystem.h"

class Task : public NoCopyable
{
public:

	enum StateInfo
	{
		Pending = 0,
		Calculating = 1,
		Finished = 2
	};

	Task(size_t sz)
	{
		field_.resize(sz);
		calc_points_.resize(sz);
		state_info_ = Pending;
	}

	std::vector<EOrHField>::iterator Begin() { return field_.begin(); }
	std::vector<EOrHField>::iterator End() { return field_.end(); }

	size_t TaskSize() const { return field_.size(); }
	void SetField(size_t i, const EOrHField& field);
	CartesianCS operator[](size_t i);

	void SetStateInfo(StateInfo si) { state_info_ = si; }
	StateInfo GetStateInfo() const { return state_info_; }

protected:
	std::vector<EOrHField> field_;	//E or H field of all target point
	std::vector<CartesianCS> calc_points_;
	StateInfo state_info_;
};

class Gain : public Task
{
public:
	Gain(const std::vector<double> &phi,
		 const std::vector<double> &theta,
		 double R = 100.)
		: Task(phi.size() * theta.size()),
		  phi_vec_(phi),
		  theta_vec_(theta),
		  input_power_(0.0),
		  r_(R)
	{
		for (size_t i = 0; i < phi.size(); i++) {
			for (size_t j = 0; j < theta.size(); j++) {
				SphericalCS pt(r_, theta[j], phi[i]);
				calc_points_[i*phi.size() + j] = pt.ToCartesian();
			}
		}
	}

	std::vector<double> PostProcess();
	void SetInputPower(double ip) { input_power_ = ip; }
	std::vector<double> GetThetaVec() const { return theta_vec_; }
	std::vector<double> GetPhiVec() const { return phi_vec_; }

protected:
	std::vector<double> phi_vec_;
	std::vector<double> theta_vec_;
	double r_;
	double input_power_;
};

// theta within [-pi/2, pi/2] by default
class Gain2D : public Gain
{
public:
	Gain2D(size_t ntheta, double phi, double R = 100.) 
		: Gain({phi}, Linspace(-M_PI / 2, M_PI / 2, ntheta)) {}
	std::vector<double> PostProcess();
};

// row : phi, col : theta
class Gain3D : public Gain
{
public:
	Gain3D(size_t ntheta, size_t nphi, double R = 100.) 
		: Gain(Linspace(0, 2*M_PI, nphi), Linspace(-M_PI / 2, M_PI / 2, ntheta)) {}
	std::vector<std::vector<double>> PostProcess();
};

// according to : https://en.wikipedia.org/wiki/Rotation_matrix
struct ObservePlane {
	enum Axis { X = 0, Y = 1, Z = 2, NONE = 3 };

	//rotate angle : rad
	ObservePlane(double rotateangle, Axis axis, double dis,
				 std::pair<double, double> wxy, std::pair<size_t, size_t> nxy)
		: RotateAngle(rotateangle), AroundAxis(axis), Distance(dis),
		  PlaneXY(wxy), NXY(nxy)
	{
		std::vector<double> xlist = Linspace(-wxy.first / 2, wxy.first / 2, nxy.first);
		std::vector<double> ylist = Linspace(-wxy.second / 2, wxy.first / 2, nxy.second);
		std::vector<std::vector<double>> conv_mtx;
		if (axis == X) 
		{
			conv_mtx = {
				{ 1, 0, 0 },
				{ 0, cos(RotateAngle), -sin(RotateAngle) },
				{ 0, sin(RotateAngle), cos(RotateAngle) }
			};
		}
		else if(axis == Y)
		{
			conv_mtx = {
				{ cos(RotateAngle), 0, sin(RotateAngle) },
				{ 0, 1, 0 },
				{ -sin(RotateAngle), 0, cos(RotateAngle) }
			};
		}
		else if(axis == Z)
		{
			conv_mtx = {
				{ cos(RotateAngle), -sin(RotateAngle), 0 },
				{ sin(RotateAngle), cos(RotateAngle), 0 },
				{ 0, 0, 1 }
			};
		}
		else
		{
			conv_mtx = {
				{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 }
			};
		}

		for (double x : xlist) {
			for (double y : ylist) {
				std::vector<double> pt = { x, y, Distance };
				std::vector<double> new_pt{ 0., 0., 0. };
				for (size_t i = 0; i < 3; i++) {
					new_pt[i] = conv_mtx[i][0] * pt[0] +
								conv_mtx[i][1] * pt[1] +
								conv_mtx[i][2] * pt[2];
				}

				AllPoints.push_back(CartesianCS(new_pt[0], new_pt[1], new_pt[2]));
			}
		}
	}
	double RotateAngle;
	Axis AroundAxis;
	double Distance;
	std::pair<double, double> PlaneXY;
	std::pair<size_t, size_t> NXY;
	std::vector<CartesianCS> AllPoints;
};

class Observe2D : public Task
{
public:
	Observe2D(const ObservePlane& op)
		: Task(op.NXY.first * op.NXY.second), plane_(op)
	{
		calc_points_ = op.AllPoints;
	}

	//x y z field
	std::vector<std::vector<EOrHField>> PostProcess();

	ObservePlane GetObservePlane() const { return plane_; }

private:
	ObservePlane plane_;
};