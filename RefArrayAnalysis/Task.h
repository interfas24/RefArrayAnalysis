#pragma once

#include <vector>
#include "Utils.h"
#include "GxxComplex.h"

class Task : public NoCopyable
{
public:
	Task(size_t sz)
	{
		field_.resize(sz);
	}

	std::vector<EOrHField>::iterator Begin() { return field_.begin(); }
	std::vector<EOrHField>::iterator End() { return field_.end(); }
	size_t TaskSize() const { return field_.size(); }

protected:
	std::vector<EOrHField> field_;	//E field of all target point
};

// theta within [-pi, pi] by default
class Gain2D : public Task
{
public:
	Gain2D(size_t ntheta, double phi, double R = 100.) : Task(ntheta) {}
	std::vector<double> PostProcess();
};

// row : phi, col : theta
class Gain3D : public Task
{
public:
	Gain3D(size_t ntheta, size_t nphi, double R = 100.) : Task(ntheta * nphi) {}
	std::vector<double> PostProcess();
};

class Gain : public Task
{
public:
	Gain(const std::vector<double> &phi,
		 const std::vector<double> &theta,
		 double R = 100.)
		: Task(phi.size() * theta.size()), n_phi_(phi.size()), n_theta_(theta.size())
	{}

	std::vector<double> PostProcess();

private:
	size_t n_phi_;
	size_t n_theta_;
};

// according to : https://en.wikipedia.org/wiki/Rotation_matrix
struct ObservePlane {
	enum Axis { X = 0, Y = 1, Z = 2 };
	ObservePlane(double rotateangle, Axis axis, double dis,
				 std::pair<double, double> wxy, std::pair<size_t, size_t> nxy)
		: RotateAngle(rotateangle), AroundAxis(axis), Distance(dis),
		  PlaneXY(wxy), NXY(nxy)
	{
	}
	double RotateAngle;
	Axis AroundAxis;
	double Distance;
	std::pair<double, double> PlaneXY;
	std::pair<size_t, size_t> NXY;
};

class Observe2D : public Task
{
public:
	Observe2D(const ObservePlane& op)
		: Task(op.NXY.first * op.NXY.second) {}
private:
	//std::vector<double> 
};