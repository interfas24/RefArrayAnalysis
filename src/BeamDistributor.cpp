#include "BeamDistributor.h"

using namespace gxx_math;

double ScalarSuperposition::AddAll(double k0, double x, double y) const
{
	double ret = 0.0;
	for (auto p : beams_) 
		ret += p->Eval(k0, x, y);
	return ret;
}

double VectorSuperposition::AddAll(double k0, double x, double y) const
{
	DoubleComplex ret;
	for (auto p : beams_)
		ret += Exp(_I_(1.0) * p->Eval(k0, x, y));
	double rad = Arg(ret);
	return rad < 0 ? rad + 2*M_PI : rad;
}

double BeamDistributor::operator()(double x, double y)
{
	double ret = feeds_->AddAll(k0_, x, y) + beams_->AddAll(k0_, x, y);
	ret = RadToDeg(ret);
	return static_cast<double>(static_cast<int>(ret) % 360);
}

double PencilBeam::Eval(double k0, double x, double y)
{
	return -k0 * sin(theta_) * cos(phi_) * x - k0 * sin(theta_) * sin(phi_) * y;
}

double SourceBeam::Eval(double k0, double x, double y)
{
	return k0 * sqrt((xf - x)*(xf - x) + (yf - y)*(yf - y) + zf*zf);
}

double OAMBeam::Eval(double k0, double x, double y)
{
	double phi_n = -k0 * sin(theta_) * cos(phi_) * x - k0 * sin(theta_) * sin(phi_) * y;
	double phi_k = atan2(-cos(theta_)*sin(phi_)*x + cos(theta_)*cos(phi_)*y, cos(phi_)*x + sin(phi_)*y);
	return phi_n + phi_k * mode_;
}
