#pragma once

#include "Utils.h"
#include <vector>
#include <memory>

class ScalarSuperposition;
class VectorSuperposition;
class Superposition;
class BeamFormer;
class SourceBeam;
class PencilBeam;
class OAMBeam;
class FocalBeam;


class BeamFormer
{
public:
	virtual double Eval(double k0, double x, double y) = 0;
};

class SourceBeam : public BeamFormer
{
public:
	SourceBeam(double x, double y, double z)
		: xf(x), yf(y), zf(z) {}
	double Eval(double k0, double x, double y) override;

private:
	double xf;
	double yf;
	double zf;
};

class PencilBeam : public BeamFormer
{
public:
	PencilBeam(double t, double p)
		: theta_(t), phi_(p) {}
	double Eval(double k0, double x, double y) override;

private:
	double theta_;
	double phi_;
};

class OAMBeam : public BeamFormer
{
public:
	OAMBeam(double theta, double phi, double mode)
		: theta_(theta), phi_(phi), mode_(mode) {}
	double Eval(double k0, double x, double y) override;

private:
	double theta_;
	double phi_;
	double mode_;
};

class FocalBeam : public BeamFormer
{
public:
	FocalBeam(double x, double y, double z);
	double Eval(double k0, double x, double y) override;
};

class Superposition : public NoCopyable
{
public:
	virtual double AddAll(double k0, double x, double y) const = 0;
	void AddBeam(std::shared_ptr<BeamFormer> p) { beams_.push_back(p); }
protected:
	std::vector<std::shared_ptr<BeamFormer>> beams_;
};

class ScalarSuperposition : public Superposition
{
public:
	double AddAll(double k0, double x, double y) const override;
};

class VectorSuperposition : public Superposition
{
public:
	VectorSuperposition() = default;
	VectorSuperposition(const PencilBeam& beam)
	{
		beams_.push_back(std::shared_ptr<BeamFormer>(new PencilBeam(beam)));
	}
	VectorSuperposition(const OAMBeam& beam)
	{
		beams_.push_back(std::shared_ptr<BeamFormer>(new OAMBeam(beam)));
	}
	VectorSuperposition(const SourceBeam& beam)
	{
		beams_.push_back(std::shared_ptr<BeamFormer>(new SourceBeam(beam)));
	}
	double AddAll(double k0, double x, double y) const override;
};

class BeamDistributor : public NoCopyable
{
public:
	BeamDistributor(double freq, 
		std::shared_ptr<Superposition> feeds,
		std::shared_ptr<Superposition> beams)
		: feeds_(feeds), beams_(beams) 
	{
		double lmbd = PhysicsConst::LightSpeed / freq;
		k0_ = 2 * M_PI / lmbd;
	}
	BeamDistributor(double freq,
		const SourceBeam& feed,
		std::shared_ptr<Superposition> beams)
		: BeamDistributor(freq,
			std::shared_ptr<Superposition>(new VectorSuperposition(feed)), beams)
	{
	}

	BeamDistributor(double freq,
		const SourceBeam& feed, const PencilBeam& beam)
		: BeamDistributor(freq,
			std::shared_ptr<Superposition>(new VectorSuperposition(feed)),
			std::shared_ptr<Superposition>(new VectorSuperposition(beam)))
	{
	}
	BeamDistributor(double freq,
		std::shared_ptr<Superposition> feeds, const PencilBeam& beam)
		: BeamDistributor(freq, feeds,
			std::shared_ptr<Superposition>(new VectorSuperposition(beam)))
	{
	}

	BeamDistributor(double freq,
		const SourceBeam& feed, const OAMBeam& beam)
		: BeamDistributor(freq,
			std::shared_ptr<Superposition>(new VectorSuperposition(feed)),
			std::shared_ptr<Superposition>(new VectorSuperposition(beam)))
	{
	}
	BeamDistributor(double freq,
		std::shared_ptr<Superposition> feeds, const OAMBeam& beam)
		: BeamDistributor(freq, feeds,
			std::shared_ptr<Superposition>(new VectorSuperposition(beam)))
	{
	}

	double operator()(double x, double y);

private:
	double k0_;
	std::shared_ptr<Superposition> feeds_;
	std::shared_ptr<Superposition> beams_;
};

inline
std::shared_ptr<VectorSuperposition>
operator&(std::shared_ptr<VectorSuperposition> vs, const PencilBeam &beam)
{
	std::shared_ptr<BeamFormer> p(new PencilBeam(beam));
	vs->AddBeam(p);
	return vs;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(std::shared_ptr<VectorSuperposition> vs, const OAMBeam &beam)
{
	std::shared_ptr<BeamFormer> p(new OAMBeam(beam));
	vs->AddBeam(p);
	return vs;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(std::shared_ptr<VectorSuperposition> vs, const FocalBeam &beam)
{
	std::shared_ptr<BeamFormer> p(new FocalBeam(beam));
	vs->AddBeam(p);
	return vs;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(std::shared_ptr<VectorSuperposition> vs, const SourceBeam &beam)
{
	std::shared_ptr<BeamFormer> p(new SourceBeam(beam));
	vs->AddBeam(p);
	return vs;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(const PencilBeam& lhs, const PencilBeam& rhs)
{
	std::shared_ptr<VectorSuperposition> p(new VectorSuperposition);
	return (p & lhs & rhs);
}
inline
std::shared_ptr<VectorSuperposition>
operator&(const PencilBeam& pb, const OAMBeam& ob)
{
	std::shared_ptr<VectorSuperposition> p(new VectorSuperposition);
	return p & pb & ob;
}
inline
std::shared_ptr<VectorSuperposition>
operator&( const OAMBeam& ob, const PencilBeam& pb)
{
	return pb & ob;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(const PencilBeam& pb, const FocalBeam& fb)
{
	std::shared_ptr<VectorSuperposition> p(new VectorSuperposition);
	return p & pb & fb;
}
inline
std::shared_ptr<VectorSuperposition>
operator&( const FocalBeam& fb, const PencilBeam& pb)
{
	return pb & fb;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(const FocalBeam& fb, const OAMBeam& ob)
{
	std::shared_ptr<VectorSuperposition> p(new VectorSuperposition);
	return p & fb & ob;
}
inline
std::shared_ptr<VectorSuperposition>
operator&( const OAMBeam& ob, const FocalBeam& fb)
{
	return fb & ob;
}

inline
std::shared_ptr<VectorSuperposition>
operator&(const SourceBeam& lhs, const SourceBeam& rhs)
{
	std::shared_ptr<VectorSuperposition> p(new VectorSuperposition);
	return p & lhs & rhs;
}