#include "Utils.h"
#include "gtest/gtest.h"
#include "CoordinateSystem.h"
#include "Source.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <functional>
#include "AntArray.h"

using namespace std;
using namespace gxx_math;

TEST(SPFunc, Fresnel) {
	double pre = 1e-6;
	EXPECT_EQ(DiffLTPrecision(SpecialFunc::Fresnel(0.0), DoubleComplex(), pre), true);
	EXPECT_EQ(DiffLTPrecision(SpecialFunc::Fresnel(-2.0), DoubleComplex(-0.48825340, 0.34341567), pre), true);
	EXPECT_EQ(DiffLTPrecision(SpecialFunc::Fresnel(-5.0), DoubleComplex(-0.56363118, 0.49919137), pre), true);
	EXPECT_EQ(DiffLTPrecision(SpecialFunc::Fresnel(2.0), DoubleComplex(0.48825340, -0.34341567), pre), true);
	EXPECT_EQ(DiffLTPrecision(SpecialFunc::Fresnel(5.0), DoubleComplex(0.56363118, -0.49919137), pre), true);
}

TEST(Utils, SmallFunc) {
	EXPECT_DOUBLE_EQ(RadToDeg(M_PI/3.), 60.);
	EXPECT_DOUBLE_EQ(DegToRad(90), M_PI / 2.0);

	EXPECT_DOUBLE_EQ(dB(10.), 20.);
	DoubleComplex z(1, 2);
	EXPECT_DOUBLE_EQ(dB(z), 20. * log10(Abs(z)));
}

TEST(CS, CSConverter) {
	CartesianCS car{ 1., 2., 3. };
	SphericalCS sph(car);

	EXPECT_DOUBLE_EQ(car.X(), sph.ToCartesian().X());
	EXPECT_DOUBLE_EQ(car.Y(), sph.ToCartesian().Y());
	EXPECT_DOUBLE_EQ(car.Z(), sph.ToCartesian().Z());

	sph.R(5.).Theta(DegToRad(30)).Phi(DegToRad(45));
	car = sph.ToCartesian();

	EXPECT_DOUBLE_EQ(sph.R(), car.ToSpherical().R());
	EXPECT_DOUBLE_EQ(sph.Theta(), car.ToSpherical().Theta());
	EXPECT_DOUBLE_EQ(sph.Phi(), car.ToSpherical().Phi());
}

TEST(Source, Horn) {
	double freq = 5e9;
	double lambda = 3e8 / freq;
	double k0 = 2 * M_PI / lambda;
	PyramidalHorn horn(3.56 * lambda,
		5.08 * lambda,
		0.762 * lambda,
		0.3386 * lambda,
		1.524 * lambda,
		1.1854 * lambda,
		freq, 10);
	vector<DoubleComplex> E_total = horn.RadiationPatternAt(SphericalCS(10, DegToRad(30), DegToRad(25)));

	double pre = 1e-6;
	EXPECT_EQ(DiffLTPrecision(E_total[0], DoubleComplex(0.0, 0.0), pre), true);
	EXPECT_EQ(DiffLTPrecision(E_total[1], DoubleComplex(-3.86551809, -1.66259236), pre), true);
	EXPECT_EQ(DiffLTPrecision(E_total[2], DoubleComplex(-8.28963030, -3.56544083), pre), true);

	//EXPECT_EQ(E_total[1], DoubleComplex(1., 1.));
	//EXPECT_EQ(E_total[2], DoubleComplex(1., 1.));
	EXPECT_EQ(horn.GetTotalRadiationPower() - 22701.26065042 < pre, true);
	//EXPECT_EQ(horn.getTotalPowerRadi(), 1.);
}

TEST(Array, distro) {
}

