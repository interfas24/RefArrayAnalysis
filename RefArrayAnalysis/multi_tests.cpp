#include "Utils.h"
#include "gtest/gtest.h"
#include "CoordinateSystem.h"
#include "Source.h"
#include <iostream>
#include <vector>
#include <fstream>

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

	EXPECT_DOUBLE_EQ(car.x(), sph.toCartesian().x());
	EXPECT_DOUBLE_EQ(car.y(), sph.toCartesian().y());
	EXPECT_DOUBLE_EQ(car.z(), sph.toCartesian().z());

	sph.r(5.).theta(DegToRad(30)).phi(DegToRad(45));
	car = sph.toCartesian();

	EXPECT_DOUBLE_EQ(sph.r(), car.toSpherical().r());
	EXPECT_DOUBLE_EQ(sph.theta(), car.toSpherical().theta());
	EXPECT_DOUBLE_EQ(sph.phi(), car.toSpherical().phi());
}

TEST(Source, Horn) {
	double freq = 5e9;
	double lambda = 3e8 / freq;
	double k0 = 2 * M_PI / lambda;
	PyramidalHorn horn(3.56, 5.08, 0.762, 0.3386, 1.524, 1.1854, 10, k0);
	vector<DoubleComplex> E_total = horn.RadiPatternAt(SphericalCS(10, DegToRad(30), DegToRad(25)));

	double pre = 1e-6;
	EXPECT_EQ(DiffLTPrecision(E_total[0], DoubleComplex(0.0, 0.0), pre), true);
	EXPECT_EQ(DiffLTPrecision(E_total[1], DoubleComplex(-3.86551809, -1.66259236), pre), true);
	EXPECT_EQ(DiffLTPrecision(E_total[2], DoubleComplex(-8.28963030, -3.56544083), pre), true);

	/*
	freq = 6e9;
	lambda = PhysicsConst::LightSpeed / freq;
	k0 = 2 * M_PI / lambda;
	double E0 = sqrt(120. * M_PI);
	PyramidalHorn h(3.56, 5.08, 0.762, 0.3386, 1.524, 1.1854, E0, k0);
	vector<double> thetas = Linspace(-90, 90, 91);
	vector<double> Ertp;
	double phi = 90;
	E_total.reserve(thetas.size());
	for (const double & t : thetas) {
		SphericalCS s(15, DegToRad(t), DegToRad(phi));
		vector<DoubleComplex> re_xyz = h.RadiPatternAt(s.toCartesian());
		DoubleComplex total = Sqrt((re_xyz[0] ^ 2.) + (re_xyz[1] ^ 2.) + (re_xyz[2] ^ 2.));
		//Ertp.push_back(dB(total));
		Ertp.push_back(dB(re_xyz[2]));
	}
	
	ofstream out("rEz_90.txt");
	for (size_t i = 0; i < thetas.size(); i++) {
		out << thetas[i] << "\t" << Ertp[i] << endl;
	}
	out.close();
	*/
}
