#include "Utils.h"
#include "gtest/gtest.h"
#include "CoordinateSystem.h"
#include "Source.h"
#include "AntArray.h"
#include "BeamDistributor.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <functional>
#include <memory>

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

TEST(Sources, Horn) {
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
	vector<DoubleComplex> E_total = horn.RadiationPatternAt(CartesianCS(0.225, -0.225, 0.4).ToSpherical());

	EXPECT_DOUBLE_EQ(E_total[1].real(), 21.333821518683511);
	EXPECT_DOUBLE_EQ(E_total[1].imag(), 89.307789828747246);
	
	EXPECT_DOUBLE_EQ(E_total[2].real(), -21.333821518683511);
	EXPECT_DOUBLE_EQ(E_total[2].imag(), -89.307789828747246);

	EXPECT_DOUBLE_EQ(horn.GetTotalRadiationPower(), 22764.839057946094);
}

TEST(Array, Distro) {
	double freq = 5e9;
	double lambda = 3e8 / freq;
	double k0 = 2 * M_PI / lambda;
	shared_ptr<Source> pSrc(
			new PyramidalHorn(
				3.56 * lambda,
				5.08 * lambda,
				0.762 * lambda,
				0.3386 * lambda,
				1.524 * lambda,
				1.1854 * lambda,
				freq, 10
			));
	pSrc->Place({ DegToRad(180), DegToRad(180), DegToRad(0), 1.0 });
	Reflectarray ref_arr(20, 20, 30 / 1000.);
	ref_arr.AddSource(pSrc);
	auto pt = ref_arr.Tests();
	EXPECT_DOUBLE_EQ(-1.0, pt[0]);
	EXPECT_DOUBLE_EQ(2.0, pt[1]);
	EXPECT_DOUBLE_EQ(-3.0, pt[2]);
}

/*
TEST(Array, IncidentField) {
	double freq = 5e9;
	double lambda = PhysicsConst::LightSpeed / freq;
	double k0 = 2 * M_PI / lambda;
	double E0 = 10.;
	shared_ptr<Source> pSrc(
			new PyramidalHorn(
				3.56 * lambda,
				5.08 * lambda,
				0.762 * lambda,
				0.3386 * lambda,
				1.524 * lambda,
				1.1854 * lambda,
				freq, E0
			));
	shared_ptr<Source> pSrc2(
			new PyramidalHorn(
				3.56 * lambda,
				5.08 * lambda,
				0.762 * lambda,
				0.3386 * lambda,
				1.524 * lambda,
				1.1854 * lambda,
				freq, E0
			));
	size_t scale = 60;
	double cell_sz = 10. / 1000.;
	double fdr = 0.8;
	double dis = scale * cell_sz * fdr;
	pSrc->Place({ DegToRad(180), DegToRad(150), DegToRad(0), dis });
	pSrc2->Place({ DegToRad(180), DegToRad(210), DegToRad(0), dis });

	SquareRefArray arr(scale, cell_sz);
	arr.AddSource(pSrc);
	arr.AddSource(pSrc2);
	arr.ResetSource();

	vector<vector<DoubleComplex>> field = arr.GetIncidentField();
	auto distro = arr.GetArrayDistro();
	string magFN = "mag2.txt";
	string phaseFN = "phase2.txt";
	ofstream out1(magFN);
	ofstream out2(phaseFN);
	for (size_t i = 0; i < scale; i++) {
		for (size_t j = 0; j < scale; j++) {
			DoubleComplex cplxMag = Hypot(field[i*scale+j][0], field[i*scale+j][1], field[i*scale+j][2]);
			double phase = Arg(field[i*scale+j][1]);
			out1 << Abs(cplxMag) << ",";
			out2 << phase << ",";
		}
		out1 << "\n";
		out2 << "\n";
	}
	out1.close();
	out2.close();
}
*/
TEST(Phase, PhaseDistro)
{
	int scale = 20;
	double cell_sz = 30.0 / 1000.;
	double freq = 5.0e9;
	double horn_z = cell_sz * scale * 0.8;

	vector<double> yl = Linspace(-scale / 2. * cell_sz + cell_sz / 2.,
								 scale / 2. * cell_sz + cell_sz / 2., scale, false);
	vector<double> xl = yl;

	shared_ptr<VectorSuperposition> feed(new VectorSuperposition(SourceBeam(0., 0., horn_z)));
	shared_ptr<VectorSuperposition> beam(new VectorSuperposition(PencilBeam(DegToRad(30), DegToRad(0))));

	ofstream out1("pencil1.txt");
	for (auto y : yl) {
		for (auto x : xl) {
			//BeamDistributor b(freq, x, y, feed, beam);
			BeamDistributor b(freq, x, y, feed,
								PencilBeam(DegToRad(20), DegToRad(0)) &
								PencilBeam(DegToRad(20), DegToRad(180)));
			out1 << b() << ',';
		}
		out1 << '\n';
	}
}