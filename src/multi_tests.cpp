#include "Utils.h"
#include "gtest/gtest.h"
#include "CoordinateSystem.h"
#include "Source.h"
#include "AntArray.h"
#include "BeamDistributor.h"
#include "Solver.h"
#include "Task.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <functional>
#include <memory>

using namespace std;
using namespace gxx_math;

#define FILE_OUT_PATH	"./output/"

string GetFilePath(const string &fn)
{
	return string(FILE_OUT_PATH) + fn;
}

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
	string magFN = GetFilePath("inc_mag.txt");
	string phaseFN = GetFilePath("inc_phase.txt");
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

TEST(Phase, PhaseDistro)
{
	int scale = 20;
	double cell_sz = 30.0 / 1000.;
	double freq = 5.0e9;
	double horn_z = cell_sz * scale * 0.8;

	vector<double> yl = Linspace(-scale / 2. * cell_sz + cell_sz / 2.,
								 scale / 2. * cell_sz + cell_sz / 2., scale, false);
	vector<double> xl = yl;

	ofstream out1(GetFilePath("pencil1.txt"));
	ofstream out2(GetFilePath("pencil2.txt"));
	ofstream out3(GetFilePath("oam1.txt"));
	ofstream out4(GetFilePath("oam2.txt"));

	BeamDistributor b1(freq, SourceBeam(0., 0., horn_z),
		PencilBeam(DegToRad(20), DegToRad(0)) &
		PencilBeam(DegToRad(20), DegToRad(180)));

	BeamDistributor b2(freq, SourceBeam(0., 0., horn_z),
		PencilBeam(DegToRad(30), DegToRad(0)));

	BeamDistributor b3(freq, SourceBeam(0., 0., horn_z),
		OAMBeam(DegToRad(0), DegToRad(0), 1));

	BeamDistributor b4(freq, SourceBeam(0., 0., horn_z),
		OAMBeam(DegToRad(0), DegToRad(0), 2));

	for (auto y : yl) {
		for (auto x : xl) {
			out1 << b1(x, y) << ',';
			out2 << b2(x, y) << ',';
			out3 << b3(x, y) << ',';
			out4 << b4(x, y) << ',';
		}
		out1 << '\n';
		out2 << '\n';
		out3 << '\n';
		out4 << '\n';
	}
}

TEST(Solver, pencilbeam)
{
	int scale = 20;
	double cell_sz = 15.0 / 1000.;
	double freq = 10.0e9;

	double lambda = PhysicsConst::LightSpeed / freq;
	double k0 = 2 * M_PI / lambda;
	double E0 = 10.;
	double R = 100.;
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
	double fdr = 1.0;
	double dis = scale * cell_sz * fdr;
	double incAngle = 0;
	pSrc->Place({ DegToRad(180), DegToRad(180 - incAngle), DegToRad(0), dis });

	auto parr = shared_ptr<Reflectarray>(new SquareRefArray(scale, cell_sz));
	parr->AddSource(pSrc);
	parr->ResetSource();

	CartesianCS pos = parr->GetSourcesPos()[0];
	
	auto pbd = shared_ptr<BeamDistributor>(new BeamDistributor(freq, SourceBeam(pos.X(), pos.Y(), pos.Z()),
		PencilBeam(DegToRad(0), DegToRad(0))));
	
	
	shared_ptr<Gain2D> ptask1(new Gain2D(181, DegToRad(0)));
	ptask1->SetInputPower(pSrc->GetTotalRadiationPower());
	shared_ptr<Gain2D> ptask2(new Gain2D(181, DegToRad(90)));
	ptask2->SetInputPower(pSrc->GetTotalRadiationPower());
	shared_ptr<Gain3D> ptask3(new Gain3D(181, 181));
	ptask3->SetInputPower(pSrc->GetTotalRadiationPower());


	Solver s;
	s.SetArray(parr);
	s.SetCell("./input/10G_square_dat.csv");
	s.SetPhaseDistributor(pbd);
	s.SetPhaseFuzzifier(PhaseStepFuzzifier);
	s.AppendTask(ptask1);
	s.AppendTask(ptask2);
	//s.AppendTask(ptask3);
	s.Run();

	vector<double> result = ptask1->PostProcess();
	vector<double> theta = ptask1->GetThetaVec();
	ofstream out(GetFilePath("t0p0.txt"));
	for (size_t i = 0; i < result.size(); i++) {
		out << theta[i] << '\t' << result[i] << '\n';
	}
	out.close();

	vector<double> result2 = ptask2->PostProcess();
	vector<double> theta2 = ptask2->GetThetaVec();
	ofstream out2(GetFilePath("t0p0_90.txt"));
	for (size_t i = 0; i < result2.size(); i++) {
		out2 << theta2[i] << '\t' << result2[i] << '\n';
	}
	out2.close();
	/*
	vector<vector<double>> result3 = ptask3->PostProcess();
	vector<double> t = ptask3->GetThetaVec();
	vector<double> p = ptask3->GetPhiVec();
	ofstream out3(GetFilePath("t0p0_gain3d.txt"));
	for (size_t i = 0; i < p.size(); i++)
	{
		for (size_t j = 0; j < t.size(); j++)
		{
			out3 << result3[i][j] << '\t';
		}
		out3 << '\n';
	}
	out3.close();*/
}

TEST(Solver, oambeam)
{
	int scale = 20;
	double cell_sz = 15.0 / 1000.;
	double freq = 10.0e9;

	double lambda = PhysicsConst::LightSpeed / freq;
	double k0 = 2 * M_PI / lambda;
	double E0 = 10.;
	double R = 100.;
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
	double fdr = 1.0;
	double dis = scale * cell_sz * fdr;
	double incAngle = 0;
	pSrc->Place({ DegToRad(180), DegToRad(180 - incAngle), DegToRad(0), dis });

	auto parr = shared_ptr<Reflectarray>(new SquareRefArray(scale, cell_sz));
	parr->AddSource(pSrc);
	parr->ResetSource();

	CartesianCS pos = parr->GetSourcesPos()[0];

	auto pbd = shared_ptr<BeamDistributor>(new BeamDistributor(freq, SourceBeam(pos.X(), pos.Y(), pos.Z()),
	OAMBeam(DegToRad(30), DegToRad(0), 1)));

	shared_ptr<Gain2D> ptask1(new Gain2D(181, DegToRad(0)));
	ptask1->SetInputPower(pSrc->GetTotalRadiationPower());
	shared_ptr<Gain3D> ptask2(new Gain3D(181, 181));
	ptask2->SetInputPower(pSrc->GetTotalRadiationPower());

	shared_ptr<Observe2D> ptask3(new Observe2D(ObservePlane(DegToRad(30),
		ObservePlane::Y,
		5.0,
		{ 1.5, 1.5 },
		{ 101, 101 })));

	Solver s;
	s.SetArray(parr);
	s.SetCell("./input/10G_square_dat.csv");
	s.SetPhaseDistributor(pbd);
	s.SetPhaseFuzzifier(PhaseStepFuzzifier);
	s.AppendTask(ptask1);
	//s.AppendTask(ptask2);
	//s.AppendTask(ptask3);
	s.Run();

	vector<double> result = ptask1->PostProcess();
	vector<double> theta = ptask1->GetThetaVec();
	ofstream out(GetFilePath("t30p0m1.txt"));
	for (size_t i = 0; i < result.size(); i++) {
		out << theta[i] << '\t' << result[i] << '\n';
	}
	out.close();

	/*
	vector<vector<double>> result2 = ptask2->PostProcess();
	vector<double> t = ptask2->GetThetaVec();
	vector<double> p = ptask2->GetPhiVec();
	ofstream out2(GetFilePath("t30p0m1_gain3d.txt"));
	for (size_t i = 0; i < p.size(); i++)
	{
		for (size_t j = 0; j < t.size(); j++)
		{
			out2 << result2[i][j] << '\t';
		}
		out2 << '\n';
	}
	out2.close();
	*/

	/*
	vector<vector<EOrHField>> result3 = ptask3->PostProcess();
	ofstream out3m(GetFilePath("t30p0m1_mag.txt"));
	ofstream out3p(GetFilePath("t30p0m1_phase.txt"));
	for (size_t i = 0; i < ptask3->GetObservePlane().NXY.first; i++)
	{
		for (size_t j = 0; j < ptask3->GetObservePlane().NXY.second; j++)
		{
			out3m << Abs(result3[i][j].B) << ',';
			out3p << RadToDeg(Arg(result3[i][j].B)) << ',';
		}
		out3m << '\n';
		out3p << '\n';
	}
	out3m.close();
	out3p.close();
	*/
}

int main(int argc, char* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}