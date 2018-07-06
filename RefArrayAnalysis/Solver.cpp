#include "Solver.h"
#include "Utils.h"

#include <fstream>
#include <sstream>

#include <fstream>

using namespace std;
using namespace gxx_math;

void Solver::SetArray(std::shared_ptr<Reflectarray> arr)
{
	ref_array_ = arr;
	each_cell_distribution_.reserve(ref_array_->TotalCells());
	px_ = ref_array_->CellSize();
	py_ = ref_array_->CellSize();
	Nx_ = ref_array_->XScale();
	Ny_ = ref_array_->YScale();
	dx_ = ref_array_->CellSize();
	dy_ = ref_array_->CellSize();
	k0_ = ref_array_->GetSource()->GetWaveNumber0();
}

void Solver::SetCell(const std::string & fn)
{
	load_cell_response_(fn);
}

void Solver::SetPhaseFuzzifier(PhaseFuzzifier fp)
{
	fp_ = fp;
	compute_each_cell_response_();
}

//calc all tasks async
//TODO : add stop flag
void Solver::Run()
{
	for (shared_ptr<Task> pTask : tasks_) {
		pTask->SetStateInfo(Task::Calculating);

		for (size_t i = 0; i < pTask->TaskSize(); i++) {
			SphericalCS pt = (*pTask)[i].ToSpherical();
			double r = pt.R();
			double t = pt.Theta();
			double p = pt.Phi();
			double u = sin(t) * cos(p);
			double v = sin(t) * sin(p);

			pair<DoubleComplex, DoubleComplex> Erxy = erxy_fft_(u, v);

			DoubleComplex E_phi = _I_(-k0_) * Exp(_I_(-k0_*r)) / (2 * M_PI * r) * cos(t) *
				(Erxy.first * sin(p) - Erxy.second * cos(p));
			DoubleComplex E_theta = _I_(k0_) * Exp(_I_(-k0_*r)) / (2 * M_PI * r) *
				(Erxy.first * cos(p) + Erxy.second * sin(p));

			pTask->SetField(i, {DoubleComplex(0., 0.), E_theta, E_phi});
		}

		pTask->SetStateInfo(Task::Finished);
	}
}

// very simple csv parser
void Solver::load_cell_response_(const string& fn)
{
	cell_response_.clear();
	ifstream in(fn);
	if (!in.good()) {
		return;
	}
	string line;
	stringstream ss;
	while (getline(in, line))
	{
		CellResponse cr;
		if (!line.empty()) {
			istringstream iss(line);
			string re, im;
			//TETE
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TETE = DoubleComplex(stod(re), stod(im));

			//TETM
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TETM = DoubleComplex(stod(re), stod(im));

			//TMTE
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TMTE = DoubleComplex(stod(re), stod(im));

			//TMTM
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TMTM = DoubleComplex(stod(re), stod(im));
		}
		cell_response_.push_back(cr);
	}
}

void Solver::compute_each_cell_response_()
{
	each_cell_distribution_.clear();
	if (!fp_ || !distributor_) {
		return;
	}
	vector<double> phase_list;
	for (auto cr : cell_response_) {
		phase_list.push_back(RadToDeg(Arg(cr.TETE)));
	}
	for (auto it = ref_array_->Begin();
		it != ref_array_->End();
		it++) 
	{
		size_t i = fp_(phase_list, (*distributor_)(it->X(), it->Y()));
		each_cell_distribution_.push_back(cell_response_[i]);
	}
}

//core func fft
std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> Solver::erxy_fft_(double u, double v)
{
	DoubleComplex K1 = Exp(_I_(-1 * k0_ / 2.0*(u*(Nx_ - 1)*dx_ + v*(Ny_ - 1)*dy_)));

	DoubleComplex d_sumx;
	DoubleComplex d_sumy;
	for (size_t n = 0; n < Ny_; n++) {
		for (size_t m = 0; m < Nx_; m++) {
			double x = m * px_ - (Nx_ - 1) * px_ / 2.;
			double y = n * py_ - (Ny_ - 1) * py_ / 2.;
			
			pair<DoubleComplex, DoubleComplex> Exymn = emn_(n * Ny_ + m);
			DoubleComplex ejk = Exp(_I_(k0_*(u*m*dx_ + v*n*dy_)));

			d_sumx += (Exymn.first * ejk);
			d_sumy += (Exymn.second * ejk);
		}
	}

	DoubleComplex A = K1 * px_*py_*Sinc(k0_*u*px_ / 2.0)*Sinc(k0_*v*py_ / 2.0);

	return { A*d_sumx, A*d_sumy };
}

//most important func
std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> Solver::emn_(size_t idx)
{
	vector<DoubleComplex> incField = ref_array_->GetIncidentField()[idx];

	DoubleComplex Erx = incField[0];
	DoubleComplex Ery = incField[1];

	//vector<vector<int>> tc = { {0, 1}, {1, 0} };
	vector<DoubleComplex> d12 = { Ery, Erx };

	CellResponse cr = each_cell_distribution_[idx];
	DoubleComplex s11 = cr.TETE;
	DoubleComplex s12 = cr.TETM;
	DoubleComplex s21 = cr.TMTE;
	DoubleComplex s22 = cr.TMTM;

	DoubleComplex a1 = s11 * d12[0] + s12 * d12[1];
	DoubleComplex a2 = s21 * d12[0] + s22 * d12[1];

	return { a2, a1 };
}
