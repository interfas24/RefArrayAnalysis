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
}

void Solver::SetCell(const std::string & fn)
{
	load_cell_response_(fn);
}

void Solver::SetPhaseFuzzifier(PhaseFuzzifier fp)
{
	fp_ = fp;
	compute_each_cell_response_();

	ofstream out1("RealPhase.txt");
	for (size_t i = 0; i < each_cell_distribution_.size(); i++) {
		if (i % ref_array_->XScale() == 0)
			out1 << '\n';
		out1 << RadToDeg(Arg(each_cell_distribution_[i].TETE)) << ',';
	}
}

void Solver::Run()
{

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

			//TETE
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TETM = DoubleComplex(stod(re), stod(im));

			//TETE
			getline(iss, re, ',');
			getline(iss, im, ',');
			cr.TMTE = DoubleComplex(stod(re), stod(im));

			//TETE
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
