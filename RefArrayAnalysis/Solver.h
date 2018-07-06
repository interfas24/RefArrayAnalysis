#pragma once

#include "Source.h"
#include "AntArray.h"
#include "Utils.h"
#include "Task.h"
#include "BeamDistributor.h"

#include <memory>
#include <list>
#include <vector>
#include <string>


// apply antenna array info and task info
// run to solve : iterate on array
class Solver : public NoCopyable
{
public:
	Solver() 
	{
		ref_array_ = nullptr;
		tasks_.clear();
	}
	
	void SetArray(std::shared_ptr<Reflectarray> arr);
	void SetCell(const std::string& fn);
	void SetPhaseDistributor(std::shared_ptr<BeamDistributor> bd) { distributor_ = bd; }
	void SetPhaseFuzzifier(PhaseFuzzifier fp);
	void AppendTask(std::shared_ptr<Task> t) { tasks_.push_back(t); }
	void Run();

private:
	std::shared_ptr<Reflectarray>		ref_array_;
	std::shared_ptr<BeamDistributor>	distributor_;
	std::vector<CellResponse>			cell_response_;
	//select cell response from cell_response_ for calculation
	std::vector<CellResponse>			each_cell_distribution_;
	PhaseFuzzifier						fp_;
	std::list<std::shared_ptr<Task>>	tasks_;

	double px_, py_;
	size_t Nx_, Ny_;
	double dx_, dy_;
	double k0_;

	void load_cell_response_(const std::string& fn);
	void compute_each_cell_response_();
	std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> erxy_fft_(double u, double v);
	std::pair<gxx_math::DoubleComplex, gxx_math::DoubleComplex> emn_(size_t idx);
};
