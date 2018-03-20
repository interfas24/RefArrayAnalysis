#pragma once

#include "Source.h"
#include "AntArray.h"
#include "Utils.h"
#include "Task.h"

#include <memory>
#include <list>
#include <vector>
#include <string>

struct TargetBeam
{
	enum BeamType
	{
		PENCIL_BEAM = 1,
		OAM_BEAM = 2,
		FOCAL_BEAM = 3
	};
};

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
	
	void SetupArray(std::shared_ptr<Reflectarray> arr) { ref_array_ = arr; }
	void SetupCell(const std::string& fn);
	void SetPhaseDistributor();
	void AppendTask(std::shared_ptr<Task> t) { tasks_.push_back(t); }
	void Run();

private:
	std::shared_ptr<Reflectarray>		ref_array_;
	std::list<std::shared_ptr<Task>>	tasks_;
	std::vector<CellResponse>			cell_response_;
	//select cell response from cell_response_ for calculation
	std::vector<CellResponse>			array_cell_distribution_;	
};
