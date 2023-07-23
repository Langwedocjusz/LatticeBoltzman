#include "Lattice.h"

#include <iostream>
#include <stdexcept>

struct ProgramArgs {
	int NumIterations;
	int SaveStep;

	LatticeSpecification LatticeSpec;
};

ProgramArgs ParseArgs(int argc, char** argv)
{
	const std::vector<std::string> args(argv + 1, argv + argc);

	//Arguments: width, height, num_iterations, save_step

	int iterations = 0, save_step = 1;

	size_t width = 0, height = 0;
	double length_unit = 1.0;
	double time_step   = 1.0;
	double mass_unit   = 1.0;

	try 
	{
		width      = std::stoi(args.at(0));
		height     = std::stoi(args.at(1));
		iterations = std::stoi(args.at(2));
		save_step  = std::stoi(args.at(3));
	}

	catch (const std::out_of_range& ex)
	{
		std::cerr << "Incorrect args\n";
	}

	return ProgramArgs{
		iterations, save_step,
		LatticeSpecification{width, height, length_unit, time_step, mass_unit}
	};
}

int main(int argc, char** argv)
{
	//Parse input
	const auto args = ParseArgs(argc, argv);

	//Initialize lattice
	Lattice lattice(args.LatticeSpec);

	auto dziabdziabdziab = [](size_t idx, size_t idy, Node& node)
	{
		const double sin_x = std::abs(std::sin(static_cast<double>(idx) / 100.0));
		const double sin_y = std::abs(std::sin(static_cast<double>(idy) / 150.0));

		node.Weights[1] = sin_x + 0.001;
		node.Weights[2] = sin_y + 0.001;
	};

	lattice.Initialize(dziabdziabdziab);

	//Setup output directory
	const std::string dir_name{ "output" };

	if (!std::filesystem::is_directory(dir_name) || !std::filesystem::exists(dir_name))
	{
		std::filesystem::create_directory(dir_name);
	}

	std::filesystem::path output_dir(dir_name);

	//Simulate dynamics and save results
	for (int i = 0; i < args.NumIterations; i++)
	{
		//Calculate next time step
		lattice.Update();

		//Once every 'SaveStep' save current configuration to new file
		if (i % args.SaveStep == 0)
		{
			std::string filename = std::to_string(i/args.SaveStep) + ".txt";

			auto filepath = output_dir / filename;

			lattice.Serialize(filepath);
		}
	}

	return 0;
}