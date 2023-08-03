#include "Lattice.h"

#include <iostream>
#include <stdexcept>

struct ProgramArgs {
	int NumIterations;
	int SaveStep;

	Lattice::Specification LatticeSpec;
};

ProgramArgs ParseArgs(int argc, char** argv)
{
	const std::vector<std::string> args(argv + 1, argv + argc);

	//Arguments: width, height, num_iterations, save_step

	int iterations = 0, save_step = 1;

	size_t width = 0, height = 0;
	double length_unit = 1.0;
	double time_step   = 0.02;
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
		std::cerr << "Supported arguments are (in order): <width> <height> <num iterations> <save step>\n";
	}

	auto boundary = Lattice::BoundaryCondition::VonNeumann;
	auto periodic = Lattice::BoundaryCondition::Periodic;

	return ProgramArgs{
		iterations, save_step,

		Lattice::Specification{
			width, height, 
			length_unit, time_step, mass_unit,
			{boundary, periodic, boundary, periodic},
			{0.09, 0.0, 0.09, 0.0},
			{0.5, 2.0, 1.5, 2.0}
		}
	};
}

int main(int argc, char** argv)
{
	//Parse input
	const auto args = ParseArgs(argc, argv);

	//Initialize lattice
	Lattice lattice(args.LatticeSpec);

	Scene scene;
	scene.PushShape<Circle>(Utils::Vec2{ 64.0, 64.0 }, 25.0);

	lattice.LoadScene(scene);

	auto dziabdziabdziab = [](size_t idx, size_t idy, Node& node)
	{
		const double sin_x = std::abs(std::sin(static_cast<double>(idx) / 100.0));
		const double sin_y = std::abs(std::sin(static_cast<double>(idy) / 150.0));

		node.Weights[0] = 0.1 + 0.9*static_cast<double>(idx)/128.0;

		//node.Weights[1] = sin_x + 0.001;
		//node.Weights[2] = sin_y + 0.001;
	};

	lattice.InitFlow(dziabdziabdziab);

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