#include "Lattice.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>

#include "nlohmann/json.hpp"

struct ProgramArgs {
	int NumIterations = 0;
	int SaveStep = 1;

	Lattice::Specification LatticeSpec;

	std::vector<AABB> AABBs;
	std::vector<Circle> Circles;
};

ProgramArgs ParseArgs(int argc, char** argv);

int main(int argc, char** argv)
{
	//Parse input
	ProgramArgs args;

	try
	{
		args = ParseArgs(argc, argv);
	}
	
	catch (std::invalid_argument ex)
	{
		std::cerr << ex.what() << '\n';
		std::cerr << "Only supported argument is a path to simulation config file.\n";
		return -1;
	}

	catch (std::runtime_error ex)
	{
		std::cerr << ex.what() << '\n';
		return -1;
	}

	//Initialize lattice
	Lattice lattice(args.LatticeSpec);

	Scene scene;
	
	for (const auto& aabb : args.AABBs)
	{
		scene.PushShape<AABB>(aabb);
	}

	for (const auto& circle : args.Circles)
	{
		scene.PushShape<Circle>(circle);
	}

	lattice.LoadScene(scene);

	//Setup output directory
	const std::string dir_name{ "output" };

	if (!std::filesystem::is_directory(dir_name) || !std::filesystem::exists(dir_name))
	{
		std::filesystem::create_directory(dir_name);
	}

	std::filesystem::path output_dir(dir_name);

	try
	{
		//Save current time
		auto start_time = std::chrono::high_resolution_clock::now();

		//Simulate dynamics and save results
		for (int i = 0; i < args.NumIterations; i++)
		{
			//Calculate next time step
			lattice.Update();

			//Once every 'SaveStep' save current configuration to new file
			if (i % args.SaveStep == 0)
			{
				std::string filename = std::to_string(i / args.SaveStep) + ".txt";

				auto filepath = output_dir / filename;

				lattice.Serialize(filepath);
			}
		}

		//Compute and print time it took to finish simulation
		auto end_time = std::chrono::high_resolution_clock::now();

		double time_seconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() * 1e-9;

		std::cout << "Simulation (" << args.NumIterations << " steps) took " << time_seconds << " seconds.\n";
	}

	catch(std::runtime_error ex)
	{
		std::cerr << ex.what() << '\n';
	}

	return 0;
}

ProgramArgs ParseArgs(int argc, char** argv)
{
	//Check if there is only one argument
	const std::vector<std::string> args(argv + 1, argv + argc);

	if (args.size() != 1)
	{
		throw std::invalid_argument("Incorrect arguments.");
	}

	//Try to open file, treating argument as a path
	auto current_path = std::filesystem::current_path();
	auto filepath = current_path / args[0];

	std::ifstream input(filepath);

	if (!input)
	{
		throw std::runtime_error("Could not open file: " + filepath.string());
	}

	//Parse the opened file as a json
	auto json = nlohmann::json::parse(input);

	ProgramArgs ret;

	ret.LatticeSpec.sizeX = json["LatticeSizeX"];
	ret.LatticeSpec.sizeY = json["LatticeSizeY"];

	ret.NumIterations = json["NumIterations"];
	ret.SaveStep = json["SaveStep"];

	ret.LatticeSpec.LengthUnit = json["LengthUnit"];
	ret.LatticeSpec.TimeStep = json["TimeStep"];
	ret.LatticeSpec.MassUnit = json["MassUnit"];

	ret.LatticeSpec.Tau = json["Tau"];

	ret.LatticeSpec.Gravity = json["Gravity"];

	std::map <std::string, Lattice::BoundaryCondition> DictionaryBC{
		{"Periodic"  , Lattice::BoundaryCondition::Periodic},
		{"VonNeumann", Lattice::BoundaryCondition::VonNeumann},
		{"Dirichlet" , Lattice::BoundaryCondition::Dirichlet},
	};

	for (int i = 0; i < 4; i++)
	{
		ret.LatticeSpec.BoundaryConditions[i] = DictionaryBC.at(json["BoundaryConditions"][i]);
		ret.LatticeSpec.VonNeumannVelocitiesNormal[i] = json["VonNeumannVelocitiesNormal"][i];
		ret.LatticeSpec.DirichletDensities[i] = json["DirichletDensities"][i];
	}

	auto aabb_data = json["AABBs"];

	for (size_t i = 0; i < aabb_data.size(); i++)
	{
		Utils::Vec2 min{
			aabb_data[i]["Min"][0],
			aabb_data[i]["Min"][1],
		};

		Utils::Vec2 max{
			aabb_data[i]["Max"][0],
			aabb_data[i]["Max"][1],
		};

		ret.AABBs.push_back(
			AABB(min, max)
		);
	}

	auto circles_data = json["Circles"];

	for (size_t i = 0; i < circles_data.size(); i++)
	{
		Utils::Vec2 pos{
			circles_data[i]["Position"][0],
			circles_data[i]["Position"][1],
		};
		
		double rad = circles_data[i]["Radius"];

		ret.Circles.push_back(
			Circle(pos, rad)
		);
	}

	return ret;
}