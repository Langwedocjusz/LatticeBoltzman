#pragma once

#include <array>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.h"
#include "Scene.h"

struct Node {
	//Distribution Weights:
	std::array<double, 9> Weights{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::array<double, 9> TmpWeights{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	//Macroscopic Quantities:
	double Density{ 0.0 };
	Utils::Vec2 Velocity = Utils::Vec2{0.0, 0.0};
	//Solid nodes are excluded from dynamics simulation
	bool IsSolid = false, IsSolidInterior;

	//Recalculates macroscopic denisty and velocity
	void UpdateMacroscopic();

	//Calculates equlibrium weight in direction labeled by 'idx'
	double Equlibrium(double base_speed, uint32_t idx);

	//Velocities in 9 base directions (zero, four axis aligned, four diagonal)
	static constexpr std::array<Utils::Vec2, 9> s_BaseVelocities{
		//Uses following convention from "Lattice Boltzman Modelling 
		// An Introduction for Geoscientists and Engineers":
		//   6 -- 2 -- 5
		//   |    |    |
		//   3 -- 0 -- 1
		//   |    |    |
		//   7 -- 4 -- 8
		Utils::Vec2{0.0, 0.0},
		Utils::Vec2{1.0, 0.0}, Utils::Vec2{ 0.0, 1.0}, Utils::Vec2{-1.0,  0.0}, Utils::Vec2{0.0, -1.0},
		Utils::Vec2{1.0, 1.0}, Utils::Vec2{-1.0, 1.0}, Utils::Vec2{-1.0, -1.0}, Utils::Vec2{1.0, -1.0}
	};
};

class Lattice {
public:

	enum Boundary {
		Right = 0, Up = 1, Left = 2, Down = 3
	};

	enum class BoundaryCondition {
		Periodic, VonNeumann, Dirichlet
	};

	struct Specification {
		size_t sizeX, sizeY;

		double LengthUnit; //[m]
		double TimeStep;   //[s]
		double MassUnit;   //[kg]

		//BoundaryCondition UpBC, DownBC, RightBC, LeftBC;
		std::array<Lattice::BoundaryCondition, 4> BoundaryConditions;

		//We assume components tangent to the boundary are zero
		std::array<double, 4> VonNeumannVelocitiesNormal;
	};

	Lattice(Specification spec);
	~Lattice();

	//Sets is_solid flag of all nodes according to a pre-defined scene
	void LoadScene(const Scene& scene);
	//Applies func(id_x, id_y, Node&) on all non-solid nodes, to initialize them
	void InitFlow(void func(size_t, size_t, Node&));
	//Calculate next step of dynamics simulation
	void Update();
	//Saves macroscopic quantities of all nodes to a text file of given path
	void Serialize(std::filesystem::path filepath);

private:
	void CalculateVonNeumann(Boundary boundary);

	Specification m_Spec;

	//Viscosity parameter in lu^2 / ts units
	double m_Tau = 1.0;
	//1 lattice (length) unit / time step (initialized in the constructor)
	double m_BaseSpeed;

	std::vector<std::vector<Node>> m_Nodes;
};