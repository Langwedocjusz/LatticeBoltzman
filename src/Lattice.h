#pragma once

#include <array>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.h"

struct Node {
	//Distribution Weights:
	std::array<double, 9> Weights, TmpWeights;
	//Macroscopic Quantities:
	double Density;
	Utils::Vec2 Velocity;
	//Solid nodes are excluded from dynamics simulation
	bool IsSolid = false;

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

struct LatticeSpecification {
	size_t sizeX, sizeY;

	double LengthUnit; //[m]
	double TimeStep;   //[s]
	double MassUnit;   //[kg]
};

class Lattice {
public:
	Lattice(LatticeSpecification spec);
	~Lattice();

	//Calculate next step of dynamics simulation
	void Update();
	//Saves macroscopic quantities of all nodes to a text file of given path
	void Serialize(std::filesystem::path filepath);

private:
	LatticeSpecification m_Spec;

	//1 lattice (length) unit / time step
	double m_BaseSpeed;
	//Viscosity parameter
	double m_Tau = 1.0;

	std::vector<std::vector<Node>> m_Nodes;
};