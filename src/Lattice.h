#pragma once

#include <array>
#include <vector>
#include <cstdint>

#include "Utils.h"

struct LatticeSpecification {
	uint32_t sizeX, sizeY;
	
	double lengthUnit; //[m]
	double timeStep;   //[s]
	double massUnit;   //[kg]
};

struct Node {
	//Distribution Weights:
	std::array<double, 9> weights;
	//Macroscopic Quantities:
	double density;
	Utils::Vec2 velocity;

	//Recalculates macroscopic denisty and velocity
	void UpdateMacroscopic();

	//Calculates equlibrium weight in direction 'idx'
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
	Lattice(LatticeSpecification spec);
	~Lattice();

	//Calculate next step of dynamics simulation
	void Update();

private:
	LatticeSpecification m_Spec;

	//1 lattice (length) unit / time step
	double m_BaseSpeed;

	std::vector<std::vector<Node>> m_Nodes;
};