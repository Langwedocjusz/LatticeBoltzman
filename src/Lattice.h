#pragma once

#include <array>
#include <vector>
#include <string>
#include <filesystem>

#include "Utils.h"
#include "Scene.h"

#define MULTITHREADED

struct Node {
	//Distribution Weights, initialized to equilibrium distribution:
	std::array<double, 9> Weights{4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	std::array<double, 9> TmpWeights{4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };

	//Macroscopic Quantities:
	double Density{ 1.0 }, InvDensity{1.0};
	Utils::Vec2 Velocity = Utils::Vec2{0.0, 0.0};
	//Solid nodes are excluded from dynamics simulation
	bool IsSolid = false, IsSolidInterior = false;

	//Recalculates macroscopic denisty and velocity
	void UpdateMacroscopic();

	//Calculates equlibrium weight in direction labeled by 'idx'
	//Includes the action of external force 'F', which also requires access to relaxation time 'tau'
	double Equlibrium(uint32_t idx, Utils::Vec2 tauF);

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

		//Relaxation time parameter, for D2Q9 lattice it is related
		//to kinematic viscosity via \nu = (\tau - 0.5)/3.0
		double Tau = 1.0;

		double Gravity;

		std::array<Lattice::BoundaryCondition, 4> BoundaryConditions;

		//We assume components tangent to the boundary are zero
		std::array<double, 4> VonNeumannVelocitiesNormal;

		std::array<double, 4> DirichletDensities;
	};

	Lattice(Specification spec);
	~Lattice();

	//Sets is_solid flag of all nodes according to a pre-defined scene
	void LoadScene(const Scene& scene);
	//Allows initializing node values with an arbitrary function
	void InitFlow(void func(size_t, size_t, Node&));
	//Calculate next step of dynamics simulation
	void Update();
	//Saves macroscopic quantities of all nodes to a text file of given path
	void Serialize(std::filesystem::path filepath);

private:
	void UpdateMacroscopic();
	void StreamingStep();
	void CollisionAndBounceback();
	void HandleBoundaries();
	void HandleBoundary(Boundary boundary);

	Specification m_Spec;

	//1 lattice (length) unit / time step
	double m_BaseSpeed;

	//Relaxation time parameter in units appropriate for including force influence on equilibrium distribution
	double m_TauForce;
	//Inverse of relaxation time in units appropriate for performin relaxation to equilibrium
	double m_Omega;

	std::vector<std::vector<Node>> m_Nodes;

#ifdef MULTITHREADED
	std::vector<size_t> m_Indices;
#endif
};