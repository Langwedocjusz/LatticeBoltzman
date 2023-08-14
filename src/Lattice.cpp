#include "Lattice.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <execution>

using namespace Utils;

void Node::UpdateMacroscopic()
{
	//Macroscopic density is the sum of all weights
	Density = 0.0;

	for (const auto& weight : Weights)
	{
		Density += weight;
	}

	//Do not attempt to calculate velocity if density is zero, as that is ill defined
	if (Density <= 0.0)
	{
		throw std::runtime_error("Numerical error: density reached non-positive value.");
	}

	//Purely for optimization purposes:
	InvDensity = 1.0 / Density;

	//Macroscopic velocity is the weighted average of base velocities
	Velocity = Vec2{ 0.0, 0.0 };

	for (size_t i = 0; i < Weights.size(); i++)
	{
		Velocity.x += Weights[i] * s_BaseVelocities[i].x;
		Velocity.y += Weights[i] * s_BaseVelocities[i].y;
	}

	Velocity.x *= InvDensity;
	Velocity.y *= InvDensity;
}

double Node::Equlibrium(uint32_t idx, Utils::Vec2 tauF)
{
	//Formula for local equilibrium weights as seen in equation (17) in
	//"Lattice Boltzman Modelling An Introduction for Geoscientists and Engineers"
	//It is actually a truncated Taylor expansion of Maxwell distribution

	constexpr std::array<double, 9> weights{
		4.0 / 9.0,
		1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
	};

	const Vec2 e = s_BaseVelocities[idx];

	//Adjusting velocity used in equilibrium calculation to include forces
	//Care needs to be taken to insure that tau used here is in units corresponding
	//to our base velocities normalization
	const Vec2 u = Velocity + InvDensity*tauF;

	const double eDotU = dot(e, u);
	const double u2 = dot(u, u);

	//Formula in the book contained divisions by characteristic velocity c
	//They can be avoided by normalizing base velocities to lengths {1, \sqrt{2}}
	//Instead of {c, \sqrt{2} c}:

	return weights[idx] * Density *
		(1.0 + 3.0 * eDotU + 4.5 * (eDotU*eDotU) - 1.5 * u2);
}

Lattice::Lattice(Specification spec)
	: m_Spec(spec), m_BaseSpeed(spec.LengthUnit/spec.TimeStep)
{
	m_Nodes.resize(m_Spec.sizeX, std::vector<Node>(m_Spec.sizeY));

#ifdef MULTITHREADED
	
	m_Indices.reserve(m_Spec.sizeX);

	for (size_t i = 0; i < m_Spec.sizeX; i++)
	{
		m_Indices.push_back(i);
	}

#endif

	//Convert tau from lu^2 / dt units:
	const double tau = m_Spec.Tau * m_Spec.LengthUnit * m_Spec.LengthUnit / m_Spec.TimeStep;
	//This is the same as above, divided by base velocity
	m_TauForce = m_Spec.Tau * m_Spec.LengthUnit;
	m_Omega = 1.0 / tau;
}

Lattice::~Lattice()
{

}

void Lattice::Update()
{
	UpdateMacroscopic();
	StreamingStep();
	HandleBoundaries();
	CollisionAndBounceback();
}

void Lattice::UpdateMacroscopic()
{
#ifdef MULTITHREADED
	std::for_each(std::execution::par, m_Indices.begin(), m_Indices.end(), [this](size_t i) 
#else
	for (size_t i = 0; i < m_Spec.sizeX; i++)
#endif
	{
		for (size_t j = 0; j < m_Spec.sizeY; j++)
		{
			auto& node = m_Nodes[i][j];

			if (node.IsSolid) continue;

			node.UpdateMacroscopic();
		}
	}
#ifdef MULTITHREADED
	);
#endif
}

void Lattice::StreamingStep()
{
#ifdef MULTITHREADED
	std::for_each(std::execution::par, m_Indices.begin(), m_Indices.end(), [this](size_t i) 
	{
		const auto& sizeX = m_Spec.sizeX;
		const auto& sizeY = m_Spec.sizeY;
#else
	const auto& sizeX = m_Spec.sizeX;
	const auto& sizeY = m_Spec.sizeY;

	for (size_t i = 0; i < sizeX; i++)
	{
#endif
		//Neighbor node ids: left, right
		const size_t l = (i != 0) ? i - 1 : sizeX - 1;
		const size_t r = (i != sizeX - 1) ? i + 1 : 0;

		for (size_t j = 0; j < sizeY; j++)
		{
			//Skip solid interior nodes
			if (m_Nodes[i][j].IsSolidInterior) continue;

			//Neightbor node ids: up, down
			const size_t u = (j != sizeY - 1) ? j + 1 : 0;
			const size_t d = (j != 0) ? j - 1 : sizeY - 1;

			//Uses following convention from "Lattice Boltzman Modelling 
			// An Introduction for Geoscientists and Engineers":
			//   6 -- 2 -- 5
			//   |    |    |
			//   3 -- 0 -- 1
			//   |    |    |
			//   7 -- 4 -- 8

			m_Nodes[i][j].TmpWeights[0] = m_Nodes[i][j].Weights[0];
			m_Nodes[r][j].TmpWeights[1] = m_Nodes[i][j].Weights[1];
			m_Nodes[i][u].TmpWeights[2] = m_Nodes[i][j].Weights[2];
			m_Nodes[l][j].TmpWeights[3] = m_Nodes[i][j].Weights[3];
			m_Nodes[i][d].TmpWeights[4] = m_Nodes[i][j].Weights[4];
			m_Nodes[r][u].TmpWeights[5] = m_Nodes[i][j].Weights[5];
			m_Nodes[l][u].TmpWeights[6] = m_Nodes[i][j].Weights[6];
			m_Nodes[l][d].TmpWeights[7] = m_Nodes[i][j].Weights[7];
			m_Nodes[r][d].TmpWeights[8] = m_Nodes[i][j].Weights[8];
		}
	}

#ifdef MULTITHREADED
	);
#endif
}

void Lattice::CollisionAndBounceback()
{
#ifdef MULTITHREADED
	std::for_each(std::execution::par, m_Indices.begin(), m_Indices.end(), [this](size_t i) 
#else
	for (size_t i = 0; i < m_Spec.sizeX; i++)
#endif
	{
		for (size_t j = 0; j < m_Spec.sizeY; j++)
		{
			auto& node = m_Nodes[i][j];

			//Skip solid interior nodes
			if (node.IsSolidInterior) continue;

			//Perform bounceback on solid boundary nodes:
			if (node.IsSolid)
			{
				auto& f_in = node.TmpWeights;
				auto& f_out = node.Weights;

				f_out[0] = f_in[0];
				f_out[1] = f_in[3];
				f_out[2] = f_in[4];
				f_out[3] = f_in[1];
				f_out[4] = f_in[2];
				f_out[5] = f_in[7];
				f_out[6] = f_in[8];
				f_out[7] = f_in[5];
				f_out[8] = f_in[6];
			}

			//Perform relaxation towards local equilibrium otherwise
			else
			{
				for (size_t a = 0; a < 9; a++)
				{
					//Relaxation time (in appropriate units, as discussed in equilibrium distribution)
					//times the force - in this case downwards pointing gravity
					const double F = node.Density * m_Spec.MassUnit * m_Spec.Gravity;
					const Utils::Vec2 tauF{0.0, - F * m_TauForce};

					const auto f_tmp = node.TmpWeights[a];
					const auto f_eq = node.Equlibrium(a, tauF);

					node.Weights[a] = f_tmp - m_Omega * (f_tmp - f_eq);
				}
			}

		}
	}

#ifdef MULTITHREADED
	);
#endif
}

void Lattice::HandleBoundaries()
{
	for (int i = 0; i < 4; i++)
	{
		//Periodic conditions need no further handling
		if (m_Spec.BoundaryConditions[i] == BoundaryCondition::Periodic)
			continue;

		HandleBoundary(static_cast<Boundary>(i));
	}
}

void Lattice::HandleBoundary(Boundary boundary)
{
	//This scheme of considering boundary conditions currently 
	// fails when applied on two consecutive edges.
	//To-do: consider some corner conditions?

	const auto& sizeX = m_Spec.sizeX;
	const auto& sizeY = m_Spec.sizeY;

	//Determine iteration range
	bool horizontal_boundary = (boundary == Up || boundary == Down);

	size_t max_id = horizontal_boundary ? sizeX : sizeY;

	//Select appropriate indices
	typedef std::array<size_t, 3> triplet;

	constexpr std::array<triplet, 4> same_options{
		triplet{1, 5, 8}, triplet{2, 5, 6}, triplet{3, 7, 6}, triplet{4, 7, 8},
	};

	constexpr std::array<triplet, 4> middle_options{
		triplet{2, 0, 4}, triplet{1, 0, 3}, triplet{2, 0, 4}, triplet{1, 0, 3}
	};

	constexpr std::array<triplet, 4> opposite_options{
		triplet{3, 7, 6}, triplet{4, 7, 8}, triplet{1, 5, 8}, triplet{2, 5, 6}
	};

	const triplet same     = same_options[boundary];
	const triplet middle   = middle_options[boundary];
	const triplet opposite = opposite_options[boundary];
	
	double sgn = (boundary == Up || boundary == Right) ? 1.0 : -1.0;

	//Iterate over appropriate edge of the lattice
	for (size_t i = 0; i < max_id; i++)
	{
		Node* node = nullptr;

		switch (boundary)
		{
			case Up:    {node = &m_Nodes[i][sizeY - 1]; break;}
			case Down:  {node = &m_Nodes[i][0];		    break;}
			case Left:  {node = &m_Nodes[0][i];		    break;}
			case Right: {node = &m_Nodes[sizeX - 1][i]; break;}
		}
		
		
		if (node)
		{
			if (node->IsSolid) continue;

			auto& fi = (*node).TmpWeights;

			const double sum_middle = fi[middle[0]] + fi[middle[1]] + fi[middle[2]];
			const double sum_opposite = fi[opposite[0]] + fi[opposite[1]] + fi[opposite[2]];

			switch (m_Spec.BoundaryConditions[boundary])
			{
			case BoundaryCondition::VonNeumann:
			{
				//Calculate Zou and He velocity BCs
				double VonNeumanVel = m_Spec.VonNeumannVelocitiesNormal[boundary];

				double rho0 = (sum_middle + 2.0 * sum_opposite) / (1.0 - sgn * VonNeumanVel);

				double ru = rho0 * VonNeumanVel;

				fi[same[0]] = fi[opposite[0]] + sgn * (2.0 / 3.0) * ru;
				fi[same[1]] = fi[opposite[1]] + sgn * (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[0]] - fi[middle[2]]);
				fi[same[2]] = fi[opposite[2]] + sgn * (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[2]] - fi[middle[0]]);

				break;
			}
			case BoundaryCondition::Dirichlet:
			{
				//Calculate Zou and He pressure BCs
				double rho0 = m_Spec.DirichletDensities[boundary];

				double u0 = -1.0 + (sum_middle + 2.0 * sum_opposite) / rho0;

				double ru = rho0 * u0;

				fi[same[0]] = fi[opposite[0]] - (2.0 / 3.0) * ru;
				fi[same[1]] = fi[opposite[1]] - (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[0]] - fi[middle[2]]);
				fi[same[2]] = fi[opposite[2]] - (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[2]] - fi[middle[0]]);

				break;
			}
			}
		}

		else
		{
			throw std::invalid_argument("Node pointer is null.");
		}
		
	}
}

void Lattice::LoadScene(const Scene& scene)
{
	//Test against the scene, if a node is inside some shape mark it as solid
	for (size_t idx = 0; idx < m_Nodes.size(); idx++)
	{
		for (size_t idy = 0; idy < m_Nodes[0].size(); idy++)
		{
			Vec2 pos{
				static_cast<double>(idx),
				static_cast<double>(idy)
			};

			m_Nodes[idx][idy].IsSolid = scene.IsInside(pos);
		}
	}

	//If a node and all its neighbors are solid mark it as solid_interior
	const auto sizeX = m_Nodes.size();
	const auto sizeY = m_Nodes[0].size();

	for (size_t idx = 0; idx < sizeX; idx++)
	{
		const size_t l = (idx != 0) ? idx - 1 : sizeX - 1;
		const size_t r = (idx != sizeX - 1) ? idx + 1 : 0;

		for (size_t idy = 0; idy < sizeY; idy++)
		{
			const size_t u = (idy != sizeY - 1) ? idy + 1 : 0;
			const size_t d = (idy != 0) ? idy - 1 : sizeY - 1;

			bool is_solid_interior = m_Nodes[idx][idy].IsSolid;

			is_solid_interior = is_solid_interior && m_Nodes[l][idy].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[r][idy].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[idx][u].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[idx][d].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[l][u].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[r][d].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[r][u].IsSolid;
			is_solid_interior = is_solid_interior && m_Nodes[l][d].IsSolid;

			m_Nodes[idx][idy].IsSolidInterior = is_solid_interior;
		}
	}
}

void Lattice::InitFlow(void func(size_t, size_t, Node&))
{
	for (size_t idx = 0; idx < m_Nodes.size(); idx++)
	{
		for (size_t idy = 0; idy < m_Nodes[0].size(); idy++)
		{
			auto& node = m_Nodes[idx][idy];

			if (node.IsSolid) continue;

			func(idx, idy, node);
		}
	}
}

void Lattice::Serialize(std::filesystem::path filepath)
{
	std::ofstream output(filepath, std::ios::trunc);

	if (output)
	{
		for (const auto& row : m_Nodes)
		{
			for (const auto& node : row)
			{
				output << std::fixed << std::setprecision(8) << node.Density    << " ";
				output << std::fixed << std::setprecision(8) << m_BaseSpeed * node.Velocity.x << " ";
				output << std::fixed << std::setprecision(8) << m_BaseSpeed * node.Velocity.y << " ";
			}

			output << '\n';
		}
	}

	else
	{
		std::cerr << "Failed to open file output stream: " << filepath << '\n';
	}
}
