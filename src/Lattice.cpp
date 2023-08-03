#include "Lattice.h"

#include <iostream>
#include <fstream>

using namespace Utils;

void Node::UpdateMacroscopic()
{
	Density = 0.0;
	Velocity = Vec2{ 0.0, 0.0 };

	//Macroscopic density is the sum of all weights
	for (const auto& weight : Weights)
	{
		Density += weight;
	}

	//Do not attempt to calculate velocity if density is zero, as that is ill defined
	if (!(Density > 0.0)) return;

	//Macroscopic velocity is the weighted average of base velocities
	for (size_t i = 0; i < Weights.size(); i++)
	{
		Velocity.x += Weights[i] * s_BaseVelocities[i].x;
		Velocity.y += Weights[i] * s_BaseVelocities[i].y;
	}

	Velocity.x /= Density;
	Velocity.y /= Density;
}

double Node::Equlibrium(double base_speed, uint32_t idx)
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

	const double eDotU = dot(e, Velocity);
	const double u2 = dot(Velocity, Velocity);

	const double c2 = base_speed * base_speed;
	const double c4 = c2 * c2;

	return weights[idx] * Density *
		(1.0 + 3.0 * eDotU / c2 + 4.5 * (eDotU*eDotU) / c4 - 1.5 * u2 /c2);
}

Lattice::Lattice(Specification spec)
	: m_Spec(spec), m_BaseSpeed(spec.LengthUnit/spec.TimeStep)
{
	m_Nodes.resize(m_Spec.sizeX, std::vector<Node>(m_Spec.sizeY));

	//Convert tau from lu^2 / dt units:
	m_Tau *= m_Spec.LengthUnit * m_Spec.LengthUnit * m_Spec.LengthUnit / m_Spec.TimeStep;
}

Lattice::~Lattice()
{

}

void Lattice::LoadScene(const Scene& scene)
{
	//Test against the scene, if a node is inside some shape mark it as solid
	for (size_t idx = 0; idx < m_Nodes.size(); idx++)
	{
		for (size_t idy = 0; idy < m_Nodes[0].size(); idy++)
		{
			const auto& lu = m_Spec.LengthUnit;

			Vec2 pos{ 
				lu * static_cast<double>(idx), 
				lu * static_cast<double>(idy) 
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

void Lattice::Update()
{
	const auto& sizeX = m_Spec.sizeX;
	const auto& sizeY = m_Spec.sizeY;

	//Streaming step
	for (size_t i = 0; i < sizeX; i++)
	{
		//Neighbor node ids: left, right
		const size_t l = (i != 0)         ? i - 1 : sizeX - 1;
		const size_t r = (i != sizeX - 1) ? i + 1 : 0;          

		for (size_t j = 0; j < sizeY; j++)
		{
			//Skip solid interior nodes
			if (m_Nodes[i][j].IsSolidInterior) continue;

			//Update Macroscopic values
			m_Nodes[i][j].UpdateMacroscopic();

			//Neightbor node ids: up, down
			const size_t u = (j != sizeY - 1) ? j + 1 : 0;
			const size_t d = (j != 0)         ? j - 1 : sizeY - 1;

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

	//Consider boundary conditions
	for (int i = 0; i < 4; i++)
	{
		const auto& boundary_condition = m_Spec.BoundaryConditions[i];

		switch (boundary_condition)
		{
			case BoundaryCondition::VonNeumann:
			{
				CalculateVonNeumann(static_cast<Boundary>(i));
				break;
			}

			case BoundaryCondition::Dirichlet:
			{
				//TO-DO: write this
				break;
			}

			default:
			{
				//Periodic conditions need no further handling
				break;
			}
		}
	}

	//Collision step
	for (size_t i = 0; i < m_Nodes.size(); i++)
	{
		for (size_t j = 0; j < m_Nodes[0].size(); j++)
		{
			auto& node = m_Nodes[i][j];

			//Skip solid interior nodes
			if (node.IsSolidInterior) continue;

			//Perform bounceback on solid boundary nodes:
			if (node.IsSolid)
			{
				double tmp;
				auto& fij = node.TmpWeights;
				
				tmp = fij[1]; fij[1] = fij[3]; fij[3] = tmp;
				tmp = fij[2]; fij[2] = fij[4]; fij[4] = tmp;
				tmp = fij[5]; fij[5] = fij[7]; fij[7] = tmp;
				tmp = fij[6]; fij[6] = fij[8]; fij[8] = tmp;

				node.Weights = node.TmpWeights;
			}

			//Perform relaxation towards local equilibrium otherwise
			else
			{
				for (size_t a = 0; a < 9; a++)
				{
					const auto f_tmp = node.TmpWeights[a];
					const auto f_eq = node.Equlibrium(m_BaseSpeed, a);

					node.Weights[a] = f_tmp - (f_tmp - f_eq) / m_Tau;
				}
			}

		}
	}
}

void Lattice::CalculateVonNeumann(Boundary boundary)
{
	//This scheme of considering VonNeumann boundary conditions currently 
	// fails when applied on two consecutive edges.
	//To-do: consider some corner conditions?

	const auto& sizeX = m_Spec.sizeX;
	const auto& sizeY = m_Spec.sizeY;

	double VonNeumanVel = m_Spec.VonNeumannVelocitiesNormal[boundary];

	//Determine iteration range
	bool horizontal = (boundary == Up || boundary == Down);

	size_t max_id = horizontal ? sizeX : sizeY;

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
			case Down:  {node = &m_Nodes[i][0];		 break;}
			case Left:  {node = &m_Nodes[0][i];		 break;}
			case Right: {node = &m_Nodes[sizeX - 1][i]; break;}
		}
		
		//Calculate Zou and He velocity BCs
		auto& fi = (*node).TmpWeights;
		
		const double sum_middle = fi[middle[0]] + fi[middle[1]] + fi[middle[2]];
		const double sum_opposite = fi[opposite[0]] + fi[opposite[1]] + fi[opposite[2]];

		double rho0 = (sum_middle + 2.0 * sum_opposite) / (1.0 - sgn * VonNeumanVel);
		
		double ru = rho0 * VonNeumanVel;

		fi[same[0]] = fi[opposite[0]] + sgn * (2.0 / 3.0) * ru;
		fi[same[1]] = fi[opposite[1]] + sgn * (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[0]] - fi[middle[2]]);
		fi[same[2]] = fi[opposite[2]] + sgn * (1.0 / 6.0) * ru - sgn * 0.5 * (fi[middle[2]] - fi[middle[0]]);

		fi[opposite[0]] = fi[opposite[0]];
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
				output << node.Density << " " << node.Velocity.x << " " << node.Velocity.y << " ";
			}

			output << '\n';
		}
	}

	else
	{
		std::cerr << "Failed to open file output stream: " << filepath << '\n';
	}
}