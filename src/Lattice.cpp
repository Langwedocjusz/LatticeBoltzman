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

Lattice::Lattice(LatticeSpecification spec)
	: m_Spec(spec), m_BaseSpeed(spec.LengthUnit/spec.TimeStep)
{
	m_Nodes.resize(m_Spec.sizeX, std::vector<Node>(m_Spec.sizeY));
}

Lattice::~Lattice()
{

}

void Lattice::Initialize(void func(size_t, size_t, Node&))
{
	for (size_t idx = 0; idx < m_Nodes.size(); idx++)
	{
		for (size_t idy = 0; idy < m_Nodes.size(); idy++)
		{
			func(idx, idy, m_Nodes[idx][idy]);
		}
	}
}

void Lattice::Update()
{
	const auto& sizeX = m_Spec.sizeX;
	const auto& sizeY = m_Spec.sizeY;

	//Streaming step

	for (size_t i = 0; i < m_Nodes.size(); i++)
	{
		//Neighbor node ids: left, right
		const size_t l = (i != 0)         ? i - 1 : sizeX - 1;
		const size_t r = (i != sizeX - 1) ? i + 1 : 0;          

		for (size_t j = 0; j < m_Nodes[0].size(); j++)
		{
			//Skip solid nodes
			if (m_Nodes[i][j].IsSolid) continue;

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

	//Collision step
	for (size_t i = 0; i < m_Nodes.size(); i++)
	{
		for (size_t j = 0; j < m_Nodes[0].size(); j++)
		{
			auto& node = m_Nodes[i][j];

			//Skip solid nodes
			if (node.IsSolid) continue;

			for (size_t a = 0; a < 9; a++)
			{
				const auto f_tmp = node.TmpWeights[a];
				const auto f_eq = node.Equlibrium(m_BaseSpeed, a);

				node.Weights[a] = f_tmp - (f_tmp - f_eq) / m_Tau;
			}

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