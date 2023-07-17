#include "Lattice.h"

using namespace Utils;

void Node::UpdateMacroscopic()
{
	density = 0.0;
	velocity = Vec2{ 0.0, 0.0 };

	//Macroscopic density is the sum of all weights
	for (const auto& weight : weights)
	{
		density += weight;
	}

	//Macroscopic velocity is the weighted average of base velocities
	for (size_t i = 0; i < weights.size(); i++)
	{
		velocity.x += weights[i] * s_BaseVelocities[i].x;
		velocity.y += weights[i] * s_BaseVelocities[i].y;
	}

	velocity.x /= density;
	velocity.y /= density;
}

double Node::Equlibrium(double base_speed, uint32_t idx)
{
	constexpr std::array<double, 9> weights{
		4.0 / 9.0,
		1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
		1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
	};

	Vec2 e = s_BaseVelocities[idx];

	const double eDotU = dot(e, velocity);
	const double u2 = dot(velocity, velocity);

	const double c2 = base_speed * base_speed;
	const double c4 = c2 * c2;

	return weights[idx] * density *
		(1.0 + 3.0 * eDotU / c2 + 4.5 * (eDotU*eDotU) / c4 - 1.5 * u2 /c2);
}

Lattice::Lattice(LatticeSpecification spec)
	: m_Spec(spec), m_BaseSpeed(spec.lengthUnit/spec.timeStep)
{
	m_Nodes.resize(m_Spec.sizeX, std::vector<Node>(m_Spec.sizeY));
}

Lattice::~Lattice()
{

}

void Lattice::Update()
{

}