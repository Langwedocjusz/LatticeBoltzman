#include "Scene.h"

bool Shape::IsInside(Utils::Vec2 p) const
{
	return false;
}

AABB::AABB(Utils::Vec2 min, Utils::Vec2 max)
	: m_Min(min), m_Max(max)
{}

bool AABB::IsInside(Utils::Vec2 p) const
{
	const bool in_x = (m_Min.x <= p.x) && (p.x <= m_Max.x);
	const bool in_y = (m_Min.y <= p.y) && (p.y <= m_Max.y);

	return in_x && in_y;
}

Circle::Circle(Utils::Vec2 pos, double rad)
	: m_Pos(pos), m_Radius(rad)
{}

bool Circle::IsInside(Utils::Vec2 p) const
{
	double dist = Utils::length(p - m_Pos);

	return (dist < m_Radius);
}

Scene::Scene()
{

}

bool Scene::IsInside(Utils::Vec2 p) const
{
	bool res = false;

	for (auto& shape : m_Shapes)
	{
		res |= shape->IsInside(p);
	}

	return res;
}