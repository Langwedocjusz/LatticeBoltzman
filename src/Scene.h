#pragma once

#include "Utils.h"

#include <vector>
#include <memory>

class Shape {
public:
	virtual bool IsInside(Utils::Vec2 p) const;
};

class AABB : public Shape {
public:
	AABB(Utils::Vec2 min, Utils::Vec2 max);

	bool IsInside(Utils::Vec2 p) const override;

private:
	Utils::Vec2 m_Min, m_Max;
};

class Circle : public Shape {
public:
	Circle(Utils::Vec2 pos, double rad);

	bool IsInside(Utils::Vec2 p) const override;

private:
	Utils::Vec2 m_Pos;
	double m_Radius;
};

class Scene {
public:
	Scene();

	template<class T, typename ... Args>
	void PushShape(Args ... args)
	{
		static_assert(std::is_base_of<Shape, T>::value,
			"Add template argument not derived from Shape"
		);

		m_Shapes.push_back(std::unique_ptr<Shape>(new T(args...)));
	}

	bool IsInside(Utils::Vec2 p) const;

private:
	std::vector<std::unique_ptr<Shape>> m_Shapes;
};