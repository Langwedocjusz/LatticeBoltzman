#include "Utils.h"

#include <cmath>


double Utils::length(Utils::Vec2 v)
{
	return std::sqrt(v.x*v.x + v.y*v.y);
}

double Utils::dot(Utils::Vec2 v1, Utils::Vec2 v2)
{
	return v1.x * v2.x + v1.y * v2.y;
}

Utils::Vec2 operator+(const Utils::Vec2 v1, const Utils::Vec2 v2)
{
	return Utils::Vec2{ v1.x + v2.x, v1.y + v2.y };
}

Utils::Vec2 operator-(const Utils::Vec2 v1, const Utils::Vec2 v2)
{
	return Utils::Vec2{ v1.x - v2.x, v1.y - v2.y };
}

Utils::Vec2 operator*(const double a, const Utils::Vec2 v)
{
	return Utils::Vec2{ a * v.x, a * v.y };
}

Utils::Vec2 operator/(const Utils::Vec2 v, const double a)
{
	return Utils::Vec2{v.x/a, v.y/a};
}