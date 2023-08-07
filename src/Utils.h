#pragma once

namespace Utils {
	//Temporary, will probably switch to some linalg library like eigen
	struct Vec2 {
		double x, y;
	};

	double length(Vec2 v);
	double dot(Vec2 v1, Vec2 v2);
};

Utils::Vec2 operator+(const Utils::Vec2 v1, const Utils::Vec2 v2);
Utils::Vec2 operator-(const Utils::Vec2 v1, const Utils::Vec2 v2);
Utils::Vec2 operator*(const double a, const Utils::Vec2 v);
Utils::Vec2 operator/(const Utils::Vec2 v, const double a);