#pragma once

namespace Utils {
	//Temporary, will probably switch to some linalg library like eigen
	struct Vec2 {
		double x, y;
	};

	double dot(Vec2 v1, Vec2 v2);
};