#pragma once

#ifndef MAIN_H_
#define MAIN_H_

#include <cmath>
#include <chrono>

#include </opt/homebrew/Cellar/boost/1.78.0_1/include/boost/math/special_functions/bessel.hpp>

class Timer
{
private:
	std::chrono::time_point <std::chrono::high_resolution_clock> start;

public:
	Timer() : start(std::chrono::high_resolution_clock::now()) {}

	void reset() { start = std::chrono::high_resolution_clock::now(); }

	double elapsed()
	{
		return (std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>
			(std::chrono::high_resolution_clock::now() - start)).count();
	}
};

#endif