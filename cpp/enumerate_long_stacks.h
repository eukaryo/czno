/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/
#ifndef ENUMERATE_LONG_STACKS_H_
#define ENUMERATE_LONG_STACKS_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <string>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <functional>
#include <complex>
#include <random>
#include <chrono>
#include <sstream>

namespace czno_cpp {

std::vector<std::array<int, 4>>
EnumerateLongStacksBruteForce(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& structure,
	const int minimum_length,
	const bool need_compatible);

std::vector<std::array<int, 4>>
EnumerateLongStacks(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& structure,
	const int minimum_length,
	const bool need_compatible);

}

#endif//ENUMERATE_LONG_STACKS_H_