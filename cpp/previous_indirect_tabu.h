/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef PREVIOUS_DIRECT_TABU_H_
#define PREVIOUS_DIRECT_TABU_H_

#include <iostream>
#include <fstream>
#include <algorithm>
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

std::pair<std::vector<std::string>, double>
Dotu2010TabuIndirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int k,
	const double w_0,
	const int max_stable,
	const double energy_diff_bound,
	const int random_seed);

}

#endif//PREVIOUS_DIRECT_TABU_H_