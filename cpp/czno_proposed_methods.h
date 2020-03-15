/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef CZNO_PROPOSED_METHODS_H_
#define CZNO_PROPOSED_METHODS_H_

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
#include <regex>
#include <queue>

namespace czno_cpp {

std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra1Turn(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2);

std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2);

std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstraOld(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2);

std::vector<std::string>
ImprovePathway(
	const std::string& sequence,
	const std::vector<std::string>& pathway,
	const int min_czno_len,
	const int max_czno_len);
}


#endif//CZNO_PROPOSED_METHODS_H_
