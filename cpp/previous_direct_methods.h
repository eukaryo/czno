/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef PREVIOUS_DIRECT_METHODS_H_
#define PREVIOUS_DIRECT_METHODS_H_

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
MorganHiggs1998Direct(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2);

std::pair<std::vector<std::string>, double>
Voss2004Direct(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int k,
	const int seed);

std::pair<std::vector<std::string>, double>
Flamm2001FindpathDirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int maxkeep);

}

#endif//PREVIOUS_DIRECT_METHODS_H_