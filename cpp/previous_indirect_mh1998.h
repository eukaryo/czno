/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef PREVIOUS_INDIRECT_MH1998_H_
#define PREVIOUS_INDIRECT_MH1998_H_

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
#include <queue>

namespace czno_cpp {
std::vector<std::string> DistinctSamples(const std::string & sequence, const std::string & structure1, const std::string & structure2, const int num);
std::pair<std::vector<std::string>, double>
MorganHiggs1998Indirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int num);
}

#endif//PREVIOUS_INDIRECT_MH1998_H_