/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef PREVIOUS_INDIRECT_EA_H_
#define PREVIOUS_INDIRECT_EA_H_

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
LiZhang2012RNAEAPathIndirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int initial_population,//最初に適当に作るパスの数。4
	const int gamma_generation,//終了条件に関わるγ 5
	const int max_generation,//終了条件に関わるMAX 10
	const int script_capital_l,//子供の数に関わる、100
	const int l1_elite_population,//エリートの数、5
	const int l2_survive_population,// 5
	const int l3_final_max_population,// 100
	const int random_seed);

}

#endif//PREVIOUS_INDIRECT_EA_H_