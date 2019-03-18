/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
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
	const int initial_population,//�ŏ��ɓK���ɍ��p�X�̐��B4
	const int gamma_generation,//�I�������Ɋւ��� 5
	const int max_generation,//�I�������Ɋւ��MAX 10
	const int script_capital_l,//�q���̐��Ɋւ��A100
	const int l1_elite_population,//�G���[�g�̐��A5
	const int l2_survive_population,// 5
	const int l3_final_max_population,// 100
	const int random_seed);

}

#endif//PREVIOUS_INDIRECT_EA_H_