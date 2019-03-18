/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef MISC_H_
#define MISC_H_

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

namespace czno_cpp {

const int TURN = 3;
std::vector<std::string> Split(const std::string& str, const char delimiter);
bool IsValidBasePair(const std::string& sequence, const std::pair<int, int>& x);
std::set<std::pair<int, int>> DotNotationToBasePairSet(const std::string& sequence);
std::string BasePairSetToDotNotation(const std::set<std::pair<int, int>>& structure, const int len);
int HammingDistance(const std::string& structure1, const std::string& structure2);
bool IsExclusive(const std::pair<int, int>& base_pair1, const std::pair<int, int>& base_pair2);
void ValidateTriple(const std::string& sequence, const std::string& structure1, const std::string& structure2);
std::vector<std::string>TrimPathway(const std::vector<std::string>& pathway);
double BarrierEnergy(const std::string& sequence, const std::vector<std::string>& pathway);
double BarrierEnergy(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list);
std::pair<int, int>Difference(const std::string& structure1, const std::string& structure2);
std::vector<std::pair<int, int>> PathwayToActionList(const std::vector<std::string>& pathway);
bool IsValidPathway(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list);
bool IsValidPathway(const std::string& A, const std::string& B, const std::vector<std::string>& pathway);

}


#endif//MISC_H_