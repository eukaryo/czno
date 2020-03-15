/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#ifndef CALL_VIENNA_RNA_H_
#define CALL_VIENNA_RNA_H_

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <cassert>
#include <regex>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

namespace czno_cpp {

double EnergyOfStructure(const std::string& sequence, const std::string& structure);
std::vector<std::string> SampleStructure(const std::string& sequence, const int num);
std::vector<std::pair<std::pair<int, int>, std::pair<double, std::string>>> CallRNA2Dfold(const std::string& sequence, const std::string& structure1, const std::string& structure2);

void ResetCount();
int GetCount();

}


#endif//CALL_VIENNA_RNA_H_
