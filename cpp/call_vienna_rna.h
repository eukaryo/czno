/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
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

namespace czno_cpp {

double EnergyOfStructure(const std::string& sequence, const std::string& structure);
std::vector<std::string> SampleStructure(const std::string& sequence, const int num);

void ResetCount();
int GetCount();

}


#endif//CALL_VIENNA_RNA_H_
