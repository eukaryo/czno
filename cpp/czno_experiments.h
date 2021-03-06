﻿/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/


#ifndef CZNO_EXPERIMENTS_H_
#define CZNO_EXPERIMENTS_H_

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <cassert>
#include <regex>
#include <fstream>

namespace czno_cpp {

void TimeExperiment1(const int i);
void TimeExperiment2(const int i);
void AccuracyExperiment1(const int i);
void IndirectExperiment1(const int i);
void RealDataExperiment1(const int i, const int j);
void GetBasicDataOfRealDataset();

void ClusteringExperiment1(const int i);
}


#endif//CZNO_EXPERIMENTS_H_
