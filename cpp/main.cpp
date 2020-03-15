/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include <iostream>
#include <string>
#include <algorithm>

#include "call_vienna_rna.h"
#include "previous_direct_methods.h"
#include "czno_proposed_methods.h"
#include "previous_indirect_ea.h"
#include "previous_indirect_mh1998.h"
#include "previous_indirect_tabu.h"
#include "previous_indirect_rna2dfold.h"
#include "misc.h"

#include "czno_experiments.h"

namespace czno_cpp {

int main_(int argc, char *argv[]) {


	if (argc == 2) {
		if (std::string(argv[1]) == std::string("hello")) {
			std::cout << czno_cpp::EnergyOfStructure(std::string("CCCCAAAAGGGG"), std::string("((((....))))")) << std::endl;
			const auto sampled = czno_cpp::SampleStructure(std::string("CCCCAAAAGGGGGGAAAACCCC"), 5);
			for (const auto x : sampled)std::cout << x << std::endl;
			return 0;
		}
		if (std::string(argv[1]) == std::string("time1")) {
			for (int i = 1; i <= 1600; ++i)TimeExperiment1(i - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("time2")) {
			for (int i = 1; i <= 1100; ++i)TimeExperiment2(i - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("accu1")) {
			for (int i = 1; i <= 1600; ++i)AccuracyExperiment1(i - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("indi1")) {
			for (int i = 1; i <= 1600; ++i)IndirectExperiment1(i - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("real1")) {
			for (int x = 1; x <= 18 * 109; ++x) {
				const int i = (x - 1) % 18;
				const int j = (x - 1) / 18;
				RealDataExperiment1(i, j);
			}
			return 0;
		}
		if (std::string(argv[1]) == std::string("getdata")) {
			GetBasicDataOfRealDataset();
			return 0;
		}
		if (std::string(argv[1]) == std::string("clus1")) {
			for (int i = 1; i <= 1600; ++i)ClusteringExperiment1(i - 1);
			return 0;
		}

		return 0;
	}
	if (argc == 3) {
		if (std::string(argv[1]) == std::string("time1")) {
			TimeExperiment1(std::stoi(std::string(argv[2])) - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("time2")) {
			TimeExperiment2(std::stoi(std::string(argv[2])) - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("accu1")) {
			AccuracyExperiment1(std::stoi(std::string(argv[2])) - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("indi1")) {
			IndirectExperiment1(std::stoi(std::string(argv[2])) - 1);
			return 0;
		}
		if (std::string(argv[1]) == std::string("real1")) {
			const int x = std::stoi(std::string(argv[2])) - 1;
			if (x < 0 || 18 * 109 <= x)return 0;
			const int i = x % 18;
			const int j = x / 18;
			RealDataExperiment1(i, j);
			return 0;
		}
		if (std::string(argv[1]) == std::string("clus1")) {
			ClusteringExperiment1(std::stoi(std::string(argv[2])) - 1);
			return 0;
		}
		return 0;
	}
	return 0;
}

}

int main(int argc, char *argv[]) {

	czno_cpp::main_(argc, argv);

	return 0;
}