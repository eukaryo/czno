/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"czno_experiments.h"

#include"misc.h"
#include"czno_proposed_methods.h"
#include"call_vienna_rna.h"
#include"previous_direct_methods.h"
#include"previous_indirect_ea.h"
#include"previous_indirect_mh1998.h"
#include"previous_indirect_tabu.h"

namespace czno_cpp {

std::vector<std::pair<std::string, std::pair<std::string, std::string>>>DatasetReader(const std::string& filename) {

	std::vector<std::pair<std::string, std::pair<std::string, std::string>>>answer;
	std::ifstream ifs(filename);
	if (ifs.fail())assert(0);
	while (1) {
		std::string line;
		if (!std::getline(ifs, line))break;
		const auto v = Split(line, ' ');
		assert(v.size() == 3);

		const static std::regex rna(R"([ACGU]+)");
		const static std::regex ss(R"([\(\.\)]+)");
		const auto func = [&](const std::string s, const std::regex r) {
			std::smatch m;
			std::regex_search(s, m, r);
			assert(m.size() == 1);
			return m[0].str();
		};
		answer.push_back(std::make_pair(func(v[0], rna), std::make_pair(func(v[1], ss), func(v[2], ss))));
		assert(answer.back().first.size() == answer.back().second.first.size());
		assert(answer.back().first.size() == answer.back().second.second.size());
	}

	return answer;
}

void TimeExperiment1(const int i) {
	//N=100, ハミング距離5～20

	const auto dataset = DatasetReader("random_dataset/dataset-iso-len.txt");
	if (i < 0 || dataset.size() <= i)return;

	const bool first_is_stable = EnergyOfStructure(dataset[i].first, dataset[i].second.first) < EnergyOfStructure(dataset[i].first, dataset[i].second.second);
	const auto stable = first_is_stable ? dataset[i].second.first : dataset[i].second.second;
	const auto instable = first_is_stable ? dataset[i].second.second : dataset[i].second.first;

	const int H = HammingDistance(dataset[i].second.first, dataset[i].second.second);
	ResetCount();
	std::cout << "LOG: start: proposed, new " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start1 = std::chrono::system_clock::now();
	const auto x1 = MinimumBarrierDirectPathDijkstra(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto end1 = std::chrono::system_clock::now();
	const int T1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
	const int C1 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T1 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=stable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start2 = std::chrono::system_clock::now();
	const auto x2 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, stable, instable);
	const auto end2 = std::chrono::system_clock::now();
	const int T2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
	const int C2 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T2 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=instable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start3 = std::chrono::system_clock::now();
	const auto x3 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, instable, stable);
	const auto end3 = std::chrono::system_clock::now();
	const int T3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();
	const int C3 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T3 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: findpath " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start4 = std::chrono::system_clock::now();
	const auto x4 = Flamm2001FindpathDirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 1024 * 1024 * 64);
	const auto end4 = std::chrono::system_clock::now();
	const int T4 = std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count();
	const int C4 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T4 << " ms" << std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x1.first));
	assert(IsValidPathway(stable, instable, x2.first));
	assert(IsValidPathway(instable, stable, x3.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x4.first));
	assert(x1.second == x2.second);
	assert(x1.second == x3.second);
	assert(x1.second == x4.second);
	std::cout << "RESULT: " << i << " " << dataset[i].first.size() << " " << H << " " <<
		T1 << " " << C1 << " " <<
		T2 << " " << C2 << " " <<
		T3 << " " << C3 << " " <<
		T4 << " " << C4 << " " << x1.second << std::endl;
	return;
}
void TimeExperiment1x(const int i) {
	//N=100, ハミング距離5～20

	const auto dataset = DatasetReader("random_dataset/dataset-iso-len.txt");
	if (i < 0 || dataset.size() <= i)return;

	const bool first_is_stable = EnergyOfStructure(dataset[i].first, dataset[i].second.first) < EnergyOfStructure(dataset[i].first, dataset[i].second.second);
	const auto stable = first_is_stable ? dataset[i].second.first : dataset[i].second.second;
	const auto instable = first_is_stable ? dataset[i].second.second : dataset[i].second.first;

	const int H = HammingDistance(dataset[i].second.first, dataset[i].second.second);
	ResetCount();
	std::cout << "LOG: start: proposed, new " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start1 = std::chrono::system_clock::now();
	const auto x1 = MinimumBarrierDirectPathDijkstra(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto end1 = std::chrono::system_clock::now();
	const int T1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
	const int C1 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T1 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=stable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start2 = std::chrono::system_clock::now();
	const auto x2 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, stable, instable);
	const auto end2 = std::chrono::system_clock::now();
	const int T2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
	const int C2 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T2 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=instable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start3 = std::chrono::system_clock::now();
	const auto x3 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, instable, stable);
	const auto end3 = std::chrono::system_clock::now();
	const int T3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();
	const int C3 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T3 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: findpath " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start4 = std::chrono::system_clock::now();
	const auto x4 = Flamm2001FindpathDirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 1024 * 1024 * 64);
	const auto end4 = std::chrono::system_clock::now();
	const int T4 = std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count();
	const int C4 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T4 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: new,coroutine " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start5 = std::chrono::system_clock::now();
	const auto x5 = MinimumBarrierDirectPathDijkstra1Turn(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto end5 = std::chrono::system_clock::now();
	const int T5 = std::chrono::duration_cast<std::chrono::milliseconds>(end5 - start5).count();
	const int C5 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T5 << " ms" << std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x1.first));
	assert(IsValidPathway(stable, instable, x2.first));
	assert(IsValidPathway(instable, stable, x3.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x4.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x5.first));
	assert(x1.second == x2.second);
	assert(x1.second == x3.second);
	assert(x1.second == x4.second);
	assert(x1.second == x5.second);
	std::cout << "RESULT: " << i << " " << dataset[i].first.size() << " " << H << " " <<
		T1 << " " << C1 << " " <<
		T2 << " " << C2 << " " <<
		T3 << " " << C3 << " " <<
		T4 << " " << C4 << " " <<
		T5 << " " << C5 << " " << x1.second << std::endl;
	return;
}
void TimeExperiment2(const int i) {
	//ハミング距離18,N＝50～150

	const auto dataset = DatasetReader("random_dataset/dataset-iso-dist.txt");
	if (i < 0 || dataset.size() <= i)return;

	const bool first_is_stable = EnergyOfStructure(dataset[i].first, dataset[i].second.first) < EnergyOfStructure(dataset[i].first, dataset[i].second.second);
	const auto stable = first_is_stable ? dataset[i].second.first : dataset[i].second.second;
	const auto instable = first_is_stable ? dataset[i].second.second : dataset[i].second.first;

	const int H = HammingDistance(dataset[i].second.first, dataset[i].second.second);
	ResetCount();
	std::cout << "LOG: start: proposed, new " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start1 = std::chrono::system_clock::now();
	const auto x1 = MinimumBarrierDirectPathDijkstra(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto end1 = std::chrono::system_clock::now();
	const int T1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
	const int C1 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T1 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=stable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start2 = std::chrono::system_clock::now();
	const auto x2 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, stable, instable);
	const auto end2 = std::chrono::system_clock::now();
	const int T2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
	const int C2 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T2 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: proposed, old, start=instable " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start3 = std::chrono::system_clock::now();
	const auto x3 = MinimumBarrierDirectPathDijkstraOld(dataset[i].first, instable, stable);
	const auto end3 = std::chrono::system_clock::now();
	const int T3 = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();
	const int C3 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T3 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: findpath " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start4 = std::chrono::system_clock::now();
	const auto x4 = Flamm2001FindpathDirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 1024 * 1024 * 64);
	const auto end4 = std::chrono::system_clock::now();
	const int T4 = std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count();
	const int C4 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T4 << " ms" << std::endl;
	ResetCount();
	std::cout << "LOG: start: new,coroutine " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	const auto start5 = std::chrono::system_clock::now();
	const auto x5 = MinimumBarrierDirectPathDijkstra1Turn(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto end5 = std::chrono::system_clock::now();
	const int T5 = std::chrono::duration_cast<std::chrono::milliseconds>(end5 - start5).count();
	const int C5 = GetCount();
	std::cout << "LOG: end: " << i << "/" << dataset.size() << ",  elapsed time = " << T5 << " ms" << std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x1.first));
	assert(IsValidPathway(stable, instable, x2.first));
	assert(IsValidPathway(instable, stable, x3.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x4.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x5.first));
	assert(x1.second == x2.second);
	assert(x1.second == x3.second);
	assert(x1.second == x4.second);
	assert(x1.second == x5.second);
	std::cout << "RESULT: " << i << " " << dataset[i].first.size() << " " << H << " " <<
		T1 << " " << C1 << " " <<
		T2 << " " << C2 << " " <<
		T3 << " " << C3 << " " <<
		T4 << " " << C4 << " " <<
		T5 << " " << C5 << " " << x1.second << std::endl;
	return;
}

void AccuracyExperiment1(const int i) {
	//direct pathwayを求める方法でバリアの値を比べる。

	const auto dataset = DatasetReader("random_dataset/dataset-iso-len.txt");
	if (i < 0 || dataset.size() <= i)return;

	const auto func = [](const std::string& sequence, const std::string& structure1, const std::string& structure2) {
		auto x1 = Voss2004Direct(sequence, structure1, structure2, 10, 123456);
		for (int j = 2; j <= 1000; ++j) {
			const auto x2 = Voss2004Direct(sequence, structure1, structure2, 10, 123456 + j);
			if (x1.second > x2.second)x1 = x2;
		}
		return x1;
	};

	const int H = HammingDistance(dataset[i].second.first, dataset[i].second.second);
	ResetCount();
	std::cout << "LOG: start: " << i << "/" << dataset.size() << ", Hamming distance = " << H << std::endl;
	const auto x1 = MinimumBarrierDirectPathDijkstra(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto x2 = MorganHiggs1998Direct(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	const auto x3 = Voss2004Direct(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 10, 12345);
	const auto x4 = Flamm2001FindpathDirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 10);
	const auto x5 = func(dataset[i].first, dataset[i].second.first, dataset[i].second.second);
	std::cout << "LOG: end: " << i << "/" << dataset.size() << std::endl;
	std::cout << "RESULT: " << i << " " << x1.second << " " << x2.second << " " << x3.second << " " << x4.second << " " << x5.second << " " << std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x1.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x2.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x3.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x4.first));
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, x5.first));
	return;
}

void IndirectExperiment1(const int i) {
	//indirect pathwayを改善させられるかやってみる

	const auto dataset = DatasetReader("random_dataset/dataset-iso-len.txt");
	if (i < 0 || dataset.size() <= i)return;

	const int H = HammingDistance(dataset[i].second.first, dataset[i].second.second);
	ResetCount();
	std::cout << "LOG: start: " << i << "/" << dataset.size() << ", length = " << dataset[i].first.size() << ", Hamming distance = " << H << std::endl;
	std::cout << "LOG: start: 5" << std::endl;
	const auto Y5 = LiZhang2012RNAEAPathIndirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 4, 5, 10, 100, 5, 5, 100, 12345);
	std::cout << "LOG: start: 6" << std::endl;
	const auto Y6 = MorganHiggs1998Indirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 10);
	std::cout << "LOG: start: 7" << std::endl;
	const auto Y7 = Dotu2010TabuIndirect(dataset[i].first, dataset[i].second.first, dataset[i].second.second, 5, 4, 7, 10000000000000.0, 12345);
	std::cout << "LOG: end: " << i << "/" << dataset.size() << std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, Y5.first) || Y5.first.size() == 0);
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, Y6.first) || Y6.first.size() == 0);
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, Y7.first) || Y7.first.size() == 0);

	const auto Improve = [&](const std::pair<std::vector<std::string>, double> first, const int num) {
		std::cout << "LOG: improve start: (data-id, method-id) = (" << i << ", " << num << ")" << std::endl;
		auto best = first;
		while (true) {
			const auto p = ImprovePathway(dataset[i].first, best.first, 2, 15);
			const auto s = BarrierEnergy(dataset[i].first, p);
			if (best.second > s) {
				std::cout << "LOG: success: (data-id, method-id) = (" << i << ", " << num << "), " << best.second << " -> " << s << std::endl;
				best = std::make_pair(p, s);
			}
			else break;
		}
		return best;
	};

	const auto YY5 = Y5.first.size() == 0 ? Y5 : Improve(Y5, 5);
	const auto YY6 = Y6.first.size() == 0 ? Y6 : Improve(Y6, 6);
	const auto YY7 = Y7.first.size() == 0 ? Y7 : Improve(Y7, 7);

	std::cout << "LOG: " << i << " " << H << " 5 " << Y5.second << " " << YY5.second << std::endl;
	std::cout << "LOG: " << i << " " << H << " 6 " << Y6.second << " " << YY6.second << std::endl;
	std::cout << "LOG: " << i << " " << H << " 7 " << Y7.second << " " << YY7.second << std::endl;
	std::cout << "RESULT: " << i << " " << H << " " <<
		Y5.second << " " << YY5.second << " " <<
		Y6.second << " " << YY6.second << " " <<
		Y7.second << " " << YY7.second <<
		std::endl;
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, YY5.first) || YY5.first.size() == 0);
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, YY6.first) || YY6.first.size() == 0);
	assert(IsValidPathway(dataset[i].second.first, dataset[i].second.second, YY7.first) || YY7.first.size() == 0);

	return;
}

}