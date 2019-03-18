/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "czno_proposed_methods.h"

#include "call_vienna_rna.h"
#include "misc.h"

namespace czno_cpp {

std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2) {

	const auto A = DotNotationToBasePairSet(structure1);
	const auto B = DotNotationToBasePairSet(structure2);

	std::vector<std::pair<int, int>>basepair_difference;
	std::set_symmetric_difference(
		A.begin(), A.end(),
		B.begin(), B.end(),
		inserter(basepair_difference, basepair_difference.end()));

	const int hamming_distance = basepair_difference.size();
	if (hamming_distance > 25)return std::make_pair(std::vector<std::string>{}, std::numeric_limits<double>::infinity());

	std::vector<std::pair<int, int>>sharing_basepairs;
	std::set_intersection(
		A.begin(), A.end(),
		B.begin(), B.end(),
		inserter(sharing_basepairs, sharing_basepairs.end()));

	std::string sharing_structure("");
	for (int i = 0; i < sequence.size(); ++i)sharing_structure += ".";
	for (const auto x : sharing_basepairs) {
		assert(x.first < x.second);
		sharing_structure[x.first] = '(';
		sharing_structure[x.second] = ')';
	}

	std::vector<int>need_to_open(hamming_distance, 0);
	for (int i = 0; i < hamming_distance; ++i) {
		need_to_open[i] = (A.find(basepair_difference[i]) != A.end()) ? 1 : 0;
	}

	const auto DotNotation = [&](const int structure_index) {
		std::string answer = sharing_structure;
		for (int i = 0; i < hamming_distance; ++i) {
			if (((structure_index & (1 << i)) ? 1 : 0) ^ need_to_open[i]) {
				assert(basepair_difference[i].first < basepair_difference[i].second);
				assert(answer[basepair_difference[i].first] == '.');
				assert(answer[basepair_difference[i].second] == '.');
				answer[basepair_difference[i].first] = '(';
				answer[basepair_difference[i].second] = ')';
			}
		}
		return answer;
	};

	std::vector<int> prerequisite_flag(hamming_distance, 0);
	for (int i = 0; i < hamming_distance; ++i) {
		if (need_to_open[i] == 0)continue;
		for (int j = 0; j < hamming_distance; ++j) {
			if (need_to_open[j] != 0)continue;

			//この時点で、
			//basepair_difference[i]はAにあってBにない、すなわちパス上で外されるべき塩基対であり、
			//basepair_difference[j]はBにあってAにない、すなわちパス上で組まれるべき塩基対である。
			//ここでもし[i]を外してからでないと[j]を組めないという依存関係があるなら、
			//prerequisite_flag[j] |= 1 << i とする。

			if (IsExclusive(basepair_difference[i], basepair_difference[j])) {
				prerequisite_flag[j] |= 1 << i;
			}
		}
	}

	//スタート構造から構造xまでの部分解。配るDPをするので、その時点で判明している最良の解を入れる。
	std::vector<double>result(1 << hamming_distance, std::numeric_limits<double>::infinity());

	//構造xのエネルギーの値
	std::vector<double>energy(1 << hamming_distance, std::numeric_limits<double>::infinity());

	//Resultが確定したかどうか
	std::vector<int>searched(1 << hamming_distance, 0);

	//DNodeは(スタート構造からその構造までのバリア値、(-時刻、その構造の添字))とする。
	//"-時刻"を入れる理由は、バリア値が同じ場合に後入れ先出しにしたいから。
	//構造の添字とは、すべての可能な中間構造をBitDP的に添字付けしたものとする。
	typedef std::pair<double, std::pair<int, int>> DNode;

	std::priority_queue<DNode, std::vector<DNode>, std::greater<DNode>>dijkstra;
	int count = 0;
	energy[0] = EnergyOfStructure(sequence, structure1);
	dijkstra.push(std::make_pair(energy[0], std::make_pair(count, 0)));
	while (!dijkstra.empty()) {
		count--;
		const auto x = dijkstra.top();
		dijkstra.pop();
		const int structure_index = x.second.second;
		if (searched[structure_index]++)continue;
		result[structure_index] = x.first;
		if (structure_index == (1 << hamming_distance) - 1)break;
		for (int i = 0; i < hamming_distance; ++i) {
			const int b = 1 << i;
			if (structure_index & b)continue;
			if ((structure_index & prerequisite_flag[i]) != prerequisite_flag[i])continue;
			const int new_structure_index = structure_index | b;
			if (searched[new_structure_index])continue;
			if (energy[new_structure_index] == std::numeric_limits<double>::infinity()) {
				energy[new_structure_index] = EnergyOfStructure(sequence, DotNotation(new_structure_index));
			}
			const double new_result = std::max<double>(result[structure_index], energy[new_structure_index]);
			if (new_result >= result[new_structure_index])continue;
			result[new_structure_index] = new_result;
			dijkstra.push(std::make_pair(new_result, std::make_pair(count, new_structure_index)));
		}
	}

	//トレースバック
	std::vector<std::string>reversed_answer{ structure2 };
	for (int structure_index = (1 << hamming_distance) - 1; structure_index;) {
		bool flag = false;
		for (int i = 0; i < hamming_distance; ++i) {
			const int b = 1 << i;
			if ((structure_index & b) == 0)continue;
			const int back_index = structure_index ^ b;
			if ((back_index & prerequisite_flag[i]) != prerequisite_flag[i])continue;
			if (searched[back_index] == 0)continue;
			if (result[structure_index] != std::max<double>(result[back_index], energy[structure_index]))continue;
			reversed_answer.push_back(DotNotation(back_index));
			structure_index = back_index;
			flag = true;
			break;
		}
		assert(flag);
	}
	assert(reversed_answer.back() == structure1);
	std::reverse(reversed_answer.begin(), reversed_answer.end());
	return std::make_pair(reversed_answer, result[(1 << hamming_distance) - 1]);
}

std::vector<std::string>
ImprovePathway(
	const std::string& sequence,
	const std::vector<std::string>& pathway,
	const int min_czno_len,
	const int max_czno_len) {

	std::vector<double>energy;
	double barrier = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < pathway.size(); ++i) {
		energy.push_back(EnergyOfStructure(sequence, pathway[i]));
		barrier = std::max<double>(barrier, energy.back());
	}

	std::vector<int>barrier_pos;
	for (int i = 0; i < pathway.size(); ++i)if (energy[i] == barrier)barrier_pos.push_back(i);

	for (int i = 0; i < pathway.size(); ++i) {
		for (int j = i + min_czno_len - 1; j < pathway.size() && j <= i + max_czno_len - 1; ++j) {
			for (const int x : barrier_pos)if (i < x && x < j) {
				const auto czno_result = MinimumBarrierDirectPathDijkstra(sequence, pathway[i], pathway[j]);
				std::vector<std::string>answer;
				for (int k = 0; k < i; ++k)answer.push_back(pathway[k]);
				for (const auto y : czno_result.first)answer.push_back(y);
				for (int k = j + 1; k < pathway.size(); ++k)answer.push_back(pathway[k]);
				return TrimPathway(answer);
			}
		}
	}
	return pathway;
}

}
