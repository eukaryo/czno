/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include "previous_indirect_rna2dfold.h"

#include "call_vienna_rna.h"
#include "misc.h"
#include "enumerate_long_stacks.h"
#include "previous_direct_methods.h"

namespace czno_cpp {


std::pair<std::vector<std::string>, double>
RNA2DFoldIndirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int maxdist) {

	//start構造はエネルギーが高い方
	if (EnergyOfStructure(sequence, structure1) < EnergyOfStructure(sequence, structure2)) {
		auto answer = RNA2DFoldIndirect(sequence, structure2, structure1, maxdist);
		std::reverse(answer.first.begin(), answer.first.end());
		return answer;
	}

	const auto mesh = CallRNA2Dfold(sequence, structure1, structure2);
	if(mesh.size() == 0)return make_pair(std::vector<std::string>{}, std::numeric_limits<double>::infinity());

	int start_index = -1, goal_index = -1;
	for (int i = 0; i < mesh.size(); ++i) {
		if (mesh[i].first.first == 0) {
			assert(start_index == -1);
			start_index = i;
		}
		if (mesh[i].first.second == 0) {
			assert(goal_index == -1);
			goal_index = i;
		}
	}

	std::vector<double>energy(mesh.size());
	for (int i = 0; i < mesh.size(); ++i)energy[i] = EnergyOfStructure(sequence, mesh[i].second.second);

	std::map<std::pair<int, int>, double>direct_energy_cache;

	const auto GetDirectBarrier = [&](const int a, const int b) {
		if (direct_energy_cache.find(std::make_pair(a, b)) != direct_energy_cache.end())return direct_energy_cache[std::make_pair(a, b)];
		direct_energy_cache[std::make_pair(a, b)] = Flamm2001FindpathDirect(sequence, mesh[a].second.second, mesh[b].second.second, 10).second;
		return direct_energy_cache[std::make_pair(a, b)];
	};

	double upper_bound = GetDirectBarrier(start_index, goal_index);

	std::vector<double>result_barrier(mesh.size(), std::numeric_limits<double>::infinity());
	std::vector<int>previous_pos(mesh.size(), -2);

	//DNodeは(スタート構造からその構造までのバリア値、(直前の構造の添字、その構造の添字))とする。
	typedef std::pair<double, std::pair<int, int>> DNode;

	std::priority_queue<DNode, std::vector<DNode>, std::greater<DNode>>dijkstra;
	dijkstra.push(std::make_pair(EnergyOfStructure(sequence, structure1), std::make_pair(-1, start_index)));
	while (!dijkstra.empty()) {
		const auto x = dijkstra.top();
		dijkstra.pop();
		const int index = x.second.second;
		if (previous_pos[index] != -2)continue;
		result_barrier[index] = x.first;
		previous_pos[index] = x.second.first;
		if (index == goal_index)break;
		for (int i = 0; i < mesh.size(); ++i)if (index != i) {
			if (previous_pos[i] != -2)continue;
			if (energy[i] > upper_bound)continue;
			if (i != goal_index && abs(mesh[index].first.first - mesh[i].first.first) + abs(mesh[index].first.second - mesh[i].first.second) > maxdist)continue;
			const double direct = GetDirectBarrier(index, i);
			const double new_result = std::max<double>(result_barrier[index], direct);
			if (new_result > upper_bound)continue;
			if (i == goal_index) {
				upper_bound = std::min<double>(upper_bound, new_result);
			}
			dijkstra.push(std::make_pair(new_result, std::make_pair(index, i)));
		}
	}

	//トレースバック
	std::vector<std::string>reversed_answer{ structure2 };
	for (int structure_index = goal_index; structure_index != start_index; structure_index = previous_pos[structure_index]) {
		const auto partial_pathway = Flamm2001FindpathDirect(sequence, mesh[previous_pos[structure_index]].second.second, mesh[structure_index].second.second, 10).first;
		for (int i = partial_pathway.size() - 2; i >= 0; --i)reversed_answer.push_back(partial_pathway[i]);
	}
	assert(reversed_answer.back() == structure1);
	reversed_answer = TrimPathway(reversed_answer);
	std::reverse(reversed_answer.begin(), reversed_answer.end());
	return std::make_pair(reversed_answer, BarrierEnergy(sequence, reversed_answer));
}

}