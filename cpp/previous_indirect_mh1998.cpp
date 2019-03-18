/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "previous_indirect_mh1998.h"

#include "call_vienna_rna.h"
#include "misc.h"
#include "enumerate_long_stacks.h"
#include "previous_direct_methods.h"

namespace czno_cpp {

std::vector<std::string>DistinctSamples(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int num) {
	const auto samples = SampleStructure(sequence, num);
	std::set<std::string>distinct_samples(samples.begin(), samples.end());
	distinct_samples.erase(structure1);
	distinct_samples.erase(structure2);
	std::vector<std::string>distinct_samples_vector(distinct_samples.begin(), distinct_samples.end());
	distinct_samples_vector.push_back(structure1);
	distinct_samples_vector.push_back(structure2);
	return distinct_samples_vector;
}

std::pair<std::vector<std::string>, double>
MorganHiggs1998Indirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int num) {

	const auto samples = DistinctSamples(sequence, structure1, structure2, num);
	std::vector<std::vector<double>>barrier(samples.size(), std::vector<double>(samples.size(), std::numeric_limits<double>::infinity()));
	for (int i = 0; i < samples.size(); ++i) {
		barrier[i][i] = EnergyOfStructure(sequence, samples[i]);
		for (int j = i + 1; j < samples.size(); ++j) {
			barrier[i][j] = barrier[j][i] = MorganHiggs1998Direct(sequence, samples[i], samples[j]).second;
		}
	}

	std::vector<double>result_barrier(samples.size(), std::numeric_limits<double>::infinity());
	std::vector<int>previous_pos(samples.size(), -2);

	//DNodeは(スタート構造からその構造までのバリア値、(直前の構造の添字、その構造の添字))とする。
	typedef std::pair<double, std::pair<int, int>> DNode;

	std::priority_queue<DNode, std::vector<DNode>, std::greater<DNode>>dijkstra;
	dijkstra.push(std::make_pair(barrier[samples.size() - 2][samples.size() - 2], std::make_pair(-1, samples.size() - 2)));
	while (!dijkstra.empty()) {
		const auto x = dijkstra.top();
		dijkstra.pop();
		const int index = x.second.second;
		if (previous_pos[index] != -2)continue;
		result_barrier[index] = x.first;
		previous_pos[index] = x.second.first;
		if (index == samples.size() - 1)break;
		for (int i = 0; i < samples.size(); ++i)if (index != i) {
			if (previous_pos[i] != -2)continue;
			const double new_result = std::max<double>(result_barrier[index], barrier[index][i]);
			if (new_result >= result_barrier[i])continue;
			result_barrier[i] = new_result;
			dijkstra.push(std::make_pair(new_result, std::make_pair(index, i)));
		}
	}

	//トレースバック
	std::vector<std::string>reversed_answer{ structure2 };
	for (int structure_index = samples.size() - 1; structure_index != samples.size() - 2; structure_index = previous_pos[structure_index]) {
		const auto partial_pathway = MorganHiggs1998Direct(sequence, samples[previous_pos[structure_index]], samples[structure_index]).first;
		for (int i = partial_pathway.size() - 2; i >= 0; --i)reversed_answer.push_back(partial_pathway[i]);
	}
	assert(reversed_answer.back() == structure1);
	reversed_answer = TrimPathway(reversed_answer);
	std::reverse(reversed_answer.begin(), reversed_answer.end());
	return std::make_pair(reversed_answer, BarrierEnergy(sequence, reversed_answer));
}

}