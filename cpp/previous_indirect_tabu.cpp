/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "previous_indirect_tabu.h"

#include "call_vienna_rna.h"
#include "misc.h"
#include "enumerate_long_stacks.h"

namespace czno_cpp {

static std::vector<std::string> EnumerateNeighborsWithoutTabu(
	const std::string& sequence,
	const std::string& structure1,
	const int count,
	const std::vector<std::vector<int>>& tabu) {

	const int n = sequence.length();
	const auto A = DotNotationToBasePairSet(structure1);
	std::vector<std::string> answer;
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			const std::pair<int, int>x = std::make_pair(i, j);
			if (!IsValidBasePair(sequence, x))continue;
			if (tabu[i][j] > count)continue;
			if (A.find(x) == A.end()) {
				bool flag = false;
				for (const auto y : A) {
					flag = IsExclusive(x, y);
					if (flag)break;
				}
				if (flag)continue;
				auto AA = A;
				AA.insert(x);
				answer.push_back(BasePairSetToDotNotation(AA, sequence.size()));
			}
			else {
				auto AA = A;
				assert(AA.erase(x) == 1);
				answer.push_back(BasePairSetToDotNotation(AA, sequence.size()));
			}
		}
	}
	return answer;
}

std::pair<std::vector<std::string>, double>
Dotu2010TabuIndirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int k,
	const double w_0,
	const int max_stable,
	const double energy_diff_bound,
	const int random_seed) {

	ValidateTriple(sequence, structure1, structure2);
	assert(1 <= k);
	assert(0.0 <= w_0);
	assert(0 <= max_stable);
	assert(0.0 <= energy_diff_bound);

	const int n = sequence.length();
	std::mt19937_64 rnd(random_seed);
	std::uniform_int_distribution<int>tabu_tenure(4, std::max<int>(4, HammingDistance(structure1, structure2) - 1));

	//変数たちの準備
	const double energy_upper_bound = std::max<double>(
		EnergyOfStructure(sequence, structure1),
		EnergyOfStructure(sequence, structure2)
		) + energy_diff_bound;
	double w_init = w_0;
	double w = w_init;
	auto S = structure1;
	double max_energy = EnergyOfStructure(sequence, S);
	std::vector<std::string>pathway{ structure1 };
	auto closest_structure = structure1;
	int closest_time = 1;
	std::vector<std::vector<int>>tabu(n, std::vector<int>(n, 0));
	int best_distance = HammingDistance(structure1, structure2);
	int no_improvement_count = 0;

	//Fitness関数は[Dotu et al.,2010]のFig10の下にある式
	const auto Fitness = [&](const std::string& x) {
		return EnergyOfStructure(sequence, x) + w * HammingDistance(x, structure2);
	};

	for (int count = 0; count < 1000 && S != structure2; ++count) {

		/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));
		std::cout << "LOG: RNATABUPATH: start: count = " << count << std::endl;
		//Sからハミング距離1違う構造たちのうち、エネルギーが一定値以下のものを
		//全てフィットネス関数の値と合わせて格納する。
		std::vector<std::pair<double, std::string>>low_energy_neighbors;
		for (const auto x : EnumerateNeighborsWithoutTabu(sequence, S, count, tabu)) {
			if (EnergyOfStructure(sequence, x) <= energy_upper_bound) {
				low_energy_neighbors.push_back(std::make_pair(Fitness(x), x));
			}
		}

		//遷移可能な構造が存在しない場合、tabu塩基対が復活するまで待つ。
		if (low_energy_neighbors.size() == 0)continue;

		//フィットネス関数の評価値が低いk個のなかからランダムに選び遷移する。
		std::sort(low_energy_neighbors.begin(), low_energy_neighbors.end());
		std::uniform_int_distribution<int>dist(0, std::min<int>(k, low_energy_neighbors.size()) - 1);
		const auto nextS = low_energy_neighbors[dist(rnd)].second;
		const auto diff = Difference(S, nextS);
		assert(tabu[diff.first][diff.second] <= count);
		S = nextS;
		pathway.push_back(S);
		tabu[diff.first][diff.second] = count + tabu_tenure(rnd);//先行研究の実装に合わせた。
		max_energy = std::max<double>(max_energy, EnergyOfStructure(sequence, S));
		/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));

		//今までで一番ゴールに近くなったかどうかで場合分け
		const int hamming_distance = HammingDistance(S, structure2);
		if (hamming_distance < best_distance) {
			best_distance = hamming_distance;
			closest_structure = S;
			closest_time = pathway.size();
			w = std::max<double>(w_0, w / 1.2);//先行研究の実装に合わせた。
			no_improvement_count = 0;
		}
		else {
			no_improvement_count++;
			if (no_improvement_count > max_stable) {

				//一定回数停滞したら、一番ゴールに近かった構造に巻き戻す。

				S = closest_structure;
				pathway.resize(closest_time);
				/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));
				for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)tabu[i][j] = 0;
				w *= 2;//先行研究の実装に合わせた。
				no_improvement_count = 0;

				if (w > 10.0 * w_init) {

					//停滞して巻き戻すが改善しないというのを一定回数繰り返したら、
					//最初の構造に巻き戻す。

					S = structure1;
					pathway.resize(1);
					for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)tabu[i][j] = 0;
					w_init += 0.5;//先行研究の実装に合わせた。
					w = w_init;

					//以下の3行は論文には書かれていないがDotu et al.,の実装にはある。必須。

					closest_structure = structure1;
					closest_time = 1;
					best_distance = HammingDistance(structure1, structure2);
				}
			}
		}
	}

	if (S != structure2) {
		//energy_upper_boundが低すぎたか、
		//w_0が小さすぎて収束する前にcountがオーバーした。
		return make_pair(std::vector<std::string>{}, std::numeric_limits<double>::infinity());
	}

	pathway = TrimPathway(pathway);
	return make_pair(pathway, BarrierEnergy(sequence, pathway));
}

}