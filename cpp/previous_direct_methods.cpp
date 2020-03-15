/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include "previous_direct_methods.h"

#include "call_vienna_rna.h"
#include "misc.h"

namespace czno_cpp {

std::pair<std::vector<std::string>, double>
MorganHiggs1998Direct(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2) {
	//Morgan S, Higgs P. Barrier heights between ground states in a model of RNA secondary structure. J. Phys. A: Math. Gen.  (1998) 31:3153–3170.

	ValidateTriple(sequence, structure1, structure2);
	const int n = sequence.length();
	std::vector<std::string>pathway{ structure1 };
	const auto B = DotNotationToBasePairSet(structure2);
	while (pathway.back() != structure2) {
		auto A = DotNotationToBasePairSet(pathway.back());

		std::vector<std::pair<int, int>>need_to_close;
		std::vector<std::pair<int, int>>need_to_open;
		for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
		for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);

		//同じ入力に対して常に同じ出力を返したいのでソートする。
		std::sort(need_to_close.begin(), need_to_close.end());
		std::sort(need_to_open.begin(), need_to_open.end());

		if (need_to_close.size() >= 1) {

			//clash[A](x)は、二次構造Aに塩基対xを組ませるために予め外さなければならない塩基対の集合を返す。
			//xは別の二次構造Bのものだとすると、そのような塩基対はAで組んでいてBで組んでいないもの、
			//すなわちneed_to_openの部分集合だと言える。
			const auto clash = [&](const std::pair<int, int>& x) {
				std::vector<std::pair<int, int>> clash_list;
				for (const auto y : need_to_open)if (IsExclusive(x, y))clash_list.push_back(y);
				return clash_list;
			};

			//clash関数の返り値の要素数が最小になるxをnext_closeとし、
			//next_opens := clash(next_close) とする。
			std::pair<int, int>next_close = need_to_close[0];
			std::vector<std::pair<int, int>>next_opens = clash(next_close);
			for (int i = 1; i < need_to_close.size(); ++i) {
				std::vector<std::pair<int, int>> x = clash(need_to_close[i]);
				if (x.size() < next_opens.size()) {
					next_opens = x;
					next_close = need_to_close[i];
				}
			}

			//next_opensを外してからnext_closeを組む。
			for (const auto x : next_opens) {
				assert(A.erase(x) == 1);
				pathway.push_back(BasePairSetToDotNotation(A, sequence.size()));
			}
			assert(A.find(next_close) == A.end());
			A.insert(next_close);
			pathway.push_back(BasePairSetToDotNotation(A, sequence.size()));
		}
		else {
			for (const auto x : need_to_open) {
				assert(A.erase(x) == 1);
				pathway.push_back(BasePairSetToDotNotation(A, sequence.size()));
			}
		}
	}

	return std::make_pair(pathway, BarrierEnergy(sequence, pathway));
}

std::pair<std::vector<std::string>, double>
Voss2004Direct(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int k,
	const int random_seed) {
	//Voss B, Meyer C, Giegerich R. Evaluating the predictability of conformational switching in RNA. Bioinformatics  (2004) 20:1573–1582.
	//Greedy(k=1) and Semi-Greedy(k>=2).

	ValidateTriple(sequence, structure1, structure2);
	const int n = sequence.length();
	std::mt19937_64 rnd(random_seed);
	std::vector<std::string>pathway{ structure1 };
	const auto B = DotNotationToBasePairSet(structure2);
	while (pathway.back() != structure2) {
		auto A = DotNotationToBasePairSet(pathway.back());
		std::vector<std::pair<int, int>>need_to_close;
		std::vector<std::pair<int, int>>need_to_open;
		for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
		for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);

		//IsClashExist[A](x)は、二次構造Aに塩基対xを組ませられないときにtrueを返す。
		const auto IsClashExist = [&](const std::pair<int, int>& x) {
			for (const auto y : need_to_open)if (IsExclusive(x, y))return true;
			return false;
		};

		//Aに対するimmediate neighborであって、Bに近づくようなものについて、
		//(エネルギー、構造)のペアを作る。
		std::vector<std::pair<double, std::string>>next_energy_and_state;
		for (const auto x : need_to_close)if (!IsClashExist(x)) {
			auto AA = A;
			assert(AA.find(x) == AA.end());
			AA.insert(x);
			const auto ss = BasePairSetToDotNotation(AA, sequence.size());
			next_energy_and_state.push_back(std::make_pair(EnergyOfStructure(sequence, ss), ss));
		}
		for (const auto x : need_to_open) {
			auto AA = A;
			assert(AA.erase(x) == 1);
			const auto ss = BasePairSetToDotNotation(AA, sequence.size());
			next_energy_and_state.push_back(std::make_pair(EnergyOfStructure(sequence, ss), ss));
		}

		//エネルギーが最も小さいk個の中からランダムに選ぶ。
		std::sort(next_energy_and_state.begin(), next_energy_and_state.end());
		std::uniform_int_distribution<int>dist(0, std::min<int>(k, next_energy_and_state.size()) - 1);
		const auto next = next_energy_and_state[dist(rnd)];
		pathway.push_back(next.second);
	}

	return std::make_pair(pathway, BarrierEnergy(sequence, pathway));
}

std::pair<std::vector<std::string>, double>
Flamm2001FindpathDirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int maxkeep) {
	//Flamm C, Hofacker IL, Maurer-Stroh S, Stadler PF, Zehl M. Design of multistable RNA molecules. RNA  (2001) 7:254–265.

	ValidateTriple(sequence, structure1, structure2);
	const int n = sequence.length();
	const std::set<std::pair<int, int>> B = DotNotationToBasePairSet(structure2);

	//EnumerateApproachingNeighbors[B](A)は、
	//Aに対するimmediate neighborであって、Bに近づくようなものをすべて求める。
	const auto EnumerateApproachingNeighbors = [&](const std::string now_structure) {
		const std::set<std::pair<int, int>> A = DotNotationToBasePairSet(now_structure);
		std::vector<std::pair<int, int>>need_to_close;
		std::vector<std::pair<int, int>>need_to_open;
		for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
		for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);
		const auto IsClashExist = [&](const std::pair<int, int>& x) {
			for (const auto y : need_to_open)if (IsExclusive(x, y))return true;
			return false;
		};
		std::vector<std::string>answer;
		for (const auto x : need_to_close)if (!IsClashExist(x)) {
			auto AA = A;
			assert(AA.find(x) == AA.end());
			AA.insert(x);
			answer.push_back(BasePairSetToDotNotation(AA, sequence.size()));
		}
		for (const auto x : need_to_open) {
			auto AA = A;
			assert(AA.erase(x) == 1);
			answer.push_back(BasePairSetToDotNotation(AA, sequence.size()));
		}
		return answer;
	};

	const double start_energy = EnergyOfStructure(sequence, structure1);

	//(double,double)は(バリアエネルギー、現在のエネルギー)
	//(int, std::string)は(直前の構造の添字、現在の構造)
	std::vector<std::vector<std::pair<std::pair<double, double>, std::pair<int, std::string>>>>beam_search{
		std::vector<std::pair<std::pair<double, double>, std::pair<int, std::string>>>{
		std::make_pair(std::make_pair(start_energy, start_energy), std::make_pair(-1, structure1))} };

	while (beam_search.back()[0].second.second != structure2) {

		//現在のビームの先端の各構造から見たapproaching neighborをkeyとして、
		//((スタートから次構造までのバリア、次構造のエネルギー)、次構造の直前の構造の添字)をvalueとするmapを作る。
		//現在のビームの先端の複数構造から同じkeyが出ることもあるが、そのときは最小のvalueだけを持つ。
		std::map<std::string, std::pair<std::pair<double, double>,int>>next_state_and_energy;
		for (int i = 0; i < maxkeep && i < beam_search.back().size(); ++i) {
			const auto next = EnumerateApproachingNeighbors(beam_search.back()[i].second.second);
			for (int j = 0; j < next.size(); ++j) {
				const auto x = next[j];
				assert(HammingDistance(x, beam_search.back()[i].second.second) == 1);
				if (next_state_and_energy.find(x) != next_state_and_energy.end())continue;
				const double now_energy = EnergyOfStructure(sequence, x);
				const double barrier_energy = std::max<double>(now_energy, beam_search.back()[i].first.first);
				const auto energy_pair = std::make_pair(std::make_pair(barrier_energy, now_energy), i);
				if (next_state_and_energy.find(x) == next_state_and_energy.end() || energy_pair < next_state_and_energy[x]) {
					next_state_and_energy[x] = energy_pair;
				}
			}
		}

		//(スタートから次構造までのバリア、次構造のエネルギー、次構造の直前の構造の添字、次構造)を並べてソートし、小さい順にk個取る。(ビームサーチ)
		std::vector<std::pair<std::pair<double, double>, std::pair<int, std::string>>>next_energy_and_state;
		for (const auto x : next_state_and_energy)next_energy_and_state.push_back(std::make_pair(x.second.first, std::make_pair(x.second.second, x.first)));
		std::sort(next_energy_and_state.begin(), next_energy_and_state.end());
		beam_search.push_back(std::vector<std::pair<std::pair<double, double>, std::pair<int, std::string>>>());
		for (int i = 0; i < maxkeep && i < next_energy_and_state.size(); ++i) {
			beam_search.back().push_back(next_energy_and_state[i]);
		}
	}

	for (int i = 1; i < beam_search.size(); ++i) {
		for (int j = 0; j < beam_search[i].size(); ++j) {
			assert(HammingDistance(beam_search[i][j].second.second, beam_search[i - 1][beam_search[i][j].second.first].second.second) == 1);
		}
	}

	std::vector<std::string>reverse_trajectory;
	int index = 0;
	for (int i = beam_search.size() - 1; i >= 0; --i) {
		reverse_trajectory.push_back(beam_search[i][index].second.second);
		index = beam_search[i][index].second.first;
	}
	std::reverse(reverse_trajectory.begin(), reverse_trajectory.end());
	return make_pair(reverse_trajectory, beam_search.back()[0].first.first);
}

}