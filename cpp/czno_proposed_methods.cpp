/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include "czno_proposed_methods.h"

#include "call_vienna_rna.h"
#include "misc.h"

namespace czno_cpp {

template<bool fine_grained>std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra1(
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
	if (hamming_distance == 0) {
		return std::make_pair(std::vector<std::string>{structure1}, EnergyOfStructure(sequence, structure1));
	}
	if (hamming_distance == 1) {
		return std::make_pair(
			std::vector<std::string>{structure1, structure2},
			std::max<double>(
				EnergyOfStructure(sequence, structure1),
				EnergyOfStructure(sequence, structure1)
				)
		);
	}

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

	std::array<std::vector<int>, 2>need_to_open{ std::vector<int>(hamming_distance, 0),std::vector<int>(hamming_distance, 0) };
	for (int i = 0; i < hamming_distance; ++i) {
		need_to_open[0][i] = (A.find(basepair_difference[i]) != A.end()) ? 1 : 0;
		need_to_open[1][i] = (B.find(basepair_difference[i]) != B.end()) ? 1 : 0;
	}

	const auto DotNotation = [&](const int universal_structure_index) {
		std::string answer = sharing_structure;
		for (int i = 0; i < hamming_distance; ++i) {
			if (((universal_structure_index & (1 << i)) ? 1 : 0) ^ need_to_open[0][i]) {
				assert(basepair_difference[i].first < basepair_difference[i].second);
				assert(answer[basepair_difference[i].first] == '.');
				assert(answer[basepair_difference[i].second] == '.');
				answer[basepair_difference[i].first] = '(';
				answer[basepair_difference[i].second] = ')';
			}
		}
		return answer;
	};

	std::array<std::vector<int>, 2>prerequisite_flag{ std::vector<int>(hamming_distance, 0),std::vector<int>(hamming_distance, 0) };
	for (int i = 0; i < hamming_distance; ++i) {
		for (int d = 0; d < 2; ++d) {
			if (need_to_open[d][i] == 0)continue;
			for (int j = 0; j < hamming_distance; ++j) {
				if (need_to_open[d][j] != 0)continue;

				//この時点で、
				//basepair_difference[d][i]はstartにあってgoalにない、すなわちパス上で外されるべき塩基対であり、
				//basepair_difference[d][j]はgoalにあってstartにない、すなわちパス上で組まれるべき塩基対である。
				//ここでもし[i]を外してからでないと[j]を組めないという依存関係があるなら、
				//prerequisite_flag[d][j] |= 1 << i とする。

				if (IsExclusive(basepair_difference[i], basepair_difference[j])) {
					prerequisite_flag[d][j] |= 1 << i;
				}
			}
		}
	}

	//スタート構造から構造xまでの部分解。
	std::array<std::vector<double>, 2>result{ std::vector<double>(1 << hamming_distance, std::numeric_limits<double>::infinity()),std::vector<double>(1 << hamming_distance, std::numeric_limits<double>::infinity()) };

	//構造xのエネルギーの値
	std::vector<double>energy(1 << hamming_distance, std::numeric_limits<double>::infinity());

	//Resultが確定したかどうか
	std::vector<int>searched(1 << hamming_distance, 0);

	//DNodeは(スタート構造からその構造までのバリア値、(-時刻、その構造の添字))とする。
	//"-時刻"を入れる理由は、バリア値が同じ場合に後入れ先出しにしたいから。
	//構造の添字とは、すべての可能な中間構造をBitDP的に添字付けしたものとする。
	typedef std::pair<double, std::pair<int, int>> DNode;
	typedef std::priority_queue<DNode, std::vector<DNode>, std::greater<DNode>> DQueue;

	std::array<DQueue, 2> dijkstra;
	int count = 0;
	energy[0] = EnergyOfStructure(sequence, structure1);
	energy[(1 << hamming_distance) - 1] = EnergyOfStructure(sequence, structure2);
	result[0][0] = energy[0];
	result[1][(1 << hamming_distance) - 1] = energy[(1 << hamming_distance) - 1];
	dijkstra[0].push(std::make_pair(energy[0], std::make_pair(count, 0)));
	dijkstra[1].push(std::make_pair(energy[(1 << hamming_distance) - 1], std::make_pair(count, (1 << hamming_distance) - 1)));
	searched[0] = 1 << 0;
	searched[(1 << hamming_distance) - 1] = 1 << 1;

	double lower_bound = -std::numeric_limits<double>::infinity();
	double upper_bound = std::numeric_limits<double>::infinity();
	int ub_universal_position = -1;

	//fine_grained版においてはコンテキストスイッチするため全てのローカル変数をここで宣言する必要がある。
	std::array<DNode, 2>x_;
	std::array<int, 2>b_;
	std::array<int, 2>i_;
	std::array<int, 2>f_;
	std::array<int, 2>directional_structure_index_;
	std::array<int, 2>universal_structure_index_;
	std::array<int, 2>new_directional_structure_index_;
	std::array<int, 2>new_universal_structure_index_;
	std::array<double, 2>barrier_if_passing_here_;
	std::array<int, 2>context_{ 0,0 };

	//const auto DijkstraOneStep = [&]<bool direction>() {
	const auto DijkstraOneStep = [&](const int direction) {
		if (fine_grained && context_[direction])goto LABEL;
		count--;
		x_[direction] = dijkstra[direction].top();
		dijkstra[direction].pop();
		universal_structure_index_[direction] = x_[direction].second.second;
		directional_structure_index_[direction] = universal_structure_index_[direction] ^ (((1 << hamming_distance) - 1) * direction);
		assert(x_[direction].first <= upper_bound);

		//真のバリアエネルギーは、priority_queueから取り出した値以上である。
		//これはどちらのpriority_queueから取り出した値であっても常に成り立つ。
		//理由はpriority_queueから出てくる順番がバリア値の低い順であることとmin-max代数であることによる。
		lower_bound = std::max<double>(lower_bound, x_[direction].first);

		//真のバリアエネルギーがupper_boundに等しいことが確定したら、あとはub_universal_positionから両側にトレースバックすればよい。
		if (lower_bound == upper_bound)return 1;

		//↓のif文が成立する状況においては、必ず↑の(lower_bound == upper_bound)が成立しているため、↓は不要である。
		//なぜなら、現在のノードに関して、
		//(1)両者がこれを発見してupper_boundが更新されるタイミング　と
		//(2)どちらか片方がこれをpriority_queueから取り出してlower_boundが更新されるタイミング　と
		//(3)両者がこれをpriority_queueから取り出すタイミング　とでは
		//必ず(3)が一番遅い。遅くとも(1)と(2)がなされた時点で真の解が確定して探索終了するが、↓は(3)の状況を意味している。
		//ちなみに、「遅くとも」と書いたが、upper_boundとlower_boundが別々の地点で確定してもよいからである。
		//if (searched[universal_structure_index] & (1 << (1 - direction)))return 0;

		searched[universal_structure_index_[direction]] |= (1 << direction);
		for (i_[direction] = 0; i_[direction] < hamming_distance; ++i_[direction]) {
			b_[direction] = 1 << i_[direction];
			if (directional_structure_index_[direction] & b_[direction])continue;
			if ((directional_structure_index_[direction] & prerequisite_flag[direction][i_[direction]]) != prerequisite_flag[direction][i_[direction]])continue;
			new_directional_structure_index_[direction] = directional_structure_index_[direction] | b_[direction];
			new_universal_structure_index_[direction] = new_directional_structure_index_[direction] ^ (((1 << hamming_distance) - 1) * direction);

			//自分の出発点からあるノードに行くのに必要なバリアは、そのノードを最初に発見した時点で確定する。
			//理由はDijkstra法であることとmin-max代数であることにより自明である。
			//ゆえに、自分がかつて発見したことのあるノードを再発見した場合、continueしてもよい。
			if (result[direction][new_universal_structure_index_[direction]] != std::numeric_limits<double>::infinity())continue;

			if (fine_grained)f_[direction] = 0;
			if (energy[new_universal_structure_index_[direction]] == std::numeric_limits<double>::infinity()) {
				energy[new_universal_structure_index_[direction]] = EnergyOfStructure(sequence, DotNotation(new_universal_structure_index_[direction]));
				if (fine_grained)f_[direction] = 1;
			}
			result[direction][new_universal_structure_index_[direction]] = std::max<double>(result[direction][universal_structure_index_[direction]], energy[new_universal_structure_index_[direction]]);
			if (upper_bound < result[direction][new_universal_structure_index_[direction]]) {
				if (fine_grained && f_[direction]) {
					context_[direction] = 1;
					return 0;
				}
				else continue;
			}
			if (result[1 - direction][new_universal_structure_index_[direction]] != std::numeric_limits<double>::infinity()) {

				//一方の側から探索していって、逆側が発見済みのノードを発見したとする。
				//自分のresultのほうが低い場合：
				//　　相手もDijkstra法であることから、今発見したノードから相手の出発点に行くのに必要なバリアは相手のresultで間違いない。
				//相手のresultのほうが低い場合：
				//　　自分はDijkstra法であることから、今発見したノードから自分の出発点に行くのに必要なバリアは自分のresultで間違いない。
				//結局、発見したノードを通る場合のバリアエネルギーは、今発見したノードにおける自分のresultと相手のresultのどちらか悪いほうで間違いない。

				barrier_if_passing_here_[direction] = std::max<double>(result[direction][new_universal_structure_index_[direction]], result[1 - direction][new_universal_structure_index_[direction]]);
				if (barrier_if_passing_here_[direction] < upper_bound) {
					upper_bound = barrier_if_passing_here_[direction];
					ub_universal_position = new_universal_structure_index_[direction];
					if (lower_bound == upper_bound)return 1;
				}
			}
			dijkstra[direction].push(std::make_pair(result[direction][new_universal_structure_index_[direction]], std::make_pair(count, new_universal_structure_index_[direction])));
			if (fine_grained && f_[direction]) {
				context_[direction] = 1;
				return 0;
			}
		LABEL:;
		}
		if (fine_grained)context_[direction] = 0;
		return 0;
	};

	const auto DijkstraExec = [&]() {
		while (true) {
			if (DijkstraOneStep(0))return;
			if (DijkstraOneStep(1))return;
		}
	};

	DijkstraExec();
	assert(ub_universal_position != -1);
	const int meeting_position = ub_universal_position;
	const std::string meeting_structure = DotNotation(meeting_position);

	//会合点からA・Bへのトレースバック
	std::vector<std::vector<std::string>>trajectory(2, std::vector<std::string>{ meeting_structure });
	const auto TraceBack = [&](const int direction) {
		for (int universal_structure_index = meeting_position; universal_structure_index != ((1 << hamming_distance) - 1) * direction;) {
			const int directional_structure_index = universal_structure_index ^ (((1 << hamming_distance) - 1) * direction);
			bool flag = false;
			for (int i = 0; i < hamming_distance; ++i) {
				const int b = 1 << i;
				if ((directional_structure_index & b) == 0)continue;
				const int back_directional_structure_index = directional_structure_index ^ b;
				const int back_universal_structure_index = back_directional_structure_index ^ (((1 << hamming_distance) - 1) * direction);
				if ((back_directional_structure_index & prerequisite_flag[direction][i]) != prerequisite_flag[direction][i])continue;
				if ((searched[back_universal_structure_index] & (1 << direction)) == 0)continue;
				trajectory[direction].push_back(DotNotation(back_universal_structure_index));
				universal_structure_index = back_universal_structure_index;
				flag = true;
				break;
			}
			assert(flag);
		}
	};
	TraceBack(0);
	TraceBack(1);
	std::reverse(trajectory[0].begin(), trajectory[0].end());
	trajectory[0].pop_back();
	for (const auto x : trajectory[1])trajectory[0].push_back(x);
	return std::make_pair(trajectory[0], upper_bound);
}

std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra1Turn(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2) {

	return MinimumBarrierDirectPathDijkstra1<true>(sequence, structure1, structure2);
}


std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstra(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2) {

	return MinimumBarrierDirectPathDijkstra1<false>(sequence, structure1, structure2);
}


std::pair<std::vector<std::string>, double>
MinimumBarrierDirectPathDijkstraOld(
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
