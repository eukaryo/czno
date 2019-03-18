/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "previous_indirect_ea.h"

#include "call_vienna_rna.h"
#include "misc.h"
#include "enumerate_long_stacks.h"

namespace czno_cpp {

std::vector<std::pair<int, int>>GenerateRandomSimplePathway(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	std::mt19937_64& rnd) {

	std::vector<std::pair<int, int>>need_to_close;
	std::vector<std::pair<int, int>>need_to_open;
	for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
	for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);
	std::sort(need_to_close.begin(), need_to_close.end());
	std::sort(need_to_open.begin(), need_to_open.end());
	std::shuffle(need_to_close.begin(), need_to_close.end(), rnd);
	std::shuffle(need_to_open.begin(), need_to_open.end(), rnd);
	std::vector<std::pair<int, int>>answer;
	auto S = A;
	for (const auto x : need_to_open) {
		assert(S.erase(x) == 1);
		answer.push_back(x);
	}
	for (const auto x : need_to_close) {
		assert(S.find(x) == S.end());
		S.insert(x);
		answer.push_back(x);
	}
	return answer;
}

std::set<std::pair<int, int>>AssosiatedIntermediateStructure(
	const std::set<std::pair<int, int>>& A,
	const std::vector<std::pair<int, int>>& action_list,
	const int t1) {

	auto S = A;
	for (int i = 0; i <= t1; ++i) {
		if (S.find(action_list[i]) == S.end()) S.insert(action_list[i]);
		else S.erase(action_list[i]);
	}

	return S;
}

std::vector<std::pair<int, int>>M1Mechanism(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	const int t1,
	const int t2) {
	//action_list[t1]を除去して、それをaction_list[t2]の直後に挿入する。

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2 && t2 < action_list.size());
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	std::vector<std::pair<int, int>>new_action_chain;
	if (t1 < t2) {
		//[t1]を抜いて[t2]の直後に入れる。
		for (int i = 0; i < t1; ++i)new_action_chain.push_back(action_list[i]);
		for (int i = t1+1; i <= t2; ++i)new_action_chain.push_back(action_list[i]);
		new_action_chain.push_back(action_list[t1]);
		for (int i = t2 + 1; i < action_list.size(); ++i)new_action_chain.push_back(action_list[i]);
	}
	else if(t2 + 1 < t1) {
		//[t1]を抜いて[t2]の直後に入れる。
		for (int i = 0; i <= t2; ++i)new_action_chain.push_back(action_list[i]);
		new_action_chain.push_back(action_list[t1]);
		for (int i = t2 + 1; i < t1; ++i)new_action_chain.push_back(action_list[i]);
		for (int i = t1 + 1; i < action_list.size(); ++i)new_action_chain.push_back(action_list[i]);
	}
	else new_action_chain = action_list;//抜いた場所に入れる、つまり何も変わらない
	return new_action_chain;
}

std::vector<std::pair<int, int>>M2Mechanism(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	std::vector<std::pair<int, int>> action_list,
	const int t1,
	const int t2) {
	//[t1]と[t2]を入れ替える。

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2 && t2 < action_list.size());
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const auto tmp = action_list[t1];
	action_list[t1] = action_list[t2];
	action_list[t2] = tmp;
	return action_list;
}

std::vector<std::pair<int, int>>M3Mechanism(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>> action_list,
	const std::pair<int, int> x,
	const int t1,
	const int t2) {
	//action_list[t1]の直後とaction_list[t2]の直後にxを挿入する。

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2 && t2 < action_list.size());
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int tmin = t1 < t2 ? t1 : t2;
	const int tmax = t1 < t2 ? t2 : t1;
	std::vector<std::pair<int, int>>new_action_chain;
	for (int i = 0; i <= tmin; ++i)new_action_chain.push_back(action_list[i]);
	new_action_chain.push_back(x);
	for (int i = tmin + 1; i <= tmax; ++i)new_action_chain.push_back(action_list[i]);
	new_action_chain.push_back(x);
	for (int i = tmax + 1; i < action_list.size(); ++i)new_action_chain.push_back(action_list[i]);
	return new_action_chain;
}

int GaussianBiasedSelection(const int size, std::mt19937_64& rnd, const bool descending) {
	//descending: 0が出やすく、size-1は出にくい。
	//ascending: size-1が出やすく、0は出にくい。
	static std::normal_distribution<double>X(0.0, sqrt(1.0 / 12.0));
	const int index = std::min<int>(size - 1, int(std::abs(X(rnd)) * double(size)));
	return descending ? index : size - 1 - index;
}

std::vector<std::pair<int, int>>M1Strategy(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	const int t1,
	const int t2lb,
	const int t2ub,
	std::mt19937_64& rnd) {
	//action_list[t1]を除去し、それをaction_list[t2lb,t2ub]の範囲内のどれかの直後に挿入する。
	//action_listがvalidになるような挿入位置のうち、有望そうな挿入位置が選ばれやすいようなbiasをかける。
	//action_listがvalidになるような挿入位置が存在しないならば空っぽのpathwayを返す。

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2lb && t2lb <= t2ub && t2ub < action_list.size());
	assert(t2ub <= t1 || t1 <= t2lb);//t2の範囲はt1の前か後ろだけ
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	std::vector<std::vector<std::pair<int, int>>>valid_mutations;

	if (t2ub <= t1) {
		//t1より前にいるt2の直後にt1を挿入するという変異パスのうち、validなパスを全列挙する。
		for (int t2 = std::max<int>(0, t2lb); t2 <= t1 && t2 <= t2ub; ++t2) {
			const auto candidate = M1Mechanism(A, B, action_list, t1, t2);
			if (IsValidPathway(A, B, candidate))valid_mutations.push_back(candidate);
		}
	}
	else {
		//t1より後ろにいるt2の直後にt1を挿入するという変異パスのうち、validなパスを全列挙する。
		for (int t2 = std::min<int>(t1, t2lb); t2 < action_list.size() && t2 <= t2ub; ++t2) {
			const auto candidate = M1Mechanism(A, B, action_list, t1, t2);
			if (IsValidPathway(A, B, candidate))valid_mutations.push_back(candidate);
		}
	}

	//validなパスが無い場合は終了
	if (valid_mutations.size() == 0)return std::vector<std::pair<int, int>>{};

	//t1が塩基対を組む操作なのか外す操作なのかを調べる。
	int count = A.find(action_list[t1]) != A.end() ? 1 : 0;
	for (int i = 0; i < t1; ++i)if (action_list[i] == action_list[t1])count++;

	//この時点でcountが奇数ならt1は偶数回目であり外す処理である。
	//ゆえにこのときはt2はt1と近いほうがよい。t2のほうが遅いならDescendingBiasが欲しい。
	//t2が前の場合(M5で使う)、biasの向きは論文に書かれていない。よって安全策としてt2とt1が近いほうがよいことにする。
	//すなわち(t1 <= t2lb)がfalseなら引数は常にfalseとする。
	return valid_mutations[GaussianBiasedSelection(valid_mutations.size(), rnd, (count % 2 == 1) & (t1 <= t2lb))];

}

std::vector<std::pair<int, int>>M1(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M1をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	if (action_list.size() <= 1)return action_list;

	const int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 2)(rnd);
	return M1Strategy(A, B, action_list, t1, t1 + 1, action_list.size() - 1, rnd);
}

std::vector<std::pair<int, int>>M2(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M2をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//スワップしうるactionの組(i,j)を全列挙してシャッフルする。
	std::vector<std::pair<int, int>>candidates;
	for (int i = 0; i < action_list.size(); ++i) {
		for (int j = i + 1; j < action_list.size(); ++j) {
			candidates.push_back(std::make_pair(i, j));
		}
	}
	std::sort(candidates.begin(), candidates.end());//結果の再現性を保つためにソートする。
	std::shuffle(candidates.begin(), candidates.end(), rnd);

	for (const auto x : candidates) {
		const auto candidate = M2Mechanism(A, B, action_list, x.first, x.second);
		if (IsValidPathway(A, B, candidate))return candidate;
	}

	//どうシャッフルしてもinvalidになるなら終了。
	return std::vector<std::pair<int, int>>{};
}

std::vector<std::pair<int, int>>M3Strategy(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	const int t1,
	const int t2lb,
	const int t2ub,
	const std::pair<int,int>x,
	const bool is_plus,
	std::mt19937_64& rnd) {
	//action_list[t1]の直後にxを挿入し、かつ
	//action_list[t2lb,t2ub]の範囲内のどれかの直後にもxを挿入する。
	//action_listがvalidになるような挿入位置のうち、有望そうな挿入位置が選ばれやすいようなbiasをかける。
	//action_listがvalidになるような挿入位置が存在しないならば空っぽのpathwayを返す。

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2lb && t2lb <= t2ub && t2ub < action_list.size());
	assert(t2ub <= t1 || t1 <= t2lb);//t2の範囲はt1の前か後ろだけ
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//action_list[t1]の直後にxを挿入すること自体がinvalidな場合、
	//IsValidPathwayがすべてfalseになることで、結局空っぽのpathwayが返される。
	//よって、それをここで検査しなくてもよい。

	std::vector<std::vector<std::pair<int, int>>>valid_mutations;

	//3. ランダムに選んだ塩基対xに対して、それをt1の直後に組んだとして、
	//それを外すタイミングとしてvalidなタイミングt2を全列挙する。
	if (t2ub <= t1) {
		for (int t2 = std::max<int>(0, t2lb); t2 + 1 < t1 && t2 <= t2ub; ++t2) {
			const auto new_action_list = M3Mechanism(A, B, action_list, x, t1, t2);
			if (IsValidPathway(A, B, new_action_list))valid_mutations.push_back(new_action_list);
		}
	}
	else {
		//t1より後ろにいるt2の直後にt1を挿入するという変異パスのうち、validなパスを全列挙する。
		for (int t2 = std::min<int>(t1 + 1, t2lb); t2 < action_list.size() && t2 <= t2ub; ++t2) {
			const auto new_action_list = M3Mechanism(A, B, action_list, x, t1, t2);
			if (IsValidPathway(A, B, new_action_list))valid_mutations.push_back(new_action_list);
		}
	}

	//論文には書かれていないが、すべてのt2でinvalidなpathwayしか生じないことが考えられる。
	if (valid_mutations.size() == 0)return std::vector<std::pair<int, int>>{};
	
	//plus, t1<t2 → false //t1で組んだらできるだけ離れて外したい＋添字が大きいほど離れる
	//plus, t2<t1 → true  //t2で組んだらできるだけ離して外したい＋添字が小さいほど離れる
	//minus,t1<t2 → true  //t1で離したらできるだけ近くで組みたい＋添字が大きいほど離れる
	//minus,t2<t1 → false //t2で離したらできるだけ近くで組みたい＋添字が小きいほど離れる
	return valid_mutations[GaussianBiasedSelection(valid_mutations.size(), rnd, (t2ub <= t1) ^ is_plus)];
}

std::vector<std::pair<int, int>>M3Plus(
	const std::string& sequence, 
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M3+をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_listから一様に選び、t1とする。
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	const int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1の直後のstructure Sを塩基対集合として求める。引数のaction_listはvalidであることを仮定する。
	//メモ：
	//t1の直後に塩基対を挿入するという動作は論文の記述通りである。
	//>introducing an addition action add i,j after at1
	//論理的には、スタート構造の直後に塩基対を挿入しても問題ないはずだが、
	//論文の記述通りに実装したのでスタート構造の直後に塩基対が挿入されることはない。
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. 現在のstructure Sに対して、validかつcompatibleに組める塩基対(i,j)を全列挙する。
	std::vector<std::pair<int, int>> candidates;
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			const std::pair<int, int>x = std::make_pair(i, j);
			if (!IsValidBasePair(sequence, x))continue;
			if (S.find(x) == S.end()) {
				bool flag = false;
				for (const auto y : S) {
					flag = IsExclusive(x, y);
					if (flag)break;
				}
				if (flag)continue;
				candidates.push_back(x);
			}
		}
	}

	//validかつcompatibleに組める塩基対全てからランダムに1つ選ぶ。
	if (candidates.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto x = candidates[std::uniform_int_distribution<int>(0, candidates.size() - 1)(rnd)];

	return M3Strategy(sequence, A, B, action_list, t1, t1, action_list.size() - 1, x, true, rnd);
}

std::vector<std::pair<int, int>>M3Minus(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M3-をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_listから一様に選び、t1とする。
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	const int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1の直後のstructure Sを塩基対集合として求める。引数のaction_listはvalidであることを仮定する。
	//
	//t1の直後に塩基対を挿入するという動作は論文の記述通りである。
	//>introducing an addition action add i,j after at1
	//論理的には、スタート構造の直後に塩基対を挿入しても問題ないはずだが、
	//論文の記述通りに実装したのでスタート構造の直後に塩基対が挿入されることはない。
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. 現在のstructure Sに対して、validに外せる塩基対(i,j)を全列挙して、その中からランダムに選ぶ。
	std::vector<std::pair<int, int>>candidates(S.begin(), S.end());
	if (candidates.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto x = candidates[std::uniform_int_distribution<int>(0, candidates.size() - 1)(rnd)];

	return M3Strategy(sequence, A, B, action_list, t1, t1, action_list.size() - 1, x, false, rnd);
}

std::vector<std::pair<int, int>>M4Strategy(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	const int t1,
	const std::pair<int, int>x,
	std::mt19937_64& rnd) {
	//action_list[t1]よりも後方にxが出現するなら、M1のストラテジーでそれをt1の直後に移動する。
	//action_list[t1]よりも後方にxが出現しないなら、action_list[t1]の直後にxを挿入し、かつ
	//action_list[t2lb,t2ub]の範囲内のどれかの直後にもxを挿入する。
	//action_listがvalidになるような挿入位置のうち、有望そうな挿入位置が選ばれやすいようなbiasをかける。
	//action_listがvalidになるような挿入位置が存在しないならば空っぽのpathwayを返す。
	//action_list[t1]の直後にxが存在するとき、そのxは塩基対形成であり除去ではなく、かつそれ自体はvalidであるということを仮定する。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));
	
	//action_listの[t1]より後方に塩基対(i,j)が出てくるかどうか調べる。
	int same_action_pos = -1;
	for (int i = t1 + 1; i < action_list.size(); ++i) {
		if (action_list[i] == x) {
			same_action_pos = i;
			break;
		}
	}

	if (same_action_pos >= 0) {
		//3.1 action_listの[t1]より後方に塩基対(i,j)が出てくるなら、
		//それを"持ち上げて"t1の直後に持ってくる。M1のストラテジーで実行する。
		//>move it up and place it after a_t using strategy M1.
		return M1Strategy(A, B, action_list, same_action_pos, t1, t1, rnd);
	}
	//3.2 action_listの[t1]より後方に塩基対(i,j)が出てこないなら、
	//[t1]の直後に(i,j)を挿入し、かつさらに後方に(i,j)を除去する処理を挿入する。
	//それはM3のストラテジーで実行する。
	return M3Strategy(sequence, A, B, action_list, t1, t1, action_list.size() - 1, x, true, rnd);
}

std::vector<std::pair<int, int>>M4(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M4をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_listから一様に選び、t1とする。
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1の直後のstructure Sを塩基対集合として求める。引数のaction_listはvalidであることを仮定する。
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. 長さ4以上の連続した塩基対(stack)のうち、
	//二次構造Sに対して即座に追加できる(compatible)ようなものを全列挙し、ランダムに1つ選ぶ。
	//ちなみに長さ4は論文の値。
	//>all possible stacks with more than 3 consecutive base pairs
	auto stacks = EnumerateLongStacks(sequence, S, 4, true);
	std::sort(stacks.begin(), stacks.end());//再現性のためソートする。
	if (stacks.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto target_stack = stacks[std::uniform_int_distribution<int>(0, stacks.size() - 1)(rnd)];

	//3. 2.で選んだステムを、RNAの内側から順番に入れていく。
	auto new_action_list = action_list;
	for (int i = target_stack[0], j = target_stack[1]; i >= target_stack[2]; --i, ++j) {
		const auto stack_bp = std::make_pair(i, j);
		new_action_list = M4Strategy(sequence, A, B, new_action_list, t1, stack_bp, rnd);
		if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
		assert(new_action_list[t1 + 1] == stack_bp);
		++t1;

	}
	return new_action_list;
}

std::vector<std::pair<int, int>>M5(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//論文の変異生成手法M5をやる。失敗したら空っぽのpathwayを返す。

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//1. action_listのうち塩基対除去するものの添字たちから一様に選び、t1とする。
	std::vector<int>t1_candidate;
	{
		auto S = A;
		for (int i = 0; i < action_list.size(); ++i) {
			if (S.find(action_list[i]) != S.end()) {
				S.erase(action_list[i]);
				t1_candidate.push_back(i);
			}
			else {
				S.insert(action_list[i]);
			}
		}
	}
	if (t1_candidate.size() == 0)return std::vector<std::pair<int, int>>{};
	int t1 = t1_candidate[std::uniform_int_distribution<int>(0, t1_candidate.size() - 1)(rnd)];

	//t1の直後の二次構造を塩基対集合として求める。引数のaction_listはvalidであることを仮定する。
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. 長さ4以上の連続した塩基対(stack)のうち、
	//二次構造Sに対して即座に追加できない(incompatible)ようなものを全列挙し、ランダムに1つ選ぶ。
	//ちなみに長さ4は論文の値。
	//>all possible stacks with more than 3 consecutive base pairs
	auto stacks = EnumerateLongStacks(sequence, S, 4, false);
	std::sort(stacks.begin(), stacks.end());//再現性のためソートする。
	if (stacks.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto target_stack = stacks[std::uniform_int_distribution<int>(0, stacks.size() - 1)(rnd)];

	//target_stack中の塩基対たちを、compatibleかどうかで分類する。
	std::vector<std::pair<int, int>>compatible_pairs, incompatible_pairs;
	for (int i = target_stack[0], j = target_stack[1]; i >= target_stack[2]; --i, ++j) {
		const auto stack_bp = std::make_pair(i, j);
		if (S.find(stack_bp) != S.end())continue;
		bool incompatibility = false;
		for (const auto y : S) {
			incompatibility = IsExclusive(stack_bp, y);
			if (incompatibility)break;
		}
		if (incompatibility)incompatible_pairs.push_back(stack_bp);
		else compatible_pairs.push_back(stack_bp);
	}

	//target_stack中の塩基対たちを、RNAの内側から順番に入れていく。
	//compatibleなやつを全部入れてから、incompatibleなやつを入れていく。
	std::vector<std::pair<int, int>>new_action_list = action_list;
	for (const auto x : compatible_pairs) {

		//3.2 compatibleなものはM4のストラテジーで入れる。
		new_action_list = M4Strategy(sequence, A, B, new_action_list, t1, x, rnd);
		if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
		assert(new_action_list[t1 + 1] == x);
		++t1;
	}
	for (const auto x : incompatible_pairs) {

		//4.1 xに対してexclusiveなS上の塩基対それぞれについて対処する。
		for (const auto y : S)if (IsExclusive(x, y)) {

			//yがaction_list[t1]よりも後方に出てくるかで場合分けする。
			int same_action_pos = -1;
			for (int i = t1 + 1; i < new_action_list.size(); ++i) {
				if (new_action_list[i] == y) {
					same_action_pos = i;
					break;
				}
			}

			if (same_action_pos >= 0) {
				//4.2. yがaction_list[t1]よりも後方に出てくるなら、それを"持ち上げて"、t1より前方に挿入する。
				//M1のストラテジーを使う。
				new_action_list = M1Strategy(A, B, new_action_list, same_action_pos, 0, t1, rnd);
			}
			else {
				//4.3. yがaction_list[t1]よりも後方に出てこないなら、action_list[t1]の直後にyを挿入し、
				//さらに後方にもyを挿入する。M3のストラテジーを使う。
				new_action_list = M3Strategy(sequence, A, B, new_action_list, t1, t1, new_action_list.size() - 1, y, false, rnd);
			}
			if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
			++t1;
		}
		//この時点で、action_state[t1]の直後の中間構造に対してxがcompatibleである。
		//よって、3.2と同じくM4のストラテジーでxを入れられる。
		new_action_list = M4Strategy(sequence, A, B, new_action_list, t1, x, rnd);
		if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
		assert(new_action_list[t1 + 1] == x);
		++t1;
	}

	return new_action_list;
}

std::pair<std::vector<std::string>, double>
LiZhang2012RNAEAPathIndirect(
	const std::string& sequence,
	const std::string& structure1,
	const std::string& structure2,
	const int initial_population,//最初に適当に作るパスの数。4
	const int gamma_generation,//終了条件に関わるγ 5
	const int max_generation,//終了条件に関わるMAX 10
	const int script_capital_l,//子供の数に関わる、100
	const int l1_elite_population,//エリートの数、5
	const int l2_survive_population,// 5
	const int l3_final_max_population,// 100
	const int random_seed) {

	ValidateTriple(sequence, structure1, structure2);

	const int n = sequence.length();
	std::mt19937_64 rnd(random_seed);

	//Fig. 2の1行目。論文の記法だとmaxではなくabsだが、maxでも等価な処理はできる。
	const double delta_energy = std::max<double>(
		EnergyOfStructure(sequence, structure1),
		EnergyOfStructure(sequence, structure2));

	//simple pathwayを作るための準備
	const auto A = DotNotationToBasePairSet(structure1);
	const auto B = DotNotationToBasePairSet(structure2);
	std::vector<std::pair<int, int>>need_to_close;
	std::vector<std::pair<int, int>>need_to_open;
	for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
	for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);
	std::sort(need_to_close.begin(), need_to_close.end());
	std::sort(need_to_open.begin(), need_to_open.end());

	//Fig. 2の3行目。一番最後のintは、それがどのストラテジーで由来か
	std::vector<std::pair<std::pair<double, std::vector<std::pair<int, int>>>, int>>population;
	for (int i = 0; i < initial_population; ++i) {
		const auto pathway = GenerateRandomSimplePathway(A, B, rnd);
		population.push_back(std::make_pair(std::make_pair(BarrierEnergy(sequence, A, B, pathway), pathway), -1));
	}
	std::sort(population.begin(), population.end());

	//Fig. 2の4行目
	std::vector<std::pair<double, std::vector<std::pair<int, int>>>>best_pathways_history{ population[0].first };

	//kはFig. 2の2行目に出てくるやつ。 論文にあるStopping conditionsの判定を行う。
	const auto IsStoppingCondition = [&](const int k) {
		assert(best_pathways_history.size() == k + 1);
		if (best_pathways_history[k].first == delta_energy)return true;
		if (k >= gamma_generation &&
			best_pathways_history[k].first ==
			best_pathways_history[k - gamma_generation].first)return true;
		if (k >= max_generation &&
			best_pathways_history[k].first ==
			best_pathways_history[k - 1].first)return true;
		return false;
	};

	std::vector<int>generation_plan(5, script_capital_l / 5) ;

	//Fig. 2の5,6行目
	for (int k = 0; !IsStoppingCondition(k); ++k) {
		std::cout << "LOG: RNAEAPath: start: k = " << k << ", pop = " << population.size() << std::endl;

		//Fig. 2の7行目。エリートを入れる。
		auto next_population = population;
		next_population.resize(std::min<int>(next_population.size(), l1_elite_population));

		//Fig. 2の8行目
		for (const auto x : population) {

			//Fig. 9行目の左辺のTを定義する。
			std::vector<std::pair<std::pair<double, std::vector<std::pair<int, int>>>, int>>kindergarten;

			const auto RegistrataBirth = [&](const std::vector<std::pair<int, int>>& pathway, const int strategy_id) {
				if (pathway.size() == 0)return;
				/*HEAVY*/assert(IsValidPathway(A, B, pathway));
				kindergarten.push_back(std::make_pair(std::make_pair(BarrierEnergy(sequence, A, B, pathway), pathway), strategy_id));
			};

			for (int i = 0; i < generation_plan[0]; ++i) {
				RegistrataBirth(M1(A, B, x.first.second, rnd), 0);
			}
			for (int i = 0; i < generation_plan[1]; ++i) {
				RegistrataBirth(M2(A, B, x.first.second, rnd), 1);
			}
			for (int i = 0; i < generation_plan[2]; i += 2) {
				RegistrataBirth(M3Plus(sequence, A, B, x.first.second, rnd), 2);
				RegistrataBirth(M3Minus(sequence, A, B, x.first.second, rnd), 2);
			}
			for (int i = 0; i < generation_plan[3]; ++i) {
				RegistrataBirth(M4(sequence, A, B, x.first.second, rnd), 3);
			}
			for (int i = 0; i < generation_plan[4]; ++i) {
				RegistrataBirth(M5(sequence, A, B, x.first.second, rnd), 4);
			}

			std::sort(kindergarten.begin(), kindergarten.end());
			for (int i = 0; i < l2_survive_population && i < kindergarten.size(); ++i) {
				next_population.push_back(kindergarten[i]);
			}
		}

		//Fig. 2の12行目
		std::sort(next_population.begin(), next_population.end());
		best_pathways_history.push_back(next_population[0].first);
		std::cout << "LOG: RNAEAPath: k = " << k << ", best barrier = "<< best_pathways_history.back().first << std::endl;

		//Fig. 2の13行目
		population = next_population;
		if (population.size() > l3_final_max_population) {
			population.resize(l3_final_max_population);
		}

		//論文の
		//>The number of offsprings produced by each mutation strategy
		//の部分の記述に従い、survive_from_the_strategyとgeneration_planを更新する。
		//論文では過去の履歴で今の値を決めるという書き方だが、
		//ここでは現在の履歴で次の値を決めるということにする。結果は等価である。

		std::vector<int>survive_from_the_strategy(5, 0);
		for (const auto p : population) {
			if (0 <= p.second && p.second < 5)++survive_from_the_strategy[p.second];
		}
		//この時点で、survive_from_the_strategyは論文のb^{k}_{y'}に等しい。

		std::vector<double>b_l_ratio(5, 0);
		for (int i = 0; i < 5; ++i)b_l_ratio[i] = double(survive_from_the_strategy[i]) / double(generation_plan[i]);
		//この時点で、b_l_ratioは論文の(b^{k}_{y'} / l^{k}_{M_{y'}})に等しい。

		//論文の式の値を計算する。l^{k+1}_{M{y}}=generation_plan[i]
		double b_l_ratio_sum = 0.0;
		for (int i = 0; i < 5; ++i)b_l_ratio_sum += b_l_ratio[i];
		for (int i = 0; i < 5; ++i)generation_plan[i] = std::max<int>(3, (b_l_ratio[i] / b_l_ratio_sum) * double(script_capital_l));
	}


	std::vector<std::string>pathway{ structure1 };
	auto S = A;
	for (const auto action : best_pathways_history.back().second) {
		if (S.find(action) != S.end()) {
			S.erase(action);
		}
		else {
			for (const auto x : S)assert(!IsExclusive(x, action));
			S.insert(action);
		}
		pathway.push_back(BasePairSetToDotNotation(S, sequence.size()));
	}
	assert(S == B);
	assert(pathway.back() == structure2);

	pathway = TrimPathway(pathway);
	return std::make_pair(pathway, BarrierEnergy(sequence, pathway));
}



}