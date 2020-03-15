/*
GNU GPL v2
Copyright (c) 2020 Hiroki Takizawa
*/

#include "enumerate_long_stacks.h"

#include "misc.h"

namespace czno_cpp {

std::vector<std::array<int,4>>
EnumerateLongStacksBruteForce(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& structure,
	const int minimum_length,
	const bool need_compatible) {
	//sequence上での長さminimum_length以上の連続する塩基対(stack)を、
	//二次構造structureに即座に追加できる(compatible)かどうかで分類して、片方を返す。
	//need_compatible==trueなら、compatibleなstackを返す。
	//need_compatible==falseなら、incompatibleなstackを返す。
	//brute forceに計算する。
	//各stack中に塩基対はO(N)個あり、structure中に塩基対はO(N)個あるので、
	//一つのstackがcompatibleかどうかの判定にはO(N^2)かかる。
	//stackはO(N^3)個あるから、トータルの時間計算量はO(N^5)である。

	const static std::regex rna(R"([ACGU]+)");
	assert(std::regex_match(sequence, rna));
	const int n = sequence.length();
	for (const auto bp : structure)assert(0 <= bp.first && bp.first < bp.second && bp.second < n);

	std::vector<std::vector<int>>is_valid_bp(n + 1, std::vector<int>(n + 1, 0));
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			is_valid_bp[i][j] = IsValidBasePair(sequence, std::make_pair(i, j)) ? 1 : 0;
		}
	}

	//stackを全列挙する。
	std::vector<std::array<int, 4>>all_long_stacks;
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			if (is_valid_bp[i][j] == 0)continue;
			for (int x = 1; i - x >= 0 && j + x < n; ++x) {
				if (is_valid_bp[i - x][j + x] == 0)break;
				if (x + 1 >= minimum_length)all_long_stacks.push_back({ i, j, i - x, j + x });
			}
		}
	}

	//各stackのcompatibilityを判定して、欲しい方を集める。
	std::vector<std::array<int, 4>>answer;
	for (const auto s : all_long_stacks) {
		bool compatibility = true;
		for (int i = s[0], j = s[1]; i >= s[2] && compatibility; --i, ++j) {
			for (const auto bp : structure) {
				if (IsExclusive(std::make_pair(i, j), bp)) {
					compatibility = false;
					break;
				}
			}
		}
		if (compatibility == need_compatible)answer.push_back(s);
	}

	const auto equivalent_answer = EnumerateLongStacks(sequence, structure, minimum_length, need_compatible);
	assert(equivalent_answer == answer);

	return answer;
}

std::vector<std::array<int, 4>>
EnumerateLongStacks(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& structure,
	const int minimum_length,
	const bool need_compatible) {

	const static std::regex rna(R"([ACGU]+)");
	assert(std::regex_match(sequence, rna));
	const int n = sequence.length();
	for (const auto bp : structure)assert(0 <= bp.first && bp.first < bp.second && bp.second < n);

	std::vector<std::vector<int>>incompatible_flag(n + 1, std::vector<int>(n + 1, 0));
	for (const auto bp : structure) {
		//塩基対(i,j)が存在するとき、
		//[(0,0),(i,j)]の矩形区間と[(i,j),(n-1,n-1)]の矩形区間に1以上の値を加算する。
		//いもす法で端の4箇所だけ加算する。
		++incompatible_flag[0][0];
		--incompatible_flag[0][bp.second + 1];
		--incompatible_flag[bp.first + 1][0];
		++incompatible_flag[bp.first + 1][bp.second + 1];
		++incompatible_flag[bp.first][bp.second];
		--incompatible_flag[bp.first][n];
		--incompatible_flag[n][bp.second];
		++incompatible_flag[n][n];
	}

	//いもす法の累積和を取る処理
	for (int i = 0; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			incompatible_flag[i][j] += incompatible_flag[i][j - 1];
		}
	}
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			incompatible_flag[i][j] += incompatible_flag[i - 1][j];
		}
	}

	//この時点で、
	//incompatible_flag[i][j] > 0 ⇔ 塩基対(i,j)とstructureとはincompatible
	//である。

	std::vector<std::vector<int>>is_valid_bp(n + 1, std::vector<int>(n + 1, 0));
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			is_valid_bp[i][j] = IsValidBasePair(sequence, std::make_pair(i, j)) ? 1 : 0;
		}
	}

	std::vector<std::array<int, 4>>answer;
	for (int i = 0; i < n - TURN - 1; ++i) {
		for (int j = i + TURN + 1; j < n; ++j) {
			if (is_valid_bp[i][j] == 0)continue;
			bool compatibility = incompatible_flag[i][j] == 0;
			for (int x = 1; i - x >= 0 && j + x < n; ++x) {
				if (is_valid_bp[i - x][j + x] == 0)break;
				if (incompatible_flag[i - x][j + x] > 0)compatibility = false;
				if (need_compatible) {
					if (!compatibility)break;
					if (x + 1 >= minimum_length)answer.push_back({ i, j, i - x, j + x });
				}
				else {
					if(!compatibility && x + 1 >= minimum_length)answer.push_back({ i, j, i - x, j + x });
				}
			}
		}
	}

	return answer;
}

}

