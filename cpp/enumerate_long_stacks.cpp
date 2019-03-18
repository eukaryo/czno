/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
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
	//sequence��ł̒���minimum_length�ȏ�̘A�����鉖���(stack)���A
	//�񎟍\��structure�ɑ����ɒǉ��ł���(compatible)���ǂ����ŕ��ނ��āA�Е���Ԃ��B
	//need_compatible==true�Ȃ�Acompatible��stack��Ԃ��B
	//need_compatible==false�Ȃ�Aincompatible��stack��Ԃ��B
	//brute force�Ɍv�Z����B
	//�estack���ɉ���΂�O(N)����Astructure���ɉ���΂�O(N)����̂ŁA
	//���stack��compatible���ǂ����̔���ɂ�O(N^2)������B
	//stack��O(N^3)���邩��A�g�[�^���̎��Ԍv�Z�ʂ�O(N^5)�ł���B

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

	//stack��S�񋓂���B
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

	//�estack��compatibility�𔻒肵�āA�~���������W�߂�B
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
		//�����(i,j)�����݂���Ƃ��A
		//[(0,0),(i,j)]�̋�`��Ԃ�[(i,j),(n-1,n-1)]�̋�`��Ԃ�1�ȏ�̒l�����Z����B
		//�������@�Œ[��4�ӏ��������Z����B
		++incompatible_flag[0][0];
		--incompatible_flag[0][bp.second + 1];
		--incompatible_flag[bp.first + 1][0];
		++incompatible_flag[bp.first + 1][bp.second + 1];
		++incompatible_flag[bp.first][bp.second];
		--incompatible_flag[bp.first][n];
		--incompatible_flag[n][bp.second];
		++incompatible_flag[n][n];
	}

	//�������@�̗ݐϘa����鏈��
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

	//���̎��_�ŁA
	//incompatible_flag[i][j] > 0 �� �����(i,j)��structure�Ƃ�incompatible
	//�ł���B

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

