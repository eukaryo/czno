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
	//action_list[t1]���������āA�����action_list[t2]�̒���ɑ}������B

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2 && t2 < action_list.size());
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	std::vector<std::pair<int, int>>new_action_chain;
	if (t1 < t2) {
		//[t1]�𔲂���[t2]�̒���ɓ����B
		for (int i = 0; i < t1; ++i)new_action_chain.push_back(action_list[i]);
		for (int i = t1+1; i <= t2; ++i)new_action_chain.push_back(action_list[i]);
		new_action_chain.push_back(action_list[t1]);
		for (int i = t2 + 1; i < action_list.size(); ++i)new_action_chain.push_back(action_list[i]);
	}
	else if(t2 + 1 < t1) {
		//[t1]�𔲂���[t2]�̒���ɓ����B
		for (int i = 0; i <= t2; ++i)new_action_chain.push_back(action_list[i]);
		new_action_chain.push_back(action_list[t1]);
		for (int i = t2 + 1; i < t1; ++i)new_action_chain.push_back(action_list[i]);
		for (int i = t1 + 1; i < action_list.size(); ++i)new_action_chain.push_back(action_list[i]);
	}
	else new_action_chain = action_list;//�������ꏊ�ɓ����A�܂艽���ς��Ȃ�
	return new_action_chain;
}

std::vector<std::pair<int, int>>M2Mechanism(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	std::vector<std::pair<int, int>> action_list,
	const int t1,
	const int t2) {
	//[t1]��[t2]�����ւ���B

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
	//action_list[t1]�̒����action_list[t2]�̒����x��}������B

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
	//descending: 0���o�₷���Asize-1�͏o�ɂ����B
	//ascending: size-1���o�₷���A0�͏o�ɂ����B
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
	//action_list[t1]���������A�����action_list[t2lb,t2ub]�͈͓̔��̂ǂꂩ�̒���ɑ}������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�̂����A�L�]�����ȑ}���ʒu���I�΂�₷���悤��bias��������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�����݂��Ȃ��Ȃ�΋���ۂ�pathway��Ԃ��B

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2lb && t2lb <= t2ub && t2ub < action_list.size());
	assert(t2ub <= t1 || t1 <= t2lb);//t2�͈̔͂�t1�̑O����낾��
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	std::vector<std::vector<std::pair<int, int>>>valid_mutations;

	if (t2ub <= t1) {
		//t1���O�ɂ���t2�̒����t1��}������Ƃ����ψكp�X�̂����Avalid�ȃp�X��S�񋓂���B
		for (int t2 = std::max<int>(0, t2lb); t2 <= t1 && t2 <= t2ub; ++t2) {
			const auto candidate = M1Mechanism(A, B, action_list, t1, t2);
			if (IsValidPathway(A, B, candidate))valid_mutations.push_back(candidate);
		}
	}
	else {
		//t1�����ɂ���t2�̒����t1��}������Ƃ����ψكp�X�̂����Avalid�ȃp�X��S�񋓂���B
		for (int t2 = std::min<int>(t1, t2lb); t2 < action_list.size() && t2 <= t2ub; ++t2) {
			const auto candidate = M1Mechanism(A, B, action_list, t1, t2);
			if (IsValidPathway(A, B, candidate))valid_mutations.push_back(candidate);
		}
	}

	//valid�ȃp�X�������ꍇ�͏I��
	if (valid_mutations.size() == 0)return std::vector<std::pair<int, int>>{};

	//t1������΂�g�ޑ���Ȃ̂��O������Ȃ̂��𒲂ׂ�B
	int count = A.find(action_list[t1]) != A.end() ? 1 : 0;
	for (int i = 0; i < t1; ++i)if (action_list[i] == action_list[t1])count++;

	//���̎��_��count����Ȃ�t1�͋�����ڂł���O�������ł���B
	//�䂦�ɂ��̂Ƃ���t2��t1�Ƌ߂��ق����悢�Bt2�̂ق����x���Ȃ�DescendingBias���~�����B
	//t2���O�̏ꍇ(M5�Ŏg��)�Abias�̌����͘_���ɏ�����Ă��Ȃ��B����Ĉ��S��Ƃ���t2��t1���߂��ق����悢���Ƃɂ���B
	//���Ȃ킿(t1 <= t2lb)��false�Ȃ�����͏��false�Ƃ���B
	return valid_mutations[GaussianBiasedSelection(valid_mutations.size(), rnd, (count % 2 == 1) & (t1 <= t2lb))];

}

std::vector<std::pair<int, int>>M1(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//�_���̕ψِ�����@M1�����B���s���������ۂ�pathway��Ԃ��B

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
	//�_���̕ψِ�����@M2�����B���s���������ۂ�pathway��Ԃ��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//�X���b�v������action�̑g(i,j)��S�񋓂��ăV���b�t������B
	std::vector<std::pair<int, int>>candidates;
	for (int i = 0; i < action_list.size(); ++i) {
		for (int j = i + 1; j < action_list.size(); ++j) {
			candidates.push_back(std::make_pair(i, j));
		}
	}
	std::sort(candidates.begin(), candidates.end());//���ʂ̍Č�����ۂ��߂Ƀ\�[�g����B
	std::shuffle(candidates.begin(), candidates.end(), rnd);

	for (const auto x : candidates) {
		const auto candidate = M2Mechanism(A, B, action_list, x.first, x.second);
		if (IsValidPathway(A, B, candidate))return candidate;
	}

	//�ǂ��V���b�t�����Ă�invalid�ɂȂ�Ȃ�I���B
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
	//action_list[t1]�̒����x��}�����A����
	//action_list[t2lb,t2ub]�͈͓̔��̂ǂꂩ�̒���ɂ�x��}������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�̂����A�L�]�����ȑ}���ʒu���I�΂�₷���悤��bias��������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�����݂��Ȃ��Ȃ�΋���ۂ�pathway��Ԃ��B

	assert(0 <= t1 && t1 < action_list.size());
	assert(0 <= t2lb && t2lb <= t2ub && t2ub < action_list.size());
	assert(t2ub <= t1 || t1 <= t2lb);//t2�͈̔͂�t1�̑O����낾��
	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//action_list[t1]�̒����x��}�����邱�Ǝ��̂�invalid�ȏꍇ�A
	//IsValidPathway�����ׂ�false�ɂȂ邱�ƂŁA���ǋ���ۂ�pathway���Ԃ����B
	//����āA����������Ō������Ȃ��Ă��悢�B

	std::vector<std::vector<std::pair<int, int>>>valid_mutations;

	//3. �����_���ɑI�񂾉����x�ɑ΂��āA�����t1�̒���ɑg�񂾂Ƃ��āA
	//������O���^�C�~���O�Ƃ���valid�ȃ^�C�~���Ot2��S�񋓂���B
	if (t2ub <= t1) {
		for (int t2 = std::max<int>(0, t2lb); t2 + 1 < t1 && t2 <= t2ub; ++t2) {
			const auto new_action_list = M3Mechanism(A, B, action_list, x, t1, t2);
			if (IsValidPathway(A, B, new_action_list))valid_mutations.push_back(new_action_list);
		}
	}
	else {
		//t1�����ɂ���t2�̒����t1��}������Ƃ����ψكp�X�̂����Avalid�ȃp�X��S�񋓂���B
		for (int t2 = std::min<int>(t1 + 1, t2lb); t2 < action_list.size() && t2 <= t2ub; ++t2) {
			const auto new_action_list = M3Mechanism(A, B, action_list, x, t1, t2);
			if (IsValidPathway(A, B, new_action_list))valid_mutations.push_back(new_action_list);
		}
	}

	//�_���ɂ͏�����Ă��Ȃ����A���ׂĂ�t2��invalid��pathway���������Ȃ����Ƃ��l������B
	if (valid_mutations.size() == 0)return std::vector<std::pair<int, int>>{};
	
	//plus, t1<t2 �� false //t1�őg�񂾂�ł��邾������ĊO�������{�Y�����傫���قǗ����
	//plus, t2<t1 �� true  //t2�őg�񂾂�ł��邾�������ĊO�������{�Y�����������قǗ����
	//minus,t1<t2 �� true  //t1�ŗ�������ł��邾���߂��őg�݂����{�Y�����傫���قǗ����
	//minus,t2<t1 �� false //t2�ŗ�������ł��邾���߂��őg�݂����{�Y�����������قǗ����
	return valid_mutations[GaussianBiasedSelection(valid_mutations.size(), rnd, (t2ub <= t1) ^ is_plus)];
}

std::vector<std::pair<int, int>>M3Plus(
	const std::string& sequence, 
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//�_���̕ψِ�����@M3+�����B���s���������ۂ�pathway��Ԃ��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_list�����l�ɑI�сAt1�Ƃ���B
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	const int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1�̒����structure S������ΏW���Ƃ��ċ��߂�B������action_list��valid�ł��邱�Ƃ����肷��B
	//�����F
	//t1�̒���ɉ���΂�}������Ƃ�������͘_���̋L�q�ʂ�ł���B
	//>introducing an addition action add i,j after at1
	//�_���I�ɂ́A�X�^�[�g�\���̒���ɉ���΂�}�����Ă����Ȃ��͂������A
	//�_���̋L�q�ʂ�Ɏ��������̂ŃX�^�[�g�\���̒���ɉ���΂��}������邱�Ƃ͂Ȃ��B
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. ���݂�structure S�ɑ΂��āAvalid����compatible�ɑg�߂鉖���(i,j)��S�񋓂���B
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

	//valid����compatible�ɑg�߂鉖��ΑS�Ă��烉���_����1�I�ԁB
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
	//�_���̕ψِ�����@M3-�����B���s���������ۂ�pathway��Ԃ��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_list�����l�ɑI�сAt1�Ƃ���B
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	const int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1�̒����structure S������ΏW���Ƃ��ċ��߂�B������action_list��valid�ł��邱�Ƃ����肷��B
	//
	//t1�̒���ɉ���΂�}������Ƃ�������͘_���̋L�q�ʂ�ł���B
	//>introducing an addition action add i,j after at1
	//�_���I�ɂ́A�X�^�[�g�\���̒���ɉ���΂�}�����Ă����Ȃ��͂������A
	//�_���̋L�q�ʂ�Ɏ��������̂ŃX�^�[�g�\���̒���ɉ���΂��}������邱�Ƃ͂Ȃ��B
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. ���݂�structure S�ɑ΂��āAvalid�ɊO���鉖���(i,j)��S�񋓂��āA���̒����烉���_���ɑI�ԁB
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
	//action_list[t1]���������x���o������Ȃ�AM1�̃X�g���e�W�[�ł����t1�̒���Ɉړ�����B
	//action_list[t1]���������x���o�����Ȃ��Ȃ�Aaction_list[t1]�̒����x��}�����A����
	//action_list[t2lb,t2ub]�͈͓̔��̂ǂꂩ�̒���ɂ�x��}������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�̂����A�L�]�����ȑ}���ʒu���I�΂�₷���悤��bias��������B
	//action_list��valid�ɂȂ�悤�ȑ}���ʒu�����݂��Ȃ��Ȃ�΋���ۂ�pathway��Ԃ��B
	//action_list[t1]�̒����x�����݂���Ƃ��A����x�͉���Ό`���ł��菜���ł͂Ȃ��A�����ꎩ�̂�valid�ł���Ƃ������Ƃ����肷��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));
	
	//action_list��[t1]������ɉ����(i,j)���o�Ă��邩�ǂ������ׂ�B
	int same_action_pos = -1;
	for (int i = t1 + 1; i < action_list.size(); ++i) {
		if (action_list[i] == x) {
			same_action_pos = i;
			break;
		}
	}

	if (same_action_pos >= 0) {
		//3.1 action_list��[t1]������ɉ����(i,j)���o�Ă���Ȃ�A
		//�����"�����グ��"t1�̒���Ɏ����Ă���BM1�̃X�g���e�W�[�Ŏ��s����B
		//>move it up and place it after a_t using strategy M1.
		return M1Strategy(A, B, action_list, same_action_pos, t1, t1, rnd);
	}
	//3.2 action_list��[t1]������ɉ����(i,j)���o�Ă��Ȃ��Ȃ�A
	//[t1]�̒����(i,j)��}�����A������Ɍ����(i,j)���������鏈����}������B
	//�����M3�̃X�g���e�W�[�Ŏ��s����B
	return M3Strategy(sequence, A, B, action_list, t1, t1, action_list.size() - 1, x, true, rnd);
}

std::vector<std::pair<int, int>>M4(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list,
	std::mt19937_64& rnd) {
	//�_���̕ψِ�����@M4�����B���s���������ۂ�pathway��Ԃ��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	const int n = sequence.length();

	//1. action_list�����l�ɑI�сAt1�Ƃ���B
	if (action_list.size() == 0)return std::vector<std::pair<int, int>>{};
	int t1 = std::uniform_int_distribution<int>(0, action_list.size() - 1)(rnd);

	//t1�̒����structure S������ΏW���Ƃ��ċ��߂�B������action_list��valid�ł��邱�Ƃ����肷��B
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. ����4�ȏ�̘A�����������(stack)�̂����A
	//�񎟍\��S�ɑ΂��đ����ɒǉ��ł���(compatible)�悤�Ȃ��̂�S�񋓂��A�����_����1�I�ԁB
	//���Ȃ݂ɒ���4�͘_���̒l�B
	//>all possible stacks with more than 3 consecutive base pairs
	auto stacks = EnumerateLongStacks(sequence, S, 4, true);
	std::sort(stacks.begin(), stacks.end());//�Č����̂��߃\�[�g����B
	if (stacks.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto target_stack = stacks[std::uniform_int_distribution<int>(0, stacks.size() - 1)(rnd)];

	//3. 2.�őI�񂾃X�e�����ARNA�̓������珇�Ԃɓ���Ă����B
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
	//�_���̕ψِ�����@M5�����B���s���������ۂ�pathway��Ԃ��B

	/*HEAVY*/assert(IsValidPathway(A, B, action_list));

	//1. action_list�̂�������Ώ���������̂̓Y�����������l�ɑI�сAt1�Ƃ���B
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

	//t1�̒���̓񎟍\��������ΏW���Ƃ��ċ��߂�B������action_list��valid�ł��邱�Ƃ����肷��B
	const auto S = AssosiatedIntermediateStructure(A, action_list, t1);

	//2. ����4�ȏ�̘A�����������(stack)�̂����A
	//�񎟍\��S�ɑ΂��đ����ɒǉ��ł��Ȃ�(incompatible)�悤�Ȃ��̂�S�񋓂��A�����_����1�I�ԁB
	//���Ȃ݂ɒ���4�͘_���̒l�B
	//>all possible stacks with more than 3 consecutive base pairs
	auto stacks = EnumerateLongStacks(sequence, S, 4, false);
	std::sort(stacks.begin(), stacks.end());//�Č����̂��߃\�[�g����B
	if (stacks.size() == 0)return std::vector<std::pair<int, int>>{};
	const auto target_stack = stacks[std::uniform_int_distribution<int>(0, stacks.size() - 1)(rnd)];

	//target_stack���̉���΂������Acompatible���ǂ����ŕ��ނ���B
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

	//target_stack���̉���΂������ARNA�̓������珇�Ԃɓ���Ă����B
	//compatible�Ȃ��S������Ă���Aincompatible�Ȃ�����Ă����B
	std::vector<std::pair<int, int>>new_action_list = action_list;
	for (const auto x : compatible_pairs) {

		//3.2 compatible�Ȃ��̂�M4�̃X�g���e�W�[�œ����B
		new_action_list = M4Strategy(sequence, A, B, new_action_list, t1, x, rnd);
		if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
		assert(new_action_list[t1 + 1] == x);
		++t1;
	}
	for (const auto x : incompatible_pairs) {

		//4.1 x�ɑ΂���exclusive��S��̉���΂��ꂼ��ɂ��đΏ�����B
		for (const auto y : S)if (IsExclusive(x, y)) {

			//y��action_list[t1]��������ɏo�Ă��邩�ŏꍇ��������B
			int same_action_pos = -1;
			for (int i = t1 + 1; i < new_action_list.size(); ++i) {
				if (new_action_list[i] == y) {
					same_action_pos = i;
					break;
				}
			}

			if (same_action_pos >= 0) {
				//4.2. y��action_list[t1]��������ɏo�Ă���Ȃ�A�����"�����グ��"�At1���O���ɑ}������B
				//M1�̃X�g���e�W�[���g���B
				new_action_list = M1Strategy(A, B, new_action_list, same_action_pos, 0, t1, rnd);
			}
			else {
				//4.3. y��action_list[t1]��������ɏo�Ă��Ȃ��Ȃ�Aaction_list[t1]�̒����y��}�����A
				//����Ɍ���ɂ�y��}������BM3�̃X�g���e�W�[���g���B
				new_action_list = M3Strategy(sequence, A, B, new_action_list, t1, t1, new_action_list.size() - 1, y, false, rnd);
			}
			if (new_action_list.size() == 0)return std::vector<std::pair<int, int>>{};
			++t1;
		}
		//���̎��_�ŁAaction_state[t1]�̒���̒��ԍ\���ɑ΂���x��compatible�ł���B
		//����āA3.2�Ɠ�����M4�̃X�g���e�W�[��x��������B
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
	const int initial_population,//�ŏ��ɓK���ɍ��p�X�̐��B4
	const int gamma_generation,//�I�������Ɋւ��� 5
	const int max_generation,//�I�������Ɋւ��MAX 10
	const int script_capital_l,//�q���̐��Ɋւ��A100
	const int l1_elite_population,//�G���[�g�̐��A5
	const int l2_survive_population,// 5
	const int l3_final_max_population,// 100
	const int random_seed) {

	ValidateTriple(sequence, structure1, structure2);

	const int n = sequence.length();
	std::mt19937_64 rnd(random_seed);

	//Fig. 2��1�s�ځB�_���̋L�@����max�ł͂Ȃ�abs�����Amax�ł������ȏ����͂ł���B
	const double delta_energy = std::max<double>(
		EnergyOfStructure(sequence, structure1),
		EnergyOfStructure(sequence, structure2));

	//simple pathway����邽�߂̏���
	const auto A = DotNotationToBasePairSet(structure1);
	const auto B = DotNotationToBasePairSet(structure2);
	std::vector<std::pair<int, int>>need_to_close;
	std::vector<std::pair<int, int>>need_to_open;
	for (const auto x : B)if (A.find(x) == A.end())need_to_close.push_back(x);
	for (const auto x : A)if (B.find(x) == B.end())need_to_open.push_back(x);
	std::sort(need_to_close.begin(), need_to_close.end());
	std::sort(need_to_open.begin(), need_to_open.end());

	//Fig. 2��3�s�ځB��ԍŌ��int�́A���ꂪ�ǂ̃X�g���e�W�[�ŗR����
	std::vector<std::pair<std::pair<double, std::vector<std::pair<int, int>>>, int>>population;
	for (int i = 0; i < initial_population; ++i) {
		const auto pathway = GenerateRandomSimplePathway(A, B, rnd);
		population.push_back(std::make_pair(std::make_pair(BarrierEnergy(sequence, A, B, pathway), pathway), -1));
	}
	std::sort(population.begin(), population.end());

	//Fig. 2��4�s��
	std::vector<std::pair<double, std::vector<std::pair<int, int>>>>best_pathways_history{ population[0].first };

	//k��Fig. 2��2�s�ڂɏo�Ă����B �_���ɂ���Stopping conditions�̔�����s���B
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

	//Fig. 2��5,6�s��
	for (int k = 0; !IsStoppingCondition(k); ++k) {
		std::cout << "LOG: RNAEAPath: start: k = " << k << ", pop = " << population.size() << std::endl;

		//Fig. 2��7�s�ځB�G���[�g������B
		auto next_population = population;
		next_population.resize(std::min<int>(next_population.size(), l1_elite_population));

		//Fig. 2��8�s��
		for (const auto x : population) {

			//Fig. 9�s�ڂ̍��ӂ�T���`����B
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

		//Fig. 2��12�s��
		std::sort(next_population.begin(), next_population.end());
		best_pathways_history.push_back(next_population[0].first);
		std::cout << "LOG: RNAEAPath: k = " << k << ", best barrier = "<< best_pathways_history.back().first << std::endl;

		//Fig. 2��13�s��
		population = next_population;
		if (population.size() > l3_final_max_population) {
			population.resize(l3_final_max_population);
		}

		//�_����
		//>The number of offsprings produced by each mutation strategy
		//�̕����̋L�q�ɏ]���Asurvive_from_the_strategy��generation_plan���X�V����B
		//�_���ł͉ߋ��̗����ō��̒l�����߂�Ƃ��������������A
		//�����ł͌��݂̗����Ŏ��̒l�����߂�Ƃ������Ƃɂ���B���ʂ͓����ł���B

		std::vector<int>survive_from_the_strategy(5, 0);
		for (const auto p : population) {
			if (0 <= p.second && p.second < 5)++survive_from_the_strategy[p.second];
		}
		//���̎��_�ŁAsurvive_from_the_strategy�͘_����b^{k}_{y'}�ɓ������B

		std::vector<double>b_l_ratio(5, 0);
		for (int i = 0; i < 5; ++i)b_l_ratio[i] = double(survive_from_the_strategy[i]) / double(generation_plan[i]);
		//���̎��_�ŁAb_l_ratio�͘_����(b^{k}_{y'} / l^{k}_{M_{y'}})�ɓ������B

		//�_���̎��̒l���v�Z����Bl^{k+1}_{M{y}}=generation_plan[i]
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