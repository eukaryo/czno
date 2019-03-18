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

	//�ϐ������̏���
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

	//Fitness�֐���[Dotu et al.,2010]��Fig10�̉��ɂ��鎮
	const auto Fitness = [&](const std::string& x) {
		return EnergyOfStructure(sequence, x) + w * HammingDistance(x, structure2);
	};

	for (int count = 0; count < 1000 && S != structure2; ++count) {

		/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));
		std::cout << "LOG: RNATABUPATH: start: count = " << count << std::endl;
		//S����n�~���O����1�Ⴄ�\�������̂����A�G�l���M�[�����l�ȉ��̂��̂�
		//�S�ăt�B�b�g�l�X�֐��̒l�ƍ��킹�Ċi�[����B
		std::vector<std::pair<double, std::string>>low_energy_neighbors;
		for (const auto x : EnumerateNeighborsWithoutTabu(sequence, S, count, tabu)) {
			if (EnergyOfStructure(sequence, x) <= energy_upper_bound) {
				low_energy_neighbors.push_back(std::make_pair(Fitness(x), x));
			}
		}

		//�J�ډ\�ȍ\�������݂��Ȃ��ꍇ�Atabu����΂���������܂ő҂B
		if (low_energy_neighbors.size() == 0)continue;

		//�t�B�b�g�l�X�֐��̕]���l���Ⴂk�̂Ȃ����烉���_���ɑI�ёJ�ڂ���B
		std::sort(low_energy_neighbors.begin(), low_energy_neighbors.end());
		std::uniform_int_distribution<int>dist(0, std::min<int>(k, low_energy_neighbors.size()) - 1);
		const auto nextS = low_energy_neighbors[dist(rnd)].second;
		const auto diff = Difference(S, nextS);
		assert(tabu[diff.first][diff.second] <= count);
		S = nextS;
		pathway.push_back(S);
		tabu[diff.first][diff.second] = count + tabu_tenure(rnd);//��s�����̎����ɍ��킹���B
		max_energy = std::max<double>(max_energy, EnergyOfStructure(sequence, S));
		/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));

		//���܂łň�ԃS�[���ɋ߂��Ȃ������ǂ����ŏꍇ����
		const int hamming_distance = HammingDistance(S, structure2);
		if (hamming_distance < best_distance) {
			best_distance = hamming_distance;
			closest_structure = S;
			closest_time = pathway.size();
			w = std::max<double>(w_0, w / 1.2);//��s�����̎����ɍ��킹���B
			no_improvement_count = 0;
		}
		else {
			no_improvement_count++;
			if (no_improvement_count > max_stable) {

				//���񐔒�؂�����A��ԃS�[���ɋ߂������\���Ɋ����߂��B

				S = closest_structure;
				pathway.resize(closest_time);
				/*HEAVY*/assert(IsValidPathway(structure1, S, pathway));
				for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)tabu[i][j] = 0;
				w *= 2;//��s�����̎����ɍ��킹���B
				no_improvement_count = 0;

				if (w > 10.0 * w_init) {

					//��؂��Ċ����߂������P���Ȃ��Ƃ����̂����񐔌J��Ԃ�����A
					//�ŏ��̍\���Ɋ����߂��B

					S = structure1;
					pathway.resize(1);
					for (int i = 0; i < n; ++i)for (int j = 0; j < n; ++j)tabu[i][j] = 0;
					w_init += 0.5;//��s�����̎����ɍ��킹���B
					w = w_init;

					//�ȉ���3�s�͘_���ɂ͏�����Ă��Ȃ���Dotu et al.,�̎����ɂ͂���B�K�{�B

					closest_structure = structure1;
					closest_time = 1;
					best_distance = HammingDistance(structure1, structure2);
				}
			}
		}
	}

	if (S != structure2) {
		//energy_upper_bound���Ⴗ�������A
		//w_0�����������Ď�������O��count���I�[�o�[�����B
		return make_pair(std::vector<std::string>{}, std::numeric_limits<double>::infinity());
	}

	pathway = TrimPathway(pathway);
	return make_pair(pathway, BarrierEnergy(sequence, pathway));
}

}