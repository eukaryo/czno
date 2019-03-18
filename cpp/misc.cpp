/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "misc.h"

#include "call_vienna_rna.h"

namespace czno_cpp {


std::vector<std::string> Split(const std::string& str, const char delimiter) {
	std::vector<std::string> ans;
	std::stringstream s(str);
	std::string tmp;
	while (std::getline(s, tmp, delimiter))ans.push_back(tmp);
	return ans;
}

bool IsValidBasePair(const std::string& sequence, const std::pair<int, int>& x) {
	const static std::string bp("AU CG GC GU UA UG");
	std::string y = std::string("..");
	y[0] = sequence[x.first];
	y[1] = sequence[x.second];
	return bp.find(y) != std::string::npos;
}
std::set<std::pair<int, int>> DotNotationToBasePairSet(const std::string& structure) {
	std::set<std::pair<int, int>>answer;
	std::stack<int>s;
	for (int i = 0; i < structure.size(); ++i) {
		if (structure[i] == '(')s.push(i);
		else if (structure[i] == ')') {
			const int x = s.top();
			s.pop();
			answer.insert(std::make_pair(x, i));
		}
		else assert(structure[i] == '.');
	}
	return answer;
}
std::string BasePairSetToDotNotation(const std::set<std::pair<int, int>>& structure, const int len) {
	std::string answer("");
	for (int i = 0; i < len; ++i)answer += ".";
	for (const auto x : structure) {
		assert(x.first < x.second && answer[x.first] == '.' && answer[x.second] == '.');
		answer[x.first] = '(';
		answer[x.second] = ')';
	}
	return answer;
}
int HammingDistance(const std::string& structure1, const std::string& structure2) {

	assert(structure1.size() == structure2.size());

	const auto A = DotNotationToBasePairSet(structure1);
	const auto B = DotNotationToBasePairSet(structure2);

	std::set<std::pair<int, int>>intersection;
	std::set_intersection(
		A.begin(), A.end(),
		B.begin(), B.end(),
		inserter(intersection, intersection.end()));

	std::set<std::pair<int, int>>symmetric_difference;
	std::set_symmetric_difference(
		A.begin(), A.end(),
		B.begin(), B.end(),
		inserter(symmetric_difference, symmetric_difference.end()));

	assert(int(A.size()) + int(B.size()) - 2 * int(intersection.size()) == int(symmetric_difference.size()));
	return int(A.size()) + int(B.size()) - 2 * int(intersection.size());
}
bool IsExclusive(const std::pair<int, int>& base_pair1, const std::pair<int, int>& base_pair2) {
	assert(base_pair1.first < base_pair1.second);
	assert(base_pair2.first < base_pair2.second);
	if (base_pair1.first == base_pair2.first)return true;
	if (base_pair1.first == base_pair2.second)return true;
	if (base_pair1.second == base_pair2.first)return true;
	if (base_pair1.second == base_pair2.second)return true;
	if (base_pair1.first < base_pair2.first &&  base_pair2.first < base_pair1.second && base_pair1.second < base_pair2.second)return true;
	if (base_pair2.first < base_pair1.first && base_pair1.first < base_pair2.second && base_pair2.second < base_pair1.second)return true;
	return false;
}
void ValidateTriple(const std::string& sequence, const std::string& structure1, const std::string& structure2) {
	const int n = sequence.length();
	assert(structure1.length() == n);
	assert(structure2.length() == n);
	const static std::regex rna(R"([ACGU]+)");
	const static std::regex ss(R"([\(\.\)]+)");
	const static std::string bp("AU CG GC GU UA UG");
	assert(std::regex_match(sequence, rna));
	assert(std::regex_match(structure1, ss));
	assert(std::regex_match(structure2, ss));
	for (const auto x : DotNotationToBasePairSet(structure1)) {
		assert(IsValidBasePair(sequence, x));
	}
	for (const auto x : DotNotationToBasePairSet(structure2)) {
		assert(IsValidBasePair(sequence, x));
	}
}
std::vector<std::string>TrimPathway(const std::vector<std::string>& pathway) {
	std::vector<std::string>answer;
	std::map<std::string, int>footprint;
	for (int i = 0; i < pathway.size(); ++i)footprint[pathway[i]] = i;
	for (int i = 0; i < pathway.size(); ++i) {
		if (i != footprint[pathway[i]])i = footprint[pathway[i]];
		answer.push_back(pathway[i]);
	}
	return answer;
}
double BarrierEnergy(const std::string& sequence, const std::vector<std::string>& pathway) {
	double answer = -std::numeric_limits<double>::infinity();
	for (const auto s : pathway)answer = std::max<double>(answer, EnergyOfStructure(sequence, s));
	return answer;
}
double BarrierEnergy(
	const std::string& sequence,
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list) {

	double answer = EnergyOfStructure(sequence, BasePairSetToDotNotation(A, sequence.size()));
	auto S = A;
	for (const auto action : action_list) {
		if (S.find(action) != S.end()) {
			S.erase(action);
		}
		else {
			for (const auto x : S)assert(!IsExclusive(x, action));
			S.insert(action);
		}
		answer = std::max<double>(answer, EnergyOfStructure(sequence, BasePairSetToDotNotation(S, sequence.size())));
	}
	assert(S == B);
	return answer;
}
std::pair<int, int>Difference(const std::string& structure1, const std::string& structure2) {
	//ハミング距離1の2構造を受け取って、違う塩基対1個を特定して返す。
	assert(HammingDistance(structure1, structure2) == 1);
	std::pair<int, int>answer = std::make_pair(-1, -1);
	for (int i = 0; i < structure1.size(); ++i)if (structure1[i] != structure2[i]) {
		if (answer.first == -1) {
			assert(
				(structure1[i] == '.' && structure2[i] == '(') ||
				(structure1[i] == '(' && structure2[i] == '.'));
			answer.first = i;
		}
		else if (answer.second == -1) {
			assert(
				(structure1[i] == '.' && structure2[i] == ')') ||
				(structure1[i] == ')' && structure2[i] == '.'));
			answer.second = i;
		}
		else assert(0);
	}
	return answer;
};
std::vector<std::pair<int, int>> PathwayToActionList(const std::vector<std::string>& pathway) {
	std::vector<std::pair<int, int>>answer;
	for (int i = 0; i < pathway.size() - 1; ++i) {
		answer.push_back(Difference(pathway[i], pathway[i + 1]));
	}
	return answer;
}
bool IsValidPathway(
	const std::set<std::pair<int, int>>& A,
	const std::set<std::pair<int, int>>& B,
	const std::vector<std::pair<int, int>>& action_list) {

	auto S = A;
	for (const auto action : action_list) {
		if (S.find(action) != S.end()) {
			S.erase(action);
		}
		else {
			for (const auto x : S)if (IsExclusive(x, action))return false;
			S.insert(action);
		}
	}
	return S == B;
}
bool IsValidPathway(const std::string& A, const std::string& B, const std::vector<std::string>& pathway) {
	return IsValidPathway(DotNotationToBasePairSet(A), DotNotationToBasePairSet(B), PathwayToActionList(pathway));
}

}