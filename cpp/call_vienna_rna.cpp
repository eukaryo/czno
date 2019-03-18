/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include "call_vienna_rna.h"

#include "misc.h"


namespace czno_cpp {

//筆者の現環境(docker for windows)でのdocker containerの名前なので、適宜変えること。
#define DOCKER_CONTAINER_ID "b5ae010ba830"

static bool ExecCmd(const char* cmd, std::string& stdOut, int& exitCode) {

#ifdef __linux
	std::shared_ptr<FILE> pipe(popen(cmd, "r"), [&](FILE* p) {exitCode = pclose(p); });
#else
	std::shared_ptr<FILE> pipe(_popen(cmd, "r"), [&](FILE* p) {exitCode = _pclose(p); });
#endif
	if (!pipe) {
		return false;
	}
	std::array<char, 256> buf;
	while (!feof(pipe.get())) {
		if (fgets(buf.data(), buf.size(), pipe.get()) != nullptr) {
			stdOut += buf.data();
		}
	}
	return true;
}

static int call_count = 0;
void ResetCount() { call_count = 0; return; }
int GetCount() { return call_count; }

double EnergyOfStructure(const std::string& sequence, const std::string& structure) {
	call_count++;
	const std::string command_string =
#ifdef __linux
		std::string(R"(echo -e ')") +
		sequence +
		std::string(R"(\n)") +
		structure +
		std::string(R"(' | RNAeval)");
#else
		std::string(R"(docker exec )") +
		std::string(DOCKER_CONTAINER_ID) +
		std::string(R"( sh -c "echo ')") +
		sequence +
		std::string(R"(\n)") +
		structure +
		std::string(R"(' | RNAeval")");
#endif

	std::string stdOut;
	int exitCode;
	assert(ExecCmd(command_string.c_str(), stdOut, exitCode));
	const static std::regex value(R"(-{0,1}\d+\.\d+)");
	std::smatch m;
	std::regex_search(stdOut, m, value);
	assert(m.size() == 1);
	return std::stod(m[0].str());
}

std::vector<std::string> SampleStructure(const std::string& sequence, const int num) {
	const std::string command_string =
#ifdef __linux
		std::string(R"(echo )") +
		sequence +
		std::string(R"( | RNAsubopt -p )") +
		std::to_string(num);
#else
		std::string(R"(docker exec )") +
		std::string(DOCKER_CONTAINER_ID) +
		std::string(R"( sh -c "echo )") +
		sequence +
		std::string(R"( | RNAsubopt -p )") +
		std::to_string(num) +
		std::string(R"(")");
#endif

	std::string stdOut;
	int exitCode;
	assert(ExecCmd(command_string.c_str(), stdOut, exitCode));

	const auto lines = Split(stdOut, '\n');
	std::vector<std::string>answer;
	for (int i = 1; i < lines.size(); ++i)answer.push_back(lines[i]);

	return answer;
}

}