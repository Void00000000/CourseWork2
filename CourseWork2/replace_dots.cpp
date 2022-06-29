#include <vector>
#include <fstream>
#include <string>
#include "replace_dots.h"

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}


void replace_dots(const std::string &fname) {
	std::ifstream fi_read;
	std::ofstream fi_write;
	fi_read.open(fname);
	fi_read.seekg(0);  // На всякий случай

	std::string line;
	std::vector<std::string> lines;
	while (std::getline(fi_read, line))
	{
		replaceAll(line, ".", ",");
		lines.push_back(line);
	}
	fi_read.close();
	fi_write.open(fname, std::ofstream::out | std::ofstream::trunc);
	for (std::string line_comma : lines) {
		fi_write << line_comma << '\n';
	}
	fi_write.close();
}