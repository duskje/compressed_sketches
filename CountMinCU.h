#include <iostream>
#include <vector>
#include <string>

#include "Sketch.h"

#ifndef COUNTMINCU
#define COUNTMINCU


class CountMinCU : public Sketch {
private:
	std::vector<std::vector<uint64_t>> freqs;

	uint32_t rows;
	uint32_t cols;

	uint64_t hash(std::string key, int seed);
public:
    CountMinCU(uint32_t d, uint32_t w) : rows(d), cols(w) {
        freqs = std::vector<std::vector<uint64_t>>(rows, std::vector<uint64_t>(cols, 0));
    }

    CountMinCU(){}

    void increment_count(std::string element) override;
	uint64_t retrieve_count(std::string element) override;
    void set_count(std::string element, uint64_t count);
};

#endif
