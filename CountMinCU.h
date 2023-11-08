#include <iostream>
#include <vector>
#include <string>
#include <sdsl/wavelet_trees.hpp>

#include "Sketch.h"
#include "MurmurHash3.h"

#ifndef COUNTMINCU
#define COUNTMINCU

template<typename t_wt>
class CompressedCountMinCU {
private:
	static std::vector<uint64_t> flatten(std::vector<std::vector<uint64_t>> &freqs){
		std::vector<uint64_t> flattened;

		for(const auto &rows: freqs){
			for(const auto val: rows){
				flattened.push_back(val);
			}
		}

		return flattened;
	}


public:
	t_wt wt;

	uint32_t rows;
	uint32_t cols;

	CompressedCountMinCU(std::vector<std::vector<uint64_t>> &freqs){
		auto flattened_freqs = this->flatten(freqs);

		this->rows = freqs.size();
		this->cols = freqs[0].size();

		sdsl::construct_im(wt, flattened_freqs, 8);
	}

	uint64_t retrieve_count(std::string element) {
		uint64_t min_value = wt[this->hash(element, 0) + 1];

		for(int d = 1; d < this->rows; d++){
			uint64_t count = wt[(d * this->cols) + this->hash(element, d) + 1];

			if(count < min_value){
				min_value = count;
			}
		}

		return min_value;
	}

	size_t size_in_bytes(){
		return size_in_bytes(wt);
	}

	uint64_t hash(std::string key, int seed){
		uint64_t out;

		char *buffer = new char[key.size() + 1];
		strcpy(buffer, key.c_str());

		MurmurHash3_x64_128(buffer, key.size(), seed, &out);

		return out % (this->cols);
	}
};

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

	template<typename t_wt>
	CompressedCountMinCU<t_wt> get_compressed(){
		return CompressedCountMinCU<t_wt>(freqs);
	}
};

#endif
