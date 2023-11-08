#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#include "CountMinCU.h"
#include "MurmurHash3.h"


uint64_t CountMinCU::hash(std::string key, int seed){
	uint64_t out;

	char *buffer = new char[key.size() + 1];
	strcpy(buffer, key.c_str());

	MurmurHash3_x64_128(buffer, key.size(), seed, &out);

	return out % (this->cols);
}

void CountMinCU::increment_count(std::string element){
	int min_i = 0;
	int min_j = hash(element, 0);

	for(int i = 1; i < this->rows; i++){
		uint64_t h = this->hash(element, i);

		if(freqs[i][h] < freqs[min_i][min_j]){
			min_i = i;
			min_j = h;
		}
	}

	for(int i = 0; i < this->rows; i++){
		uint64_t h = this->hash(element, i);

		if(freqs[i][h] < freqs[min_i][min_j] + 1){
			freqs[i][h] = freqs[min_i][min_j] + 1;
		}
	}
}

void CountMinCU::set_count(std::string element, uint64_t count) {
	int min_i = 0;
	int min_j = hash(element, 0);

	for(int i = 1; i < (this->rows - 1); i++){
		uint64_t h = hash(element, i);

		if(freqs[i][h] < freqs[min_i][min_j]){
			min_i = i;
			min_j = h;
		}
	}

	freqs[min_i][min_j] = count;
}

uint64_t CountMinCU::retrieve_count(std::string element){
	uint64_t min_value = freqs[0][this->hash(element, 0)];

	for(int i = 1; i < this->rows; i++){
		int count = freqs[i][this->hash(element, i)];

		if(count < min_value){
			min_value = count;
		}
	}

	return min_value;
}

