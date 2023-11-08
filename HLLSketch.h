#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cstring>
#include <sdsl/wavelet_trees.hpp>

#include "MurmurHash3.h"

constexpr uint64_t p = 14;
constexpr uint64_t m = 1 << p;
constexpr double alpha_A = 0.7213 / (1 + 1.079 / m);

template<typename t_wt>
class CompressedHLLSketch{
private:
	uint64_t get_zero_count(){
		uint64_t count = 0;

		for(int i = 0; i < wt.size(); i++){
			if(wt[i] == 0){
				count++;
			}
		}
		std::cout << "crash!" << std::endl;

		return count;
	}
public:
	t_wt wt;

	CompressedHLLSketch(std::array<uint64_t, m> &M){
		sdsl::construct_im(wt, M, 8);
	}

	size_t size(){
		return size_in_bytes(wt);
	}

	void merge(CompressedHLLSketch &other_sketch){
		std::vector<uint64_t> new_wt;

		for(int i = 0; i < this->wt.size(); i++){
			if(this->wt[i] < other_sketch.wt[i]){
				new_wt.push_back(other_sketch.wt[i]);
			} else {
				new_wt.push_back(this->wt[i]);
			}
		}

		sdsl::construct_im(wt, new_wt, 8);
	}

	uint64_t get_cardinality(){
		double Z = 0;

		for(int i = 0; i < m; i++){
			Z += static_cast<double>(1) / (1 << wt[i]);
		}


		double C_HLL = alpha_A * (pow(m, 2)) / Z;

		if(C_HLL <= (2.5 * m)){
			return m * log2(static_cast<double>(m) / this->get_zero_count());
		} 

		return C_HLL;
	}
};

class HLLSketch {
private:
	uint64_t get_zero_count(){
		uint64_t count = 0;

		for(auto &x: M){
			if(x == 0){
				count++;
			}
		}

		return count;
	}

	static uint8_t get_leading_zeros(uint64_t value){
		uint8_t count = 0;

		for(int offset = 64; offset; offset--){
			if( (value & (1 << offset)) ){
				return count;
			}

			count++;
		}

		return count;
	}

	uint64_t canonical_kmer(uint64_t original_kmer){
		const uint64_t complement = (~original_kmer) & 0x3fffffffffffffff;
		uint64_t reverse_complement = 0;

		for(int i = 0; i < 31; i++){
			const uint64_t source_base = (complement & (0b11 << (2 * i)));
			const uint64_t dest_base = source_base << (2 * (31 - i));

			reverse_complement |= dest_base; 
		}

		if(reverse_complement > original_kmer){
			return original_kmer;
		} else {
			return reverse_complement;
		}
	}

	uint64_t hash(uint64_t kmer_seq){
		uint8_t kmer_bytes[sizeof(kmer_seq)];

		std::memcpy(kmer_bytes, &kmer_seq, sizeof(kmer_bytes));

		uint64_t h;

 		MurmurHash3_x64_128(kmer_bytes, 4, 0, &h);

		return h;
	}

public:
	std::array<uint64_t, m> M = {};

	void insert_kmer(uint64_t kmer_seq){
		constexpr uint64_t mask_14_bits = 0xfffc000000000000;

		const uint64_t h = this->hash(kmer_seq);

		const uint64_t bucket = (h & mask_14_bits) >> (64 - p);
		const uint64_t value = (~mask_14_bits) & h;

		const uint8_t leading_zeros = this->get_leading_zeros(value);

		if(M[bucket] < leading_zeros + 1){
			M[bucket] = leading_zeros + 1;
		}
	}

	uint64_t get_cardinality(){
		double Z = 0;

		for(int i = 0; i < m; i++){
			Z += static_cast<double>(1) / (1 << M[i]);
		}

		double C_HLL = alpha_A * (pow(m, 2)) / Z;

		if(C_HLL <= (2.5 * m)){
			return m * log2(static_cast<double>(m) / this->get_zero_count());
		} 

		return C_HLL;
	}

	void merge(HLLSketch &other_sketch){
		for(int i = 0; i < this->M.size(); i++){
			if(this->M[i] < other_sketch.M[i]){
				this->M[i] = other_sketch.M[i];
			}
		}
	}

	uint64_t size_in_bytes(){
		return this->M.size() * sizeof(uint64_t);
	}

	template <typename t_wt>
	CompressedHLLSketch<t_wt> get_compressed(){
		return CompressedHLLSketch<t_wt>(M);
	}
};
