#include <iostream>
#include <deque>
#include <fstream>
#include <vector>

#include "HLLSketch.h"

void drop_first_line(std::ifstream &f){
	char curr_char;

	f.get(curr_char);

	while(curr_char != '\n'){
		f.get(curr_char);
	}
}

uint64_t encode_kmer(std::deque<char> kmer){
	uint64_t encoded = 0;

	for(uint64_t k = 0; k < kmer.size(); k++){
		uint8_t base_2_bits;

		if(kmer[k] == 'A'){
			base_2_bits = 0b00;
		} else if(kmer[k] == 'C'){
			base_2_bits = 0b01;
		} else if(kmer[k] == 'G'){
			base_2_bits = 0b10;
		} else if(kmer[k] == 'T'){
			base_2_bits = 0b11;
		} 

		encoded |= base_2_bits << (2 * k);
	}

	return encoded;
}

std::vector<uint64_t> fetch_kmers(std::ifstream &f){
	std::vector<uint64_t> kmers;

	std::deque<char> kmer;
	char curr_char;

	for(int k = 0; k < 31; k++){
		if(f.eof()){
			break;
		}

		if(curr_char == '>'){
			break;
		}

		if(curr_char == '\n'){
			continue;
		}

		kmer.push_back(curr_char);

		f.get(curr_char);
	}
	
	kmers.push_back(encode_kmer(kmer));

	while(!f.eof() && curr_char != '>'){
	 	kmer.pop_front();
		kmer.push_back(curr_char);

		kmers.push_back(encode_kmer(kmer));

		f.get(curr_char);
	}

	return kmers;
}

int main(){
	std::ifstream f("GCF_000717965.1_ASM71796v1_genomic.fna", std::ios_base::in);

	drop_first_line(f);

	HLLSketch hll_sketch;

	for(auto &kmer: fetch_kmers(f)){
		hll_sketch.insert_kmer(kmer);
	}

	auto compressed_hll_wm_int = hll_sketch.get_compressed<sdsl::wm_int<>>();

	std::cout << "Compression ratio wm_int: " << static_cast<double>(compressed_hll_wm_int.size_in_bytes()) / hll_sketch.size() << std::endl;

	auto compressed_hll_wt_huff = hll_sketch.get_compressed<sdsl::wt_huff<>>();

	std::cout << "Compression ratio wt_huff: " << static_cast<double>(compressed_hll_wt_huff.size_in_bytes()) / hll_sketch.size() << std::endl;

	return 0;
}
