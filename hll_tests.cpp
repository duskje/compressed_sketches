#include <iostream>
#include <cassert>

#include "HLLSketch.h"

int main(){
	HLLSketch hll_sketch1;
	hll_sketch1.insert_kmer((uint64_t) 17);
	hll_sketch1.insert_kmer((uint64_t) 0b010011);

	HLLSketch hll_sketch2;
	hll_sketch2.insert_kmer((uint64_t) 0b110011);
	hll_sketch2.insert_kmer((uint64_t) 0b101100);

	auto compressed_hll_sketch1 = hll_sketch1.get_compressed<sdsl::wm_int<>>();
	auto compressed_hll_sketch2 = hll_sketch2.get_compressed<sdsl::wm_int<>>();

	assert(hll_sketch1.get_cardinality() == 2);

	hll_sketch1.merge(hll_sketch2);

	assert(hll_sketch1.get_cardinality() == 5);

	compressed_hll_sketch1.merge(compressed_hll_sketch2);

	assert(compressed_hll_sketch1.get_cardinality() == 5);

	return 0;
}
