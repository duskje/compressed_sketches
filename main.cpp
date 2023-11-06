#include <iostream>
#include "HLLSketch.h"

int main(){
	HLLSketch hll_sketch1;
	hll_sketch1.insert_kmer((uint64_t) 17);
	hll_sketch1.insert_kmer((uint64_t) 0b010011);

	HLLSketch hll_sketch2;
	hll_sketch2.insert_kmer((uint64_t) 0b110011);
	hll_sketch2.insert_kmer((uint64_t) 0b101100);

	std::cout << hll_sketch1.get_cardinality() << std::endl;

	hll_sketch1.merge(hll_sketch2);

	std::cout << hll_sketch1.get_cardinality() << std::endl;

	return 0;
}
