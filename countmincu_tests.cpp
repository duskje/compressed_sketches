#include <cassert>
#include <iostream>
#include <sdsl/wavelet_trees.hpp>

#include "CountMinCU.h"
#include "FetchKmers.h"
#include "MurmurHash3.h"

int main(){
	auto sketch = CountMinCU(4, 1024);

	sketch.increment_count("hola");
	sketch.increment_count("hola");
	sketch.increment_count("hola");

	sketch.increment_count("the smiths");

	sketch.increment_count("kind of blue");
	sketch.increment_count("kind of blue");

	assert(sketch.retrieve_count("hola") == 3);
	assert(sketch.retrieve_count("the smiths") == 1);
	assert(sketch.retrieve_count("kind of blue") == 2);

	auto compressed_sketch_wm_int = sketch.get_compressed<sdsl::wm_int<>>();

	assert(compressed_sketch_wm_int.retrieve_count("hola") == 3);
	assert(compressed_sketch_wm_int.retrieve_count("the smiths") == 1);
	assert(compressed_sketch_wm_int.retrieve_count("kind of blue") == 2);

	auto compressed_sketch_wt_huff = sketch.get_compressed<sdsl::wt_huff<>>();

	assert(compressed_sketch_wt_huff.retrieve_count("hola") == 3);
	assert(compressed_sketch_wt_huff.retrieve_count("the smiths") == 1);
	assert(compressed_sketch_wt_huff.retrieve_count("kind of blue") == 2);
}
