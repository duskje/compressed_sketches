CC=g++
CFLAGS=-std=c++11 -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 -O0 -DNDEBUG

runtests:
	@echo "Running CountMinCU tests..."
	@$(CC) countmincu_tests.cpp CountMinCU.cpp MurmurHash3.cpp -o countmincu_tests $(CFLAGS)
	@./countmincu_tests
	@echo "Done."
	@rm countmincu_tests
	@echo "Running HyperLogLog tests..."
	@$(CC) MurmurHash3.cpp hll_tests.cpp -o hll_tests $(CFLAGS)
	@./hll_tests
	@echo "Done."
	@rm hll_tests

