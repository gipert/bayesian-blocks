CFLAGS = -Wall -Wextra -Werror \
         -std=c++11 -O3 -L. -I.. \
         -fprofile-arcs -ftest-coverage \
         $(shell root-config --cflags)

LIBS = $(shell root-config --libs)

run_test : run_test.cc ../include/bayesian_blocks.hpp ../include/bayesian_blocks_root.hpp
	$(CXX) $(CFLAGS) -o $@ $< $(LIBS)

test : run_test
	./run_test

clean :
	rm -f *.so *.gcov *.gcda *.gcno run_test

.PHONY : test clean
