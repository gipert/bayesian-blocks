CFLAGS = -Wall -Wextra -Werror \
         -std=c++11 -O3 -L. -I.. \
         -fprofile-arcs -ftest-coverage \
         $(shell root-config --cflags)

LIBS = $(shell root-config --libs)

run_test : run_test.cc libbayesian_blocks.so
	$(CXX) $(CFLAGS) -o $@ $< $(LIBS) -lbayesian_blocks

libbayesian_blocks.so : ../bayesian_blocks.cc ../bayesian_blocks.hh
	$(CXX) $(CFLAGS) -fPIC -shared -o $@ $< $(LIBS)

test : run_test
	./run_test

clean :
	rm -f *.so *.gcov *.gcda *.gcno run_test

.PHONY : test clean
